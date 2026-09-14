#!/usr/bin/env python3
"""
F2 feature assembly — per target junction: cis splice-site strength + eCLIP RBP binding proximity.

The mechanistic grammar features: a junction's usage should depend on (a) how strong its
splice sites are (cis), and (b) which RBPs bind near it (cis binding) INTERACTED with those RBPs'
cell-type expression (trans). This script builds the junction-level (cell-type-invariant) part:
  - pwm_donor  : 5'SS PWM log-odds score at donor site a
  - pwm_acceptor: 3'SS PWM log-odds score at acceptor site b
  - eCLIP binding: for each RBP with a peak within +/-WIN bp of a or b -> (junction_id, rbp) pair
The (junction x RBP binding) x (cell-type RBP expression) interaction is formed in the probe.

Main chromosomes only (chr1-22,X,Y; scaffolds excluded — tiny fraction, naming mismatch with eCLIP).
Outputs (results/paper1/f2_grammar/features/):
  f2_junction_cis.parquet   : junction_id, chrom, strand, a, b, pwm_donor, pwm_acceptor
  f2_junction_rbp.parquet   : junction_id, rbp  (long form; RBP binds within +/-WIN of a or b)
Usage: f2_features.py [--win 50]
"""
import sys, os, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2_splice_pwm as sp

WIN = int(sys.argv[sys.argv.index("--win") + 1]) if "--win" in sys.argv else 50
RESULTS_ROOT = os.environ.get("SCTHREAD_RESULTS_ROOT", "results")
FEAT = os.environ.get("SCTHREAD_F2_FEATURES", f"{RESULTS_ROOT}/paper1/f2_grammar/features")
ECLIP = os.environ.get("SCTHREAD_ECLIP_FEATURES", f"{RESULTS_ROOT}/paper1/f2_grammar/rbp_peaks_union.parquet")
PWM = json.load(open(os.environ.get("SCTHREAD_SPLICE_PWM", f"{RESULTS_ROOT}/paper1/f2_grammar/splice_pwm.json")))
don_freq = np.array(PWM["donor"]); acc_freq = np.array(PWM["acceptor"])
FA = os.environ.get("SCTHREAD_GENOME_FASTA", "reference/genome.fa")

MAIN = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
tgt = pd.read_parquet(f"{FEAT}/f2_target_junctions.parquet")
tgt = tgt[tgt.chrom.isin(MAIN)].copy()
print(f"target junctions on main chroms: {len(tgt):,}")

# ---- cis: PWM scores (vectorized loop; genome.fa via pyfaidx, junction chrom matches) ----
from pyfaidx import Fasta
fa = Fasta(FA, rebuild=False)
def sc(seq, freq):
    if len(seq) != freq.shape[0] or any(b not in "ACGT" for b in seq):
        return np.nan
    return float(sum(np.log2(freq[i, sp.BIDX[b]] / 0.25) for i, b in enumerate(seq)))
chrom = tgt.chrom.values; strand = tgt.strand.values
av = tgt.a.values.astype(int); bv = tgt.b.values.astype(int)
pd_don = np.full(len(tgt), np.nan); pd_acc = np.full(len(tgt), np.nan)
for i in range(len(tgt)):
    try:
        pd_don[i] = sc(sp.donor_seq(fa, chrom[i], av[i], strand[i]), don_freq)
        pd_acc[i] = sc(sp.acceptor_seq(fa, chrom[i], bv[i], strand[i]), acc_freq)
    except Exception:
        pass
tgt["pwm_donor"] = pd_don; tgt["pwm_acceptor"] = pd_acc
tgt[["junction_id", "chrom", "strand", "a", "b", "pwm_donor", "pwm_acceptor"]].to_parquet(
    f"{FEAT}/f2_junction_cis.parquet", index=False)
print(f"PWM: donor mean={np.nanmean(pd_don):.2f} acceptor mean={np.nanmean(pd_acc):.2f} "
      f"| NaN donor={np.isnan(pd_don).mean():.1%}")

# ---- eCLIP RBP binding proximity via duckdb range join ----
import duckdb
tmp = f"/tmp/f2feat_{os.getpid()}"; os.makedirs(tmp, exist_ok=True)
con = duckdb.connect()
con.execute(f"PRAGMA threads=4; PRAGMA temp_directory='{tmp}'; PRAGMA memory_limit='40GB'")
# pre-filter eCLIP to main chroms + real RBPs once (big reduction before range join)
con.execute(f"""CREATE TEMP TABLE ec AS
  SELECT chrom, pstart, pend, rbp FROM read_parquet('{ECLIP}')
  WHERE rbp <> '.' AND chrom IN ({','.join(chr(39)+c+chr(39) for c in MAIN)})""")
con.execute("CREATE INDEX ec_idx ON ec(chrom)")
# chunk junctions (bounds the range-join intermediate); two clean overlap joins per chunk
jt_all = tgt[["junction_id", "chrom", "a", "b"]].reset_index(drop=True)
CH = 20000
parts = []
for s0 in range(0, len(jt_all), CH):
    chunk = jt_all.iloc[s0:s0 + CH]
    con.register("jt", chunk)
    part = con.execute(f"""
    WITH sites AS (
      SELECT junction_id, chrom, a-{WIN} AS lo, a+{WIN} AS hi FROM jt
      UNION ALL
      SELECT junction_id, chrom, b-{WIN} AS lo, b+{WIN} AS hi FROM jt
    )
    SELECT DISTINCT s.junction_id, e.rbp
    FROM sites s JOIN ec e ON s.chrom = e.chrom
    WHERE e.pend >= s.lo AND e.pstart <= s.hi
    """).fetchdf()
    parts.append(part)
    con.unregister("jt")
    print(f"  eCLIP chunk {s0//CH + 1}/{-(-len(jt_all)//CH)}: {len(part):,} pairs", flush=True)
rbp = pd.concat(parts, ignore_index=True).drop_duplicates()
con.close(); os.system(f"rm -rf {tmp}")
rbp.to_parquet(f"{FEAT}/f2_junction_rbp.parquet", index=False)
print(f"eCLIP binding pairs: {len(rbp):,} | junctions with >=1 RBP: {rbp.junction_id.nunique():,}/{len(tgt):,} "
      f"| distinct RBP: {rbp.rbp.nunique()}")
print(f"RBP/junction: mean={rbp.groupby('junction_id').size().mean():.1f}")
