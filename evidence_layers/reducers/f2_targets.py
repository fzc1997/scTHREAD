#!/usr/bin/env python3
"""
F2 target-junction definition + response matrix (the predictive-grammar substrate).

From the per-run junction x cell_type reduction (agg_jct/<SRR>.jct_ct.parquet), define high-confidence
RECURRENT NOVEL junctions and their donor(5')-anchored competition PSI response per (donor, cell_type).

Target junction criteria:
  - NOVEL: never observed in an FSM read anywhere (sum n_fsm == 0 across all runs)
  - recurrent: total molecules >= MIN_MOL (default 5) AND seen in >= MIN_STUDIES (default 2) studies
  - competing splice choice: its 5' donor site (chrom,strand,a) has >= 2 distinct acceptors b
    (so junction usage is a real CHOICE, not the only option)
  - canonical GT-AG at the intron boundaries (added as is_canonical; targets filtered to canonical)

Response (PSI) = donor(a)-anchored competition usage:
  PSI(junction j, donor d, cell_type c) = molecules(j | d,c) / molecules(all junctions at 5' site a_j | d,c)
Requires the site to have >= MIN_SITE_MOL molecules in that (donor,cell_type) to be defined.

Outputs (outputs/f2_grammar/features/):
  f2_target_junctions.parquet : junction_id, chrom, strand, a, b, n_mol, n_studies, n_donors, is_canonical
  f2_response_psi.parquet     : junction_id, gse, donor, cell_type, platform, tgt_mol, site_mol, psi
Usage: f2_targets.py [--min-mol N] [--min-studies N] [--min-site-mol N] [--no-canonical]
"""
import sys, glob, os, json
import numpy as np, pandas as pd

MIN_MOL = int(sys.argv[sys.argv.index("--min-mol") + 1]) if "--min-mol" in sys.argv else 5
MIN_STUDIES = int(sys.argv[sys.argv.index("--min-studies") + 1]) if "--min-studies" in sys.argv else 2
MIN_SITE_MOL = int(sys.argv[sys.argv.index("--min-site-mol") + 1]) if "--min-site-mol" in sys.argv else 10
DO_CANON = "--no-canonical" not in sys.argv

FA = os.environ.get("SCTHREAD_GENOME_FASTA", "reference/genome.fa")  # pyfaidx-indexed
LED = os.environ.get("SCTHREAD_INPUT_LEDGER", "metadata/science_input_ledger_v2.json")
RESULTS_ROOT = os.environ.get("SCTHREAD_RESULTS_ROOT", "results")
OUT = os.environ.get("SCTHREAD_F2_FEATURES", f"{RESULTS_ROOT}/paper1/f2_grammar/features")
os.makedirs(OUT, exist_ok=True)
recs = json.load(open(LED))["records"]
# platform harmonization: OXFORD_NANOPORE and ONT are the same platform (all data is Nanopore; NO PacBio
# anywhere in the cohort -> leave-one-platform-out is not possible, only leave-one-study-out).
_PLATMAP = {"OXFORD_NANOPORE": "ONT", "ONT": "ONT"}
r2plat = {r["run_accession"]: _PLATMAP.get(r.get("platform"), r.get("platform", "NA")) for r in recs}

files = glob.glob(os.environ.get("SCTHREAD_JUNCTION_INPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/agg_jct/*.jct_ct.parquet"))
if not files:
    sys.exit("no agg_jct output yet")
print(f"loading {len(files)} run reductions")
d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
d = d[d.cell_type != "Unlabeled"].copy()
d["donor"] = d.run
d["platform"] = d.run.map(r2plat).fillna("NA")
d["jid"] = d.chrom.astype(str) + ":" + d.strand + ":" + d.a.astype(str) + "-" + d.b.astype(str)

# junction-level recurrence / novelty
jl = d.groupby(["jid", "chrom", "strand", "a", "b"], as_index=False).agg(
    n_mol=("n_molecules", "sum"), n_fsm=("n_fsm", "sum"), n_novelcls=("n_novelcls", "sum"),
    n_studies=("gse", "nunique"), n_donors=("donor", "nunique"))
jl["is_novel"] = jl.n_fsm == 0
# competing: 5' donor site (chrom,strand,a) has >=2 distinct acceptors
site = jl.groupby(["chrom", "strand", "a"]).b.nunique().reset_index(name="n_acc")
jl = jl.merge(site, on=["chrom", "strand", "a"])
jl["competing"] = jl.n_acc >= 2
tgt = jl[(jl.is_novel) & (jl.n_mol >= MIN_MOL) & (jl.n_studies >= MIN_STUDIES) & (jl.competing)].copy()
print(f"junctions total={len(jl):,} | novel={jl.is_novel.sum():,} | "
      f"novel+recurrent+competing targets(pre-canonical)={len(tgt):,}")

# canonical GT-AG — EXACT convention from gtag_junctions.py (validated pipeline):
#   donor = fa[chrom][a:a+2], acceptor = fa[chrom][b-3:b-1]; minus strand -> swap + reverse-complement
if DO_CANON and len(tgt):
    from pyfaidx import Fasta
    fa = Fasta(FA, rebuild=False)
    _CMP = str.maketrans("ACGTN", "TGCAN")
    def rc(x): return x.translate(_CMP)[::-1]
    # vectorized: raw numpy arrays + tight loop (avoids per-row Series overhead of df.apply)
    chrom = tgt.chrom.values; strand = tgt.strand.values
    av = tgt.a.values.astype(int); bv = tgt.b.values.astype(int)
    out = np.zeros(len(tgt), dtype=bool)
    for i in range(len(tgt)):
        try:
            donor = str(fa[chrom[i]][av[i]:av[i] + 2]).upper()
            acc = str(fa[chrom[i]][bv[i] - 3:bv[i] - 1]).upper()
        except Exception:
            continue
        if strand[i] == "-":
            donor, acc = rc(acc), rc(donor)
        out[i] = (donor == "GT" and acc == "AG")
    tgt["is_canonical"] = out
    print(f"canonical GT-AG among targets: {tgt.is_canonical.sum():,}/{len(tgt):,} "
          f"({tgt.is_canonical.mean():.1%})")
    tgt = tgt[tgt.is_canonical].copy()
else:
    tgt["is_canonical"] = np.nan

tgt[["jid", "chrom", "strand", "a", "b", "n_mol", "n_studies", "n_donors", "is_canonical"]].rename(
    columns={"jid": "junction_id"}).to_parquet(f"{OUT}/f2_target_junctions.parquet", index=False)

# response PSI: per (jid, donor, cell_type) tgt_mol / site_mol
tgt_keys = set(tgt.jid)
tgt_sites = tgt[["chrom", "strand", "a"]].drop_duplicates()
dd = d.merge(tgt_sites, on=["chrom", "strand", "a"])           # all junctions at target donor sites
site_mol = dd.groupby(["chrom", "strand", "a", "gse", "donor", "cell_type", "platform"], as_index=False).n_molecules.sum().rename(columns={"n_molecules": "site_mol"})
tgt_mol = dd[dd.jid.isin(tgt_keys)].groupby(["jid", "chrom", "strand", "a", "gse", "donor", "cell_type", "platform"], as_index=False).n_molecules.sum().rename(columns={"n_molecules": "tgt_mol"})
resp = tgt_mol.merge(site_mol, on=["chrom", "strand", "a", "gse", "donor", "cell_type", "platform"])
resp = resp[resp.site_mol >= MIN_SITE_MOL].copy()
resp["psi"] = resp.tgt_mol / resp.site_mol
resp = resp.rename(columns={"jid": "junction_id"})
resp[["junction_id", "gse", "donor", "cell_type", "platform", "tgt_mol", "site_mol", "psi"]].to_parquet(
    f"{OUT}/f2_response_psi.parquet", index=False)
print(f"response rows={len(resp):,} | target junctions w/ response={resp.junction_id.nunique():,} | "
      f"donors={resp.donor.nunique()} | studies={resp.gse.nunique()} | platforms={sorted(resp.platform.unique())}")
print(f"PSI dist: mean={resp.psi.mean():.3f} median={resp.psi.median():.3f} "
      f"frac>0.1={np.mean(resp.psi>0.1):.1%} frac>0.5={np.mean(resp.psi>0.5):.1%}")
