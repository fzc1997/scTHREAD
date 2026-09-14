#!/usr/bin/env python3
"""
APA (analysis #2) — per-run 3'-end → PolyASite assignment + per (gene, cell_type, PAS) counts.

For each admitted read: assign observed_3p_raw to a PolyASite 2.0 cluster within +/-WIN bp on the
same chromosome+strand (PolyASite chrom lacks 'chr' prefix -> add it). Join CB -> cell_type.
Emit per (gene, cell_type, pas_id) molecule counts, plus an endpoint_evidence_tier_3p breakdown
for the internal-priming / unsupported-end negative control (PAS-supported vs not).

Output: results/paper1/f2_grammar/agg_apa/<SRR>.apa.parquet
        columns: gene_id, cell_type, pas_id, pas_strand, n_molecules, n_reads
Usage: apa_reduce.py <SRR> [LIMIT]
"""
import sys, os
from pathlib import Path
from species_resources import ANNOTATION, annotation_pattern, resolved_read_facts, run_config

RUN = sys.argv[1]
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 0
WIN = 25
REGISTRY_ROW, RESOURCE = run_config(RUN)
SPECIES = REGISTRY_ROW["species"].strip().lower()
GSE = REGISTRY_ROW["gse"]
RRF = resolved_read_facts(RUN)
PAS = Path(RESOURCE["polyasite"])
PAS_COLUMN_COUNT = 11 if SPECIES == "human" else 6
PAS_COLUMNS_SQL = "{" + ",".join(
    f"'column{i}':'VARCHAR'" for i in range(PAS_COLUMN_COUNT)
) + "}"
OUT = os.environ.get("SCTHREAD_APA_OUTPUT", "outputs/f2_grammar/agg_apa")
os.makedirs(OUT, exist_ok=True)
ANN_PAT = annotation_pattern(RUN, REGISTRY_ROW)

# DuckDB spill dir. MUST be a REAL DISK: if /tmp is RAM-backed (tmpfs) on your machine,
# set REDUCER_TMP to a real disk path or spilling silently eats RAM and fills up ("No space left").
tmp = f"{os.environ.get('REDUCER_TMP', '/tmp')}/apa_{RUN}_{os.getpid()}"
os.makedirs(tmp, exist_ok=True)
import duckdb
con = duckdb.connect()
_MEM = os.environ.get("APA_DUCKDB_MEM", "40GB")  # raise on RAM-rich machines to spill less
con.execute(f"PRAGMA threads=6; PRAGMA temp_directory='{tmp}'; PRAGMA memory_limit='{_MEM}'")
# PolyASite tables differ in width between species. Only BED columns 1–6 are
# required; chromosome names are normalized without changing coordinates.
con.execute(f"""
CREATE TEMP TABLE pas AS
SELECT CASE WHEN starts_with(column0, 'chr') THEN column0 ELSE 'chr'||column0 END AS pchrom,
       CAST(column1 AS BIGINT) AS pstart, CAST(column2 AS BIGINT) AS pend,
       column3 AS pas_id, column5 AS pstrand
FROM read_csv(
  '{PAS}', delim='\t', header=false, compression='auto',
  columns={PAS_COLUMNS_SQL}
)
""")
con.execute(f"""
CREATE TEMP TABLE ann AS
SELECT barcode_norm AS bc, active_label AS ct
FROM read_parquet('{ANNOTATION}') WHERE gse='{GSE}' AND source_cell_id LIKE '{ANN_PAT}'
""")
src = f"(SELECT * FROM read_parquet('{RRF}') LIMIT {LIMIT})" if LIMIT else f"read_parquet('{RRF}')"
con.execute(f"""
CREATE TEMP VIEW rf AS
SELECT r.chromosome AS chrom, r.transcript_strand AS strand, r.gene_id,
       r.observed_3p_raw AS p3, r.cell_barcode||'|'||r.umi AS mol,
       COALESCE(a.ct,'Unlabeled') AS cell_type,
       r.endpoint_evidence_tier_3p AS etier
FROM {src} r
LEFT JOIN ann a ON r.cell_barcode=a.bc
WHERE r.read_admitted_primary AND r.gene_id IS NOT NULL AND r.observed_3p_raw IS NOT NULL
""")
# range join: 3' end within +/-WIN of a PAS cluster, same chrom+strand
df = con.execute(f"""
SELECT rf.gene_id, rf.cell_type, pas.pas_id, pas.pstrand,
       count(DISTINCT rf.mol) AS n_molecules, count(*) AS n_reads
FROM rf JOIN pas
  ON rf.chrom = pas.pchrom AND rf.strand = pas.pstrand
     AND rf.p3 BETWEEN pas.pstart - {WIN} AND pas.pend + {WIN}
GROUP BY 1,2,3,4
""").fetchdf()
# also: PAS-supported vs unsupported (negative control denominator)
support = con.execute(f"""
SELECT count(DISTINCT mol) AS total_mol,
       count(DISTINCT mol) FILTER (WHERE EXISTS(
         SELECT 1 FROM pas WHERE rf.chrom=pas.pchrom AND rf.strand=pas.pstrand
           AND rf.p3 BETWEEN pas.pstart-{WIN} AND pas.pend+{WIN})) AS pas_supported_mol
FROM rf
""").fetchone()
df["run"] = RUN; df["gse"] = GSE
if not LIMIT:
    if df.empty:
        raise SystemExit(f"NO_PAS_ROWS {RUN} {GSE}")
    if not df.gene_id.str.startswith(str(RESOURCE["gene_prefix"])).all():
        raise SystemExit(f"GENE_ASSEMBLY_MISMATCH {RUN} {SPECIES}")
    target = Path(OUT) / f"{RUN}.apa.parquet"
    # per-writer partial name so concurrent double-runners never corrupt each other's
    # partial; both atomically replace the same target -> whoever finishes writes it (identical content).
    temporary = target.with_suffix(f".parquet.{os.getpid()}.partial")
    df.to_parquet(temporary, index=False)
    temporary.replace(target)
con.close(); os.system(f"rm -rf {tmp}")
lab = df[df.cell_type != "Unlabeled"]
print(f"[{RUN}] {GSE} PAS-assigned rows={len(df):,} | labeled={len(lab):,} "
      f"| PAS support rate={support[1]/support[0]:.1%} ({support[1]:,}/{support[0]:,} mol) "
      f"| genes with >=2 PAS={df[df.cell_type!='Unlabeled'].groupby('gene_id').pas_id.nunique().ge(2).sum():,}")
