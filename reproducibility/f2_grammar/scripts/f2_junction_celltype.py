#!/usr/bin/env python3
"""
F2.0 — per-run junction × cell_type reduction (atomic table for the predictive grammar).

For each distinct splice junction (chrom,strand,donor_a,acceptor_b) in a run, count molecules
(CB+UMI) per cell_type, plus how many of its reads are FSM vs novel-class (NIC/NNC). This lets
downstream:
  - classify a junction as KNOWN (n_fsm>0) vs NOVEL (n_fsm==0, seen only in NIC/NNC reads)
  - build donor(a)-anchored competition PSI = novel-junction molecules / all-junction molecules at donor a
  - all resolved per (cell_type, donor=run)

Cell-type labels join collision-safe on (run, barcode) via release_v1 source_cell_id.
Output: results/paper1/f2_grammar/agg_jct/<SRR>.jct_ct.parquet
        columns: chrom,strand,a,b,cell_type,n_molecules,n_reads,n_fsm,n_novelcls
Usage: f2_junction_celltype.py <SRR>
"""
import sys, os, shutil
from pathlib import Path
from species_resources import ANNOTATION, annotation_pattern, resolved_read_facts, run_config

RUN = sys.argv[1]
REGISTRY_ROW, RESOURCE = run_config(RUN)
SPECIES = REGISTRY_ROW["species"].strip().lower()
GSE = REGISTRY_ROW["gse"]
RRF = resolved_read_facts(RUN)
OUT = os.environ.get("SCTHREAD_JUNCTION_OUTPUT", "results/paper1/f2_grammar/agg_jct")
os.makedirs(OUT, exist_ok=True)
ANN_PAT = annotation_pattern(RUN, REGISTRY_ROW)

# DuckDB spill dir. MUST be a REAL DISK: HPC /tmp is XFS (ok default), but Tower /tmp is tmpfs=RAM ->
# set REDUCER_TMP=/mnt/extssd2 on Towers or spilling silently eats RAM and fills up ("No space left").
tmp = f"{os.environ.get('REDUCER_TMP', '/tmp')}/f2jct_{RUN}_{os.getpid()}"
os.makedirs(tmp, exist_ok=True)
import duckdb
con = duckdb.connect()
_MEM = os.environ.get("F2_DUCKDB_MEM", "24GB")  # Tower has huge RAM -> raise to spill less
con.execute(f"PRAGMA threads=4; PRAGMA temp_directory='{tmp}'; PRAGMA memory_limit='{_MEM}'")
# per-run label table (collision-safe: filter to this run)
con.execute(f"""
CREATE TEMP TABLE ann AS
SELECT barcode_norm AS bc, active_label AS ct
FROM read_parquet('{ANNOTATION}')
WHERE gse='{GSE}' AND source_cell_id LIKE '{ANN_PAT}'
""")
# unnest chain -> per (junction, cell_type): distinct molecules + reads + class breakdown
df = con.execute(f"""
WITH J AS (
  SELECT r.chromosome AS chrom, r.transcript_strand AS strand,
         CAST(split_part(t.j,'-',1) AS BIGINT) AS a,
         CAST(split_part(t.j,'-',2) AS BIGINT) AS b,
         r.classification_status AS cls,
         r.cell_barcode || '|' || r.umi AS mol,
         COALESCE(ann.ct,'Unlabeled') AS cell_type
  FROM read_parquet('{RRF}') r
  LEFT JOIN ann ON r.cell_barcode = ann.bc,
  UNNEST(r.splice_chain_transcript_order) AS t(j)
  WHERE r.read_admitted_primary AND r.exon_class='MULTI_EXON' AND t.j LIKE '%-%'
)
SELECT chrom, strand, a, b, cell_type,
       count(DISTINCT mol) AS n_molecules,
       count(*) AS n_reads,
       count(*) FILTER (WHERE cls='FSM') AS n_fsm,
       count(*) FILTER (WHERE cls IN ('NIC','NNC')) AS n_novelcls
FROM J GROUP BY 1,2,3,4,5
""").fetchdf()
df["run"] = RUN; df["gse"] = GSE
if df.empty:
    raise SystemExit(f"NO_JUNCTION_ROWS {RUN} {GSE}")
target = Path(OUT) / f"{RUN}.jct_ct.parquet"
# per-writer partial name so concurrent double-runners (Tower1/Tower2/HPC) never corrupt each other's
# partial; both atomically replace the same target -> whoever finishes writes it (identical content).
temporary = target.with_suffix(f".parquet.{os.getpid()}.partial")
df.to_parquet(temporary, index=False)
temporary.replace(target)
con.close()
shutil.rmtree(tmp, ignore_errors=True)

lab = df[df.cell_type != "Unlabeled"]
print(f"[{RUN}] {GSE} junction×ct rows={len(df):,} | labeled rows={len(lab):,} "
      f"| distinct junctions={df[['chrom','a','b']].drop_duplicates().shape[0]:,} "
      f"| cell types={lab.cell_type.nunique()}")
