#!/usr/bin/env python3
"""
DIU (analysis #3) — per-run isoform expression pseudobulk by cell type.

Loads the IsoQuant transcript x CB count matrix, joins CB -> release_v1 cell_type, aggregates to
per (cell_type, transcript) counts, maps ENST -> ENSG, and emits per (gene, transcript, cell_type)
pseudobulk isoform counts. Downstream: per-gene isoform usage per cell type + donor-aware DIU
(differential isoform usage) for multi-isoform genes.

Output: results/paper1/f2_grammar/agg_diu/<SRR>.diu.parquet
        columns: gene_id, transcript_id, cell_type, count
Usage: diu_reduce.py <SRR>
"""
import sys, os, gzip
from pathlib import Path
import numpy as np, pandas as pd
import scipy.io as sio, scipy.sparse as sp
from species_resources import ANNOTATION, annotation_pattern, isoquant_matrix_files, run_config

RUN = sys.argv[1]
REGISTRY_ROW, RESOURCE = run_config(RUN)
SPECIES = REGISTRY_ROW["species"].strip().lower()
GSE = REGISTRY_ROW["gse"]
GTF = Path(RESOURCE["gtf"])
RESULTS_ROOT = Path(os.environ.get("SCTHREAD_RESULTS_ROOT", "results"))
T2G_CACHE = Path(
    RESULTS_ROOT
    / Path("paper1/f2_grammar")
    / f"transcript2gene.{SPECIES}.{RESOURCE['assembly']}.parquet"
)
OUT = os.environ.get("SCTHREAD_DIU_OUTPUT", "outputs/f2_grammar/agg_diu")
os.makedirs(OUT, exist_ok=True)
MATRIX_FILES = isoquant_matrix_files(RUN, REGISTRY_ROW)
MTX = MATRIX_FILES["matrix"]
FEAT = MATRIX_FILES["features"]
BC = MATRIX_FILES["barcodes"]

# Assembly-matched reference transcript -> gene map.
if os.path.exists(T2G_CACHE):
    t2g = pd.read_parquet(T2G_CACHE)
else:
    rows = []
    import re
    opener = gzip.open if str(GTF).endswith(".gz") else open
    for line in opener(GTF, "rt"):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "transcript":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]); g = re.search(r'gene_id "([^"]+)"', f[8])
        if t and g:
            rows.append((t.group(1).split(".")[0], g.group(1)))
    t2g = pd.DataFrame(rows, columns=["transcript_id", "gene_id"]).drop_duplicates(
        "transcript_id"
    )
    t2g.to_parquet(T2G_CACHE, index=False)
t2g_map = dict(zip(t2g.transcript_id, t2g.gene_id))

# load matrix (transcript x CB)
m = sio.mmread(MTX).tocsr()   # features x barcodes
feats = [l.split("\t")[0].split(".")[0] for l in open(FEAT).read().splitlines()]
bcs = [l.split("\t")[0].split("-")[0] for l in open(BC).read().splitlines()]  # 16bp core
_pattern = annotation_pattern(RUN, REGISTRY_ROW)
_sub = _pattern.strip("%")
lab = pd.read_parquet(ANNOTATION, columns=["gse", "source_cell_id", "barcode_norm", "active_label"])
lab = lab[lab.gse == GSE]
if _sub:
    lab = lab[lab.source_cell_id.str.contains(_sub, regex=False)]
# defensive: some releases carry a degenerate barcode_norm (e.g. GSE295353 all "1"); the true
# 16bp core lives in the computational key source_cell_id = gse::run::BARCODE_suffix. Recover it
# and restrict to this run so barcodes are collision-safe. No-op when barcode_norm is normal.
_bn = lab.barcode_norm.astype(str)
if len(_bn) and (_bn.str.len() <= 3).mean() > 0.5:
    _scid = lab.source_cell_id.astype(str)
    lab = lab[_scid.str.contains(f"::{RUN}::", regex=False)].copy()
    lab["barcode_norm"] = (
        lab.source_cell_id.astype(str).str.split("::").str[-1]
        .str.split("_").str[0].str.split("-").str[0]
    )
bc2ct = dict(zip(lab.barcode_norm, lab.active_label))
cts = np.array([bc2ct.get(b, "Unlabeled") for b in bcs])

rows = []
for ct in pd.unique(cts):
    if ct in {"Unlabeled", "Unassigned"}:
        continue
    colmask = np.where(cts == ct)[0]
    if len(colmask) < 5:
        continue
    v = np.asarray(m[:, colmask].sum(axis=1)).ravel()   # per-transcript total in this cell type
    nz = np.where(v > 0)[0]
    for ti in nz:
        g = t2g_map.get(feats[ti])
        if g:
            rows.append((g, feats[ti], ct, float(v[ti])))
df = pd.DataFrame(rows, columns=["gene_id", "transcript_id", "cell_type", "count"])
df["run"] = RUN; df["gse"] = GSE
if df.empty:
    sys.exit(f"NO_LABELED_DIU_ROWS {RUN} {GSE}")
if not df.gene_id.str.startswith(str(RESOURCE["gene_prefix"])).all():
    sys.exit(f"GENE_ASSEMBLY_MISMATCH {RUN} {SPECIES}")
target = Path(OUT) / f"{RUN}.diu.parquet"
temporary = target.with_suffix(".parquet.partial")
df.to_parquet(temporary, index=False)
temporary.replace(target)
# feasibility: multi-isoform genes with >=2 cell types
multi = df.groupby(["gene_id", "cell_type"]).transcript_id.nunique().reset_index()
multigenes = df.groupby("gene_id").transcript_id.nunique()
print(f"[{RUN}] {GSE} rows={len(df):,} cell_types={df.cell_type.nunique()} "
      f"| genes={df.gene_id.nunique():,} | multi-isoform genes={ (multigenes>=2).sum():,}")
