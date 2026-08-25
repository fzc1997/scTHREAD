#!/usr/bin/env python3
"""
APA #2 — cell-type differential poly(A)-site usage (donor-aware, permutation-calibrated).

Structurally identical to DIU (usage vectors over a categorical set): here the set is PolyASite 2.0
clusters (pas_id) and the weight is molecule count. For each multi-PAS gene: per (cell_type, donor)
build the PAS-usage vector (fraction of the gene's molecules at each PAS). Statistic = donor-adjusted
between-cell-type divergence of PAS usage (weighted L2 of donor-centered usage). Within-donor cell_type
permutation -> empirical p. Effect = max pairwise cell-type total-variation distance (donor-pooled).
Gate q<0.05 & TV>=0.20.

Internal-priming note: PAS are PolyASite 2.0 annotated clusters (already filtered for internal priming),
NOT de novo 3'-end peaks, so A-rich internal-priming artifacts are controlled at the annotation level.

Marrow studies (GSE307660, GSE276974), harmonized labels. Min gene depth per (ct,donor) >=20 molecules;
>=2 PAS; >=2 cell types; >=3 donors. Exclude 'Unassigned'.
Usage: apa_test.py [K] [--unit N] [--negctrl]
"""
import sys, glob, os
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ct_harmonize import harmonize

K = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1000
UNIT = int(sys.argv[sys.argv.index("--unit") + 1]) if "--unit" in sys.argv else None
NEGCTRL = "--negctrl" in sys.argv
RESULTS_ROOT = os.environ.get("SCTHREAD_RESULTS_ROOT", "results")
INPUT_GLOB = os.environ.get("SCTHREAD_APA_INPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/agg_apa/*.apa.parquet")
OUTPUT = os.environ.get("SCTHREAD_APA_TEST_OUTPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/figdata/apa_celltype.tsv")

d = pd.concat([pd.read_parquet(f) for f in glob.glob(
    INPUT_GLOB)])
d = d[d.gse.isin(["GSE307660", "GSE276974"])].copy()
d["cell_type"] = harmonize(d.cell_type)
d = d[d.cell_type != "Unassigned"].copy()
d = d.rename(columns={"n_molecules": "count"})
d = d.groupby(["gene_id", "pas_id", "cell_type", "run", "gse"], as_index=False)["count"].sum()
d["donor"] = d.run
if NEGCTRL:
    rng0 = np.random.RandomState(7)
    # shuffle cell_type within donor across the (gene,pas) rows -> destroy cell-type structure
    d["cell_type"] = d.groupby("donor")["cell_type"].transform(lambda x: rng0.permutation(x.values))
gd = d.groupby(["gene_id", "cell_type", "donor"])["count"].sum().reset_index(name="gdepth")
d = d.merge(gd, on=["gene_id", "cell_type", "donor"])
d = d[d.gdepth >= 20]
ni = d.groupby("gene_id").pas_id.nunique()
nc = d.groupby("gene_id").apply(lambda g: g[["cell_type", "donor"]].drop_duplicates().cell_type.nunique())
nd = d.groupby("gene_id").apply(lambda g: g[["cell_type", "donor"]].drop_duplicates().donor.nunique())
genes = [g for g in ni.index if ni[g] >= 2 and nc.get(g, 0) >= 2 and nd.get(g, 0) >= 3]
if UNIT:
    genes = genes[:UNIT]
print(f"testing {len(genes):,} multi-PAS genes, K={K}{' [NEGCTRL]' if NEGCTRL else ''}")

rng = np.random.RandomState(0)
rows = []
d["cd_key"] = d.cell_type + "||" + d.donor
for g in genes:
    sub = d[d.gene_id == g]
    piv = sub.pivot_table(index="cd_key", columns="pas_id", values="count", fill_value=0)
    order = list(piv.index)
    meta = sub[["cd_key", "cell_type", "donor"]].drop_duplicates().set_index("cd_key").loc[order]
    dn0 = meta.donor.values
    dep = piv.values.sum(axis=1).astype(float)
    U = piv.values.astype(float) / dep[:, None]
    Uc = U.copy(); dgroups = []
    for dn in np.unique(dn0):
        m = np.where(dn0 == dn)[0]
        Uc[m] = U[m] - np.average(U[m], axis=0, weights=dep[m])
        dgroups.append(m)
    ccode, cuniq = pd.factorize(meta.cell_type.values)
    C = len(cuniq); W = dep.sum(); depUc = Uc * dep[:, None]
    def T_codes(codes):
        num = np.zeros((C, Uc.shape[1])); np.add.at(num, codes, depUc)
        den = np.bincount(codes, weights=dep, minlength=C)
        mr = np.divide(num, den[:, None], out=np.zeros_like(num), where=den[:, None] > 0)
        return float(np.sum((den / W) * np.sum(mr ** 2, axis=1)))
    Tobs = T_codes(ccode)
    ge = 0
    for _ in range(K):
        cp = ccode.copy()
        for m in dgroups:
            cp[m] = rng.permutation(ccode[m])
        if T_codes(cp) >= Tobs:
            ge += 1
    p = (1 + ge) / (K + 1)
    pooled = sub.groupby(["cell_type", "pas_id"])["count"].sum().reset_index()
    pv = pooled.pivot_table(index="cell_type", columns="pas_id", values="count", fill_value=0)
    pvn = pv.div(pv.sum(axis=1), axis=0)
    tv = 0.0; cts = list(pvn.index)
    for i in range(len(cts)):
        for j in range(i + 1, len(cts)):
            tv = max(tv, 0.5 * np.abs(pvn.iloc[i] - pvn.iloc[j]).sum())
    rows.append({"gene": g, "pval": p, "effect": float(tv), "n_pas": sub.pas_id.nunique()})
res = pd.DataFrame(rows)
res["qval"] = multipletests(res.pval, method="fdr_bh")[1]
res["sig"] = (res.qval < 0.05) & (res.effect >= 0.20)
if not UNIT and not NEGCTRL:
    res.sort_values("pval").to_csv(
        OUTPUT, sep="\t", index=False)
print(f"genes tested={len(res):,} | sig APA (q<0.05 & TV>=0.20)={res.sig.sum():,}")
print(f"frac raw p<0.05 = {(res.pval<0.05).mean():.1%} (calibration check)")
print("top:", res.sort_values('pval').head(6)[['gene','effect','pval','qval','n_pas']].to_string(index=False))
