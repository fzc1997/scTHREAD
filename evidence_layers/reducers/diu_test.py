#!/usr/bin/env python3
"""
DIU #3 — cell-type differential isoform usage (donor-aware, permutation-calibrated).

Applies the ASE lesson proactively: NO parametric over-detection — use a directly-permutable
statistic + within-donor cell_type permutation for empirical p.

For each multi-isoform gene: per (transcript, cell_type, donor) counts -> isoform usage vector
(fraction of the gene's counts) per (cell_type, donor). Statistic = donor-adjusted between-cell-type
divergence in isoform composition (weighted L2 of donor-centered usage vectors). Permute cell_type
within donor -> null -> empirical p. Effect = max pairwise cell-type total-variation distance in
donor-pooled isoform composition. Gate: q<0.05 AND effect>=0.20.

Marrow only (shared cell types). Min gene depth per (cell_type,donor) >=20; >=2 isoforms; >=2
cell types; >=3 donors. Exclude 'Unassigned'.
Usage: diu_test.py [K] [--unit N]
"""
import sys, glob, os
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ct_harmonize import harmonize

K = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1000
UNIT = int(sys.argv[sys.argv.index("--unit") + 1]) if "--unit" in sys.argv else None
RESULTS_ROOT = os.environ.get("SCTHREAD_RESULTS_ROOT", "results")
INPUT_GLOB = os.environ.get("SCTHREAD_DIU_INPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/agg_diu/*.diu.parquet")
OUTPUT = os.environ.get("SCTHREAD_DIU_TEST_OUTPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/figdata/diu_celltype.tsv")

d = pd.concat([pd.read_parquet(f) for f in glob.glob(
    INPUT_GLOB)])
d = d[d.gse.isin(["GSE307660", "GSE276974"])].copy()
d["cell_type"] = harmonize(d.cell_type)          # cross-study vocab harmonization
d = d[d.cell_type != "Unassigned"].copy()
# re-aggregate: harmonization can merge multiple raw types (CD4 T + CD8 T -> T cell) within a run
d = d.groupby(["gene_id", "transcript_id", "cell_type", "run", "gse"], as_index=False)["count"].sum()
d["donor"] = d.run
# gene depth per (cell_type,donor)
gd = d.groupby(["gene_id", "cell_type", "donor"])["count"].sum().reset_index(name="gdepth")
d = d.merge(gd, on=["gene_id", "cell_type", "donor"])
d = d[d.gdepth >= 20]
# eligible genes
ni = d.groupby("gene_id").transcript_id.nunique()
nc = d.groupby(["gene_id"]).apply(lambda g: g[["cell_type", "donor"]].drop_duplicates().cell_type.nunique())
nd = d.groupby(["gene_id"]).apply(lambda g: g[["cell_type", "donor"]].drop_duplicates().donor.nunique())
genes = [g for g in ni.index if ni[g] >= 2 and nc.get(g, 0) >= 2 and nd.get(g, 0) >= 3]
if UNIT:
    genes = genes[:UNIT]
print(f"testing {len(genes):,} multi-isoform genes, K={K}")

def usage_stat(sub, ct):
    # build (cell_type,donor) x transcript usage; donor-adjust; between-ct weighted divergence
    piv = sub.pivot_table(index=["cd_key"], columns="transcript_id", values="count", fill_value=0)
    dep = piv.sum(axis=1).values
    U = piv.div(piv.sum(axis=1), axis=0).values           # usage fractions per (ct,donor)
    keys = piv.index.to_frame(index=False)
    keys["ct"] = ct; keys["donor"] = [k.split("||")[1] for k in piv.index]
    # donor-center usage (weighted by depth within donor)
    Uc = U.copy()
    for dn in np.unique(keys.donor):
        m = (keys.donor.values == dn)
        wm = np.average(U[m], axis=0, weights=dep[m])
        Uc[m] = U[m] - wm
    # between-ct weighted mean residual, L2
    W = dep.sum(); T = 0.0
    for c in np.unique(ct):
        m = (ct == c)
        if m.sum() == 0: continue
        mr = np.average(Uc[m], axis=0, weights=dep[m])
        T += (dep[m].sum() / W) * float(np.sum(mr ** 2))
    return T

rng = np.random.RandomState(0)
rows = []
d["cd_key"] = d.cell_type + "||" + d.donor
for g in genes:
    sub = d[d.gene_id == g]
    # pivot ONCE: (cd_key unit) x transcript; donor-centered usage Uc is FIXED across ct permutations
    piv = sub.pivot_table(index="cd_key", columns="transcript_id", values="count", fill_value=0)
    order = list(piv.index)
    meta = sub[["cd_key", "cell_type", "donor"]].drop_duplicates().set_index("cd_key").loc[order]
    dn0 = meta.donor.values
    dep = piv.values.sum(axis=1).astype(float)
    U = piv.values.astype(float) / dep[:, None]
    Uc = U.copy()
    dgroups = []
    for dn in np.unique(dn0):
        m = np.where(dn0 == dn)[0]
        Uc[m] = U[m] - np.average(U[m], axis=0, weights=dep[m])   # donor-center (FIXED)
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
    # effect: max pairwise TV distance in donor-pooled usage
    pooled = sub.groupby(["cell_type","transcript_id"])["count"].sum().reset_index()
    pv = pooled.pivot_table(index="cell_type", columns="transcript_id", values="count", fill_value=0)
    pvn = pv.div(pv.sum(axis=1), axis=0)
    tv = 0.0
    cts = list(pvn.index)
    for i in range(len(cts)):
        for j in range(i + 1, len(cts)):
            tv = max(tv, 0.5 * np.abs(pvn.iloc[i] - pvn.iloc[j]).sum())
    rows.append({"gene": g, "pval": p, "effect": float(tv), "n_iso": sub.transcript_id.nunique()})
res = pd.DataFrame(rows)
res["qval"] = multipletests(res.pval, method="fdr_bh")[1]
res["sig"] = (res.qval < 0.05) & (res.effect >= 0.20)
if not UNIT:
    res.sort_values("pval").to_csv(OUTPUT, sep="\t", index=False)
print(f"genes tested={len(res):,} | sig DIU (q<0.05 & TV>=0.20)={res.sig.sum():,}")
print(f"frac raw p<0.05 = {(res.pval<0.05).mean():.1%} (calibration check)")
print("top:", res.sort_values('pval').head(6)[['gene','effect','pval','qval','n_iso']].to_string(index=False))
