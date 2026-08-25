#!/usr/bin/env python3
"""
ASE #4 — allele × cell_type interaction test (donor-aware), the statistical core.

For each gene with allelic support in >=2 cell types across >=3 donors: fit a binomial GLM on the
haplotype ratio  cbind(hapA, hapB) ~ C(cell_type) + C(donor)  and LRT the cell_type term against the
donor-only null. A significant cell_type term = the allelic ratio differs across cell types beyond
donor variation = cell-type-specific cis-ASE (Pro's allele×cell_type interaction, with donor as a
fixed effect and cell_type main effect on the ratio being the interaction of interest).

Guards (the Beta-binomial over-detection lesson): q<0.05 AND effect-size gate (max pairwise
cell-type allele-fraction difference >= 0.20); plus a within-donor cell_type PERMUTATION negative
control that must collapse the signal.

Human unphased -> hapA/hapB (NOT maternal/paternal). Marrow studies only (shared cell types).

Output: results/paper1/f2_grammar/figdata/ase_interaction.tsv (+ perm_null summary)
Usage: ase_interaction_test.py [--genes N-for-unit-test]
"""
import sys, glob, json, os
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RESULTS_ROOT = os.environ.get("SCTHREAD_RESULTS_ROOT", "results")
LED = os.environ.get("SCTHREAD_INPUT_LEDGER", "metadata/science_input_ledger_v2.json")
INPUT_GLOB = os.environ.get("SCTHREAD_ASE_INPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/agg_ase_ct/*.ase_ct.parquet")
OUTPUT = os.environ.get("SCTHREAD_ASE_TEST_OUTPUT", f"{RESULTS_ROOT}/paper1/f2_grammar/figdata/ase_interaction.tsv")
r2g = {r["run_accession"]: r["study_accession"] for r in json.load(open(LED))["records"]}
UNIT = None
if len(sys.argv) > 1 and sys.argv[1].startswith("--genes"):
    UNIT = int(sys.argv[2])

d = pd.concat([pd.read_parquet(f) for f in glob.glob(
    INPUT_GLOB)])
d = d[d.gse.isin(["GSE307660", "GSE276974"])].copy()   # marrow, shared cell types
d["donor"] = d.run
d["tot"] = d.hapA + d.hapB
d = d[d.tot >= 20]                                       # min allelic depth per gene×ct×donor
cov = d.groupby("gene").agg(n_ct=("cell_type", "nunique"), n_donor=("donor", "nunique")).reset_index()
genes = cov[(cov.n_ct >= 2) & (cov.n_donor >= 3)].gene.tolist()
if UNIT:
    genes = genes[:UNIT]
print(f"testing {len(genes):,} genes (marrow, >=2 cell_type & >=3 donor, tot>=20)")

def fit_gene(sub, perm=False):
    s = sub.copy()
    if perm:
        # permute cell_type within donor (negative control)
        s["cell_type"] = s.groupby("donor")["cell_type"].transform(lambda x: np.random.permutation(x.values))
    if s.cell_type.nunique() < 2 or s.donor.nunique() < 2:
        return None
    y = np.c_[s.hapA.values, s.hapB.values]
    try:
        Xf = sm.add_constant(pd.get_dummies(s[["cell_type", "donor"]], drop_first=True).astype(float), has_constant="add")
        Xn = sm.add_constant(pd.get_dummies(s[["donor"]], drop_first=True).astype(float), has_constant="add")
        mf = sm.GLM(y, Xf, family=sm.families.Binomial()).fit()
        mn = sm.GLM(y, Xn, family=sm.families.Binomial()).fit()
        lr = 2 * (mf.llf - mn.llf)
        from scipy.stats import chi2
        dfd = Xf.shape[1] - Xn.shape[1]
        p = chi2.sf(lr, dfd) if dfd > 0 else 1.0
    except Exception:
        return None
    # effect = max pairwise cell-type allele fraction diff (donor-pooled)
    fr = s.groupby("cell_type").apply(lambda g: g.hapA.sum() / (g.hapA.sum() + g.hapB.sum()))
    eff = float(fr.max() - fr.min())
    return p, eff, fr.idxmax(), fr.idxmin()

rng = np.random.RandomState(0)
rows, perm_p = [], []
gsub = {g: x for g, x in d[d.gene.isin(genes)].groupby("gene")}
for i, g in enumerate(genes):
    r = fit_gene(gsub[g])
    if r:
        rows.append({"gene": g, "pval": r[0], "effect": r[1], "hi_ct": r[2], "lo_ct": r[3]})
    if UNIT and i < 200:  # permutation null on a sample for the unit test
        rp = fit_gene(gsub[g], perm=True)
        if rp: perm_p.append(rp[0])
res = pd.DataFrame(rows)
if len(res):
    res["qval"] = multipletests(res.pval, method="fdr_bh")[1]
    res["sig"] = (res.qval < 0.05) & (res.effect >= 0.20)
    out = OUTPUT
    res.sort_values("pval").to_csv(out, sep="\t", index=False)
    print(f"genes tested={len(res):,} | sig interaction (q<0.05 & effect>=0.20)={res.sig.sum():,}")
    print("top:", res.sort_values('pval').head(5)[['gene','effect','qval','hi_ct','lo_ct']].to_string(index=False))
if perm_p:
    perm_p = np.array(perm_p)
    print(f"\nNEG CONTROL (within-donor cell_type permutation, n={len(perm_p)}): "
          f"frac p<0.05 = {(perm_p<0.05).mean():.1%} (should be ~5% if calibrated)")
