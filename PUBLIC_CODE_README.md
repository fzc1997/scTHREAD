# scTHREAD reproducibility release staging

This clean staging tree is based on the public `fzc1997/scTHREAD` data-only
repository and adds the reproducibility and figure entry points required by the
NAR revision. It is intended for a reviewed GitHub branch/tag and subsequent
Zenodo version archive; it has not yet been pushed.

## Included scope

- `reproducibility/`: DIU, APA, ASE, junction and biological-unit permutation
  reducers/validators, real examples, the CD45 CB/UB sensitivity analysis and
  environment/manifest definitions;
- `results/paper1/figdata/` and `results/paper1/f2_grammar/figdata/`: frozen
  source tables used by the included NAR figure renderers;
- `figures/scripts/`: path-configurable NAR Figure 1–3 and graphical-abstract
  renderers plus QA helpers;
- `tables/processing_manifest_public.tsv`: the redacted 453-run layer contract;
- frozen release tables, P0 observed/null artifacts, CD45 UMI and representative
  platform-sensitivity summaries.

The staging tree contains no internal `/gpfs`, `/Users` or `/home` paths in its
text/source files. Raw reads and private reference bundles are not included.
Run-local novel models and cross-platform checks retain their stated scope and
do not imply cross-run equivalence or platform invariance.

## Validation entry point

```bash
./workflows/run_public_validation.sh
```

The external GitHub/Zenodo publication step must be completed only after a
human review of the staged file whitelist, license and author metadata.
