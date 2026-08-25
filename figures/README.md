# Figure renderers

The included renderers use `SCTHREAD_PROJECT_ROOT` to locate the frozen
release tables. From this repository root:

```bash
export SCTHREAD_PROJECT_ROOT="$PWD"
python figures/scripts/render_nar_graphical_abstract.py --stem figures/graphical_abstract
python figures/scripts/render_nar_fig2.py
python figures/scripts/render_nar_fig3.py
```

The figures are evidence renderers, not a claim that all database records have
identical assay resolution. Their source tables and scope contracts are shipped
alongside the scripts.
