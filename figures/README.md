# Part 2 — How every figure was made

The manuscript has **three main figures** and **seven supplementary figures**
plus a graphical abstract. Each folder below holds the delivered figure file
(PDF) together with its renderer and QA script where code exists. Every
renderer reads only frozen inputs from [`../source_data/`](../source_data/)
and fails closed when an input is missing.

## Figure index

| Figure | Type | Delivered file | Renderer | Notes |
|---|---|---|---|---|
| **Figure 1** | Main — schematic | `main/fig1_schematic/NAR_Fig1_r15.pdf` | none | Hand-edited vector schematic; there is no generating code by design. The frozen numbers it prints (453 / 34 / 923,389) are validated by the SF1 QA instead. |
| **Figure 2** | Main | `main/fig2/NAR_Fig2_v4.pdf` | `main/fig2/render_nar_bio.py` (QA: `qa_nar_fig2_v3.py`) | Render: `python render_nar_bio.py --fig 2 --stem2 <out> --observed-suffix _9999`. Panel f is the three-axis evidence inventory. |
| **Figure 3** | Main | `main/fig3/NAR_Fig3_v8.pdf` | `main/fig3/render_nar_fig3_v8.py` (module chain v3–v7; QA: `qa_nar_fig3_v8.py`; data builders in `main/fig3/builders/`) | Render: `python render_nar_fig3_v8.py --stem <out> --observed-suffix _9999`. CD45 variable-exon evidence; panel b embeds a portal screenshot captured 2026-08-07 (shipped in `source_data/`). The builders that produced the CD45 figdata from BAMs are included for provenance; they need the source BAMs and are not re-runnable from this package. |
| **SF1–SF4** | Supplementary | `supplementary/sf1_sf4_sf7/NAR_SF{1,2,3,4}.pdf` | `supplementary/sf1_sf4_sf7/render_nar_sf_compact.py` (helper modules `render_nar_sf_all.py`, `render_nar_sf10_v3.py`, `render_nar_fig3_mouse_scont.py`; QA: `qa_nar_supplementary_refresh.py`) | Render all four: `python render_nar_sf_compact.py --outdir <dir>`. SF3 embeds six portal screenshots (2026-08-05, shipped in `source_data/portal_screenshots/`). |
| **SF5** | Supplementary | `supplementary/sf5_discovery/NAR_SF5_discovery_screen_20260828.pdf` (+svg/png) | `supplementary/sf5_discovery/render_nar_sf5_discovery_screen.py` (QA: `qa_nar_sf5_v2.py`) | Nine-locus discovery screen. CLI-driven; see repository root README for the exact command. |
| **SF6** | Supplementary | `supplementary/sf6_processing_flow/scTHREAD_Processing_Route_Flow.{pdf,svg,png}` | none | Hand-assembled vector flow diagram of the four processing routes and the statistical decision flow; the SVG is the editable source. |
| **SF7** | Supplementary | `supplementary/sf1_sf4_sf7/NAR_SF7.pdf` (+svg) | `supplementary/sf1_sf4_sf7/render_cd45_umi_sensitivity_figure.py` | CD45 read-versus-UMI sensitivity. Reads `evidence_layers/cd45_sensitivity/cd45_umi_sensitivity_20260826.tsv`; the recount itself is reproduced by the scripts in that directory. |
| **Graphical abstract** | — | `graphical_abstract/NAR_GraphicalAbstract_v2_20260805.pdf` | `graphical_abstract/render_nar_graphical_abstract.py` (QA: `qa_nar_graphical_abstract.py`) | Renders the schematic version of the abstract. |

## Conventions

- Renderers write to the directory given on the command line (default
  `figures/output/`, which is not tracked). Delivered PDFs are committed
  alongside their code.
- Shared figure style (palette, fonts, panel labels) lives in
  [`style/nar_style.py`](style/nar_style.py).
- Outputs are byte-reproducible from the frozen inputs except where a panel
  embeds a dated screenshot; those screenshots ship in `source_data/` so the
  render still reproduces the delivered figure.
