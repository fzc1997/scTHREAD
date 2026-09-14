# Reproducibility environment

The scripts in this package do not redistribute raw sequencing data or large
reference bundles. Accessions in `database/frozen_release/release_manifest.tsv`
identify the public inputs. Set the following variables before running a
reducer or a validation workflow:

| Variable | Required content |
|---|---|
| `SCTHREAD_SCLONG_ROOT` | Local root containing the source registry, annotations and IsoQuant inputs |
| `SCTHREAD_EXTSSD_ROOT` | Workstation root hosting the IsoQuant tool env, GRCm39 reference and Fluidigm C1 inputs used by `reproducibility/processing_routes/run_fluidigm_c1_mouse.sh` |
| `SCTHREAD_REFERENCE_ROOT` | Root containing `10x_ref/` and `endpoint_atlas/` reference resources |
| `SCTHREAD_FIGURE_OUTPUT` | Output directory for figure renderers (default `figures/output/`) |
| `SCTHREAD_SCLONG_ROOT` | Root providing the sample-registry snapshot and the two processing-route data roots (defaults resolve inside this repository or its documented machine layout) |
| `SCTHREAD_RESULTS_ROOT` | Results root used by feature/reducer scripts |
| `SCTHREAD_GENOME_FASTA` | Assembly-matched indexed FASTA for splice-site features |
| `SCTHREAD_REFERENCE_JUNCTIONS` | Assembly-matched reference junction BED for PWM construction |
| `SCTHREAD_ANALYSIS_UNIT_REGISTRY` | Validated study-by-donor analysis-unit manifest |
| `SCTHREAD_SIGNED_ASE_ROOT` | ANCHOR signed allele-aware matrix root for the ASE reducer |

The default values are relative package paths where possible. An unset external
resource variable should fail with a missing-resource error rather than silently
using a researcher-specific absolute path.
