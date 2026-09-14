# Processing-route entry scripts

Entry scripts for the two frozen-release processing routes that use in-house
software instead of a third-party single-cell demultiplexer. The third route
(10x scONT droplet runs, `wf-single-cell` -> IsoQuant) uses published
third-party tools and is documented in the release manifest rather than
re-distributed here.

| Script | Route | Scope in the frozen release |
|---|---|---|
| `pacbio_sc_extract.py` | PacBio barcode-resolved | Single-adapter 10x 3' cDNA CCS reads; splits R1 adapter -> 16 bp CB -> 12 bp UMI -> cDNA per read (35 runs) |
| `pacbio_sc_extract_kinnex.py` | PacBio barcode-resolved | Kinnex/MAS concatenated-array CCS reads; splits at every adapter and treats each segment as one molecule (7 runs) |
| `pacbio_sc_tag.py` | PacBio barcode-resolved | Tags extracted reads with `CB`/`UB` SAM tags for IsoQuant `--read_group tag:CB` |
| `run_smartseq2_isoquant.sh` | Sample/library-file-grouped | Plate-based long-read Smart-seq2: one FASTQ = one cell; per-cell `minimap2 -ax splice` then IsoQuant `--read_group file_name` (read counts, no UMI dedup) |
| `run_fluidigm_c1_mouse.sh` | Sample/library-file-grouped | Fluidigm C1 PacBio (non-10x, one SRR = one cell); `minimap2 -ax splice:hq` then IsoQuant `--read_group file_name`, mouse GRCm39 |

## PacBio cell-barcode correction (v2)

Both PacBio extractors implement 1-Hamming-distance whitelist correction of
the 16 bp cell barcode following HIT-scISOseq (Shi et al. 2023, Nat Commun
14:2385), which delegates to the Cell Ranger `correct_bc_error` posterior
model:

- prior = exact-hit whitelist frequency in the run; likelihood =
  prior x 10^(-min(Q_diff, 33)/10); posterior normalized over candidates;
- a barcode is corrected only when the argmax candidate exceeds
  `--cb-posterior-threshold` (default 0.975), otherwise the read/segment
  stays unassigned and is not written out;
- `--correct-cb` is on by default; `--exact-only` reproduces the legacy
  exact-match behaviour byte for byte;
- UMI correction is deliberately not performed.

`pacbio_sc_extract.py` additionally retries the reverse strand when the
forward-strand barcode misses the whitelist (first exact hit wins; in
correction mode the first strand passing the posterior threshold wins).

## Tool versions used for the frozen release

- minimap2 2.31-r1302, IsoQuant 3.13.1 (PacBio and sample/library-file routes);
  minimap2 2.24-r1122 for the 10x scONT droplet route inside `wf-single-cell`.
- 10x whitelists per run family: 737K-august-2016, 3M-february-2018,
  3M-3pgex-may-2023 (not redistributed here; see 10x Genomics).
- Accessions, per-run route assignment, whitelist assignment and output
  locations: `tables/processing_manifest_public.tsv` and the frozen release
  manifest.

## Environment variables

The two shell scripts resolve their data roots through variables defined in
`ENVIRONMENT.md (repository root)` (`SCTHREAD_SCLONG_ROOT`,
`SCTHREAD_EXTSSD_ROOT`). Defaults point at the internal machine layout used
for the frozen release; override them to run against your own copy of the
public accessions. See `PATCH_NOTES.md` for the exact differences between
these copies and the in-house sources.
