# Changelog

## Unreleased development fixes

- Recognize `5UTR candidate sequence` and the observed
  `5UTR cadidate sequence` header automatically, including when an existing
  configuration still requests `UTR_sequence`.
- Treat SampleSheet and Top Unknown Barcodes inputs as optional selections;
  answering `no` no longer aborts project configuration.
- Separate result run names from software versions. New result folders use
  `YYYYMMDD_1pct`, `YYYYMMDD_full`, and `_runN` suffixes, while the exact
  pipeline version remains recorded in `manifest.json` and `run.log`.

## v1.3 — 2026-07-26

- Fixed Illumina section-style SampleSheet parsing when rows contain surplus
  trailing columns, preventing the `list.strip()` crash before index QC.
- Added automatic project inventory and recoverable cleanup with yes/no
  selection; no filenames need to be typed.
- Replaced the repeated numbered workspace layout with one simple
  `01_5UTR_Plasmid` project containing raw data, reference, 1% subset, results,
  config, and pipeline.
- Added one-command `1pct` and `full` runs with elapsed time, peak memory,
  output-size logging, `[7/7]` verification, and full-run resource estimates.
- Added base-R library coverage, rank-abundance, count-distribution, length,
  and GC plots.
- Standardized future version labels to two parts (`v1.3`, `v1.4`; major
  redesigns use `v2.0`).

## v1.2 — 2026-07-24

- Added direct PhiX174 classification from bundled reference.
- Added observed PhiX percentage and non-PhiX:PhiX ratio.
- Added target mapping rate with PhiX excluded from the denominator.
- Added PhiX-adjusted Undetermined fraction for index-failure assessment.
- Added 5′UTR/NextSeq QC and trimming interpretation guide.
- Simplified the repository structure and excluded generated test outputs.

## v1.1 — 2026-07-24

- Added configuration-driven reference-defined amplicon analysis.
- Added paired-end target reconstruction and exact/near assignment.
- Added library representation, dropout, coverage-uniformity, and index QC.
- Added HTML, TSV, manifest, and shareable text reports.
