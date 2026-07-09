# NGS_LibraryQC

CLI-based NGS library QC tool for pooled amplicon libraries such as 5'UTR, promoter, enhancer, barcode, and mutagenesis libraries.

The tool extracts the sequence between two user-defined anchor sequences and reports:

- observed insert sequence counts
- optional reference/library member counts
- detected reference members
- dropout reference members
- non-reference sequences
- library representation metrics
- count distribution and ranked abundance plots

## Typical workflow

```text
FASTQ / FASTQ.gz
  -> find left/right anchors
  -> extract insert sequence
  -> count unique observed inserts
  -> compare to optional reference table
  -> report library representation QC
```

This is intended for library QC/counting, not conventional genomic variant calling.

## Install

```bash
conda create -n ngs_libraryqc python=3.11 -y
conda activate ngs_libraryqc
pip install -r requirements.txt
```

For offline workstations, download wheels on an internet-connected PC:

```bash
mkdir wheels
pip download -r requirements.txt -d wheels
```

Move the repository and `wheels/` directory to the workstation, then install offline:

```bash
conda create -n ngs_libraryqc python=3.11 -y
conda activate ngs_libraryqc
pip install --no-index --find-links wheels -r requirements.txt
```

## Recommended config-based usage

Edit `configs/example_5utr_config.json`, then run:

```bash
python ngs_libraryqc.py --config configs/example_5utr_config.json
```

or:

```bash
bash scripts/run_5utr_qc.sh
```

## Direct command-line usage

```bash
python ngs_libraryqc.py \
  --fastq sample_R1.fastq.gz sample_R2.fastq.gz \
  --left CTATAAAAGAGCTCACAACCCCTCA \
  --right GGAGGCCACACCCGCCACTCACCTG \
  --ref reference_5UTR.csv \
  --id-col utr_id \
  --seq-col utr_seq \
  --orientation both \
  --min-len 20 \
  --max-len 300 \
  --out-prefix results/sample_5utr
```

## Reference CSV format

```csv
utr_id,utr_seq
UTR_0001,ACGTTG...
UTR_0002,GGCTAA...
```

The column names can be changed in the config file or CLI.

## Output files

Without a reference CSV:

- `<prefix>.insert_counts.csv`
- `<prefix>.summary.txt`
- `<prefix>.count_distribution.png`
- `<prefix>.ranked_abundance.png`

With a reference CSV:

- `<prefix>.reference_counts.csv`
- `<prefix>.non_reference.csv`
- `<prefix>.insert_counts.csv`
- `<prefix>.summary.txt`
- `<prefix>.count_distribution.png`
- `<prefix>.ranked_abundance.png`

## Key QC metrics

- `detected_reference_count`: number of expected reference members detected
- `dropout_reference_count`: number of expected reference members not detected
- `anchor_found_rate`: fraction of reads with both anchors detected
- `exact_ref_match_rate_of_total`: fraction of total reads exactly matching a reference member
- `non_reference_rate_of_anchor_found`: fraction of anchor-found reads not matching reference
- `p90_p10_ratio_nonzero`: representation skew among detected members
- `gini_index`: inequality of representation; lower is more even
- `shannon_entropy`: diversity of the counted reference library

## Notes for staggered primers

If staggered N bases are outside the anchor sequence, no special handling is required. The parser finds the anchor first and extracts only the sequence between the two anchors.

For paired-end reads, use `--orientation both` so both the read sequence and its reverse complement are scanned.
