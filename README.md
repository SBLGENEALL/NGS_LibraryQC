# NGS_LibraryQC

Anchor-based NGS library QC tool for pooled amplicon libraries such as 5'UTR, promoter, enhancer, barcode, and mutagenesis libraries.

The tool extracts the sequence between two user-defined anchor sequences and reports:

- observed insert sequence counts
- optional reference/library member counts
- dropout variants
- non-reference sequences
- library representation metrics
- count distribution plots in the interactive app

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

## Interactive app

```bash
streamlit run app.py
```

In the sidebar, provide:

- FASTQ or FASTQ.gz files
- left anchor sequence
- right anchor sequence
- optional reference CSV
- optional insert length range
- read orientation: forward, reverse-complement, or both

## Command-line usage

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
  --out-prefix sample_5utr
```

## Reference CSV format

```csv
utr_id,utr_seq
UTR_0001,ACGTTG...
UTR_0002,GGCTAA...
```

The column names can be changed in the app or CLI.

## Output files

Without a reference CSV:

- `<prefix>.insert_counts.csv`
- `<prefix>.summary.txt`

With a reference CSV:

- `<prefix>.reference_counts.csv`
- `<prefix>.non_reference.csv`
- `<prefix>.insert_counts.csv`
- `<prefix>.summary.txt`

## Notes for staggered primers

If staggered N bases are outside the anchor sequence, no special handling is required. The parser finds the anchor first and extracts only the sequence between the two anchors.

For paired-end reads, use `--orientation both` so both the read sequence and its reverse complement are scanned.
