#!/usr/bin/env python3
"""NGS_LibraryQC: anchor-based insert extraction and library representation counting."""

from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def clean_seq(x: str) -> str:
    """Normalize a DNA sequence string."""
    return str(x).strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")


def open_fastq(path: str):
    """Open plain or gzipped FASTQ."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def hamming_distance(a: str, b: str) -> Optional[int]:
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))


def find_anchor(seq: str, anchor: str, max_mismatch: int = 0) -> int:
    """Find an anchor sequence with optional Hamming mismatch allowance."""
    if max_mismatch == 0:
        return seq.find(anchor)

    k = len(anchor)
    for i in range(0, len(seq) - k + 1):
        dist = hamming_distance(seq[i : i + k], anchor)
        if dist is not None and dist <= max_mismatch:
            return i
    return -1


def extract_between_anchors(
    seq: str,
    left: str,
    right: str,
    max_mismatch: int = 0,
    min_len: int = 0,
    max_len: int = 10000,
) -> Optional[str]:
    """Return the sequence between left and right anchors."""
    left_pos = find_anchor(seq, left, max_mismatch)
    if left_pos < 0:
        return None

    start = left_pos + len(left)
    right_rel = find_anchor(seq[start:], right, max_mismatch)
    if right_rel < 0:
        return None

    insert = seq[start : start + right_rel]
    if len(insert) < min_len or len(insert) > max_len:
        return None
    return insert


def scan_seq(
    seq: str,
    left: str,
    right: str,
    max_mismatch: int,
    orientation: str,
    min_len: int,
    max_len: int,
) -> Tuple[Optional[str], str]:
    """Scan one read in forward and/or reverse-complement orientation."""
    if orientation in {"forward_only", "both"}:
        insert = extract_between_anchors(seq, left, right, max_mismatch, min_len, max_len)
        if insert is not None:
            return insert, "forward"

    if orientation in {"reverse_complement_only", "both"}:
        rc = str(Seq(seq).reverse_complement())
        insert = extract_between_anchors(rc, left, right, max_mismatch, min_len, max_len)
        if insert is not None:
            return insert, "reverse_complement"

    return None, "anchor_not_found"


def gini(values: Iterable[int]) -> float:
    arr = np.asarray(list(values), dtype=float)
    if arr.size == 0 or arr.sum() == 0:
        return float("nan")
    arr = np.sort(arr)
    n = arr.size
    return float((2 * np.sum(np.arange(1, n + 1) * arr) / (n * arr.sum())) - (n + 1) / n)


def shannon_entropy(values: Iterable[int]) -> float:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[arr > 0]
    if arr.size == 0:
        return float("nan")
    p = arr / arr.sum()
    return float(-np.sum(p * np.log2(p)))


def run_counting(
    fastq_paths: list[str],
    left: str,
    right: str,
    ref_path: Optional[str] = None,
    id_col: str = "utr_id",
    seq_col: str = "utr_seq",
    orientation: str = "both",
    anchor_mismatch: int = 0,
    min_len: int = 0,
    max_len: int = 10000,
) -> tuple[pd.DataFrame, Optional[pd.DataFrame], pd.DataFrame, dict]:
    left = clean_seq(left)
    right = clean_seq(right)

    ref_df = None
    seq_to_id = {}
    if ref_path:
        ref_df = pd.read_csv(ref_path)
        ref_df[id_col] = ref_df[id_col].astype(str)
        ref_df[seq_col] = ref_df[seq_col].map(clean_seq)
        seq_to_id = dict(zip(ref_df[seq_col], ref_df[id_col]))

    stats = Counter()
    insert_counts = Counter()
    ref_counts = Counter()
    nonref_counts = Counter()

    for fastq in fastq_paths:
        with open_fastq(fastq) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                stats["total_reads"] += 1
                seq = clean_seq(str(record.seq))
                insert, status = scan_seq(seq, left, right, anchor_mismatch, orientation, min_len, max_len)

                if insert is None:
                    stats["anchor_not_found"] += 1
                    continue

                stats["anchor_found"] += 1
                stats[f"orientation_{status}"] += 1
                insert_counts[insert] += 1

                if ref_df is not None:
                    if insert in seq_to_id:
                        ref_counts[seq_to_id[insert]] += 1
                        stats["exact_ref_match"] += 1
                    else:
                        nonref_counts[insert] += 1
                        stats["non_reference"] += 1

    insert_df = pd.DataFrame(insert_counts.items(), columns=["insert_seq", "count"])
    if insert_df.empty:
        insert_df = pd.DataFrame(columns=["insert_seq", "count", "frequency"])
    else:
        insert_df = insert_df.sort_values("count", ascending=False).reset_index(drop=True)
        insert_df["frequency"] = insert_df["count"] / insert_df["count"].sum()

    nonref_df = pd.DataFrame(nonref_counts.items(), columns=["sequence", "count"])
    if nonref_df.empty:
        nonref_df = pd.DataFrame(columns=["sequence", "count", "frequency"])
    else:
        nonref_df = nonref_df.sort_values("count", ascending=False).reset_index(drop=True)
        nonref_df["frequency"] = nonref_df["count"] / nonref_df["count"].sum()

    summary = {
        "total_reads": int(stats["total_reads"]),
        "anchor_found": int(stats["anchor_found"]),
        "anchor_found_rate": stats["anchor_found"] / stats["total_reads"] if stats["total_reads"] else 0,
        "anchor_not_found": int(stats["anchor_not_found"]),
        "unique_extracted_inserts": int(len(insert_counts)),
        "orientation_forward": int(stats["orientation_forward"]),
        "orientation_reverse_complement": int(stats["orientation_reverse_complement"]),
    }

    if ref_df is not None:
        ref_df["raw_count"] = ref_df[id_col].map(ref_counts).fillna(0).astype(int)
        total_ref_counts = int(ref_df["raw_count"].sum())
        ref_df["frequency"] = ref_df["raw_count"] / total_ref_counts if total_ref_counts else 0
        ref_df["detected"] = ref_df["raw_count"] > 0

        nonzero = ref_df.loc[ref_df["raw_count"] > 0, "raw_count"]
        p10 = np.percentile(nonzero, 10) if len(nonzero) > 10 else np.nan
        p90 = np.percentile(nonzero, 90) if len(nonzero) > 10 else np.nan

        summary.update(
            {
                "exact_ref_match": int(stats["exact_ref_match"]),
                "exact_ref_match_rate_of_total": stats["exact_ref_match"] / stats["total_reads"] if stats["total_reads"] else 0,
                "non_reference": int(stats["non_reference"]),
                "non_reference_rate_of_anchor_found": stats["non_reference"] / stats["anchor_found"] if stats["anchor_found"] else 0,
                "detected_reference_count": int(ref_df["detected"].sum()),
                "dropout_reference_count": int((~ref_df["detected"]).sum()),
                "median_count": float(ref_df["raw_count"].median()),
                "min_count": int(ref_df["raw_count"].min()),
                "max_count": int(ref_df["raw_count"].max()),
                "p90_p10_ratio_nonzero": float(p90 / p10) if p10 and not np.isnan(p10) else np.nan,
                "gini_index": gini(ref_df["raw_count"]),
                "shannon_entropy": shannon_entropy(ref_df["raw_count"]),
            }
        )

    return insert_df, ref_df, nonref_df, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract and count sequences between two anchors from FASTQ.")
    parser.add_argument("--fastq", nargs="+", required=True, help="FASTQ/FASTQ.gz files")
    parser.add_argument("--left", required=True, help="Left anchor sequence")
    parser.add_argument("--right", required=True, help="Right anchor sequence")
    parser.add_argument("--ref", default=None, help="Optional reference CSV")
    parser.add_argument("--id-col", default="utr_id", help="Reference ID column")
    parser.add_argument("--seq-col", default="utr_seq", help="Reference sequence column")
    parser.add_argument("--orientation", choices=["both", "forward_only", "reverse_complement_only"], default="both")
    parser.add_argument("--anchor-mismatch", type=int, default=0)
    parser.add_argument("--min-len", type=int, default=0)
    parser.add_argument("--max-len", type=int, default=10000)
    parser.add_argument("--out-prefix", default="ngs_libraryqc")
    args = parser.parse_args()

    insert_df, ref_df, nonref_df, summary = run_counting(
        fastq_paths=args.fastq,
        left=args.left,
        right=args.right,
        ref_path=args.ref,
        id_col=args.id_col,
        seq_col=args.seq_col,
        orientation=args.orientation,
        anchor_mismatch=args.anchor_mismatch,
        min_len=args.min_len,
        max_len=args.max_len,
    )

    prefix = Path(args.out_prefix)
    insert_df.to_csv(f"{prefix}.insert_counts.csv", index=False)
    if ref_df is not None:
        ref_df.to_csv(f"{prefix}.reference_counts.csv", index=False)
        nonref_df.to_csv(f"{prefix}.non_reference.csv", index=False)

    with open(f"{prefix}.summary.txt", "w") as handle:
        for key, value in summary.items():
            handle.write(f"{key}: {value}\n")

    print("Done.")
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
