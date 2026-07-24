#!/usr/bin/env python3
"""
Configuration-driven amplicon NGS QC and designed-library quantification.

The script intentionally uses only the Python standard library so that it can
run on an offline corporate Linux server.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import gzip
import hashlib
import html
import json
import math
import re
import shutil
import sqlite3
import statistics
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import zip_longest
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


DNA_RE = re.compile(r"^[ACGTN+_-]+$", re.IGNORECASE)
FASTQ_SUFFIX_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)
READ_TOKEN_RE = re.compile(
    r"^(?P<prefix>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?"
    r"[_\.](?P<read>R?[12]|I[12])(?:[_\.](?P<chunk>\d+))?$",
    re.IGNORECASE,
)

PIPELINE_NAME = "Amplicon Library QC"
PIPELINE_VERSION = "1.2.0"


def normalize_dna(value: str, allow_n: bool = False) -> str:
    seq = re.sub(r"[^A-Za-z]", "", value).upper().replace("U", "T")
    allowed = set("ACGTN" if allow_n else "ACGT")
    if not seq:
        raise ValueError("empty DNA sequence")
    invalid = sorted(set(seq) - allowed)
    if invalid:
        raise ValueError(f"invalid DNA character(s): {','.join(invalid)}")
    return seq


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def open_text(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return cleaned or "sample"


def fmt_int(value: float | int | None) -> str:
    if value is None:
        return "NA"
    return f"{int(round(value)):,}"


def fmt_float(value: float | None, digits: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def pct(numerator: float, denominator: float) -> Optional[float]:
    if denominator == 0:
        return None
    return 100.0 * numerator / denominator


def percentile(sorted_values: Sequence[float], p: float) -> Optional[float]:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    position = (len(sorted_values) - 1) * p
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return float(sorted_values[low])
    fraction = position - low
    return float(sorted_values[low] * (1 - fraction) + sorted_values[high] * fraction)


def gini(values: Sequence[int]) -> Optional[float]:
    if not values:
        return None
    sorted_values = sorted(max(0, int(x)) for x in values)
    total = sum(sorted_values)
    if total == 0:
        return None
    n = len(sorted_values)
    weighted = sum((i + 1) * x for i, x in enumerate(sorted_values))
    return (2 * weighted) / (n * total) - (n + 1) / n


def rankdata(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for i in order[start:end]:
            ranks[i] = rank
        start = end
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) != len(y) or len(x) < 3:
        return None
    mx = statistics.fmean(x)
    my = statistics.fmean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denominator == 0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / denominator


def spearman(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) != len(y) or len(x) < 3:
        return None
    return pearson(rankdata(x), rankdata(y))


def shannon_effective_size(values: Sequence[int]) -> Optional[float]:
    total = sum(values)
    if total <= 0:
        return None
    entropy = 0.0
    for value in values:
        if value > 0:
            p = value / total
            entropy -= p * math.log(p)
    return math.exp(entropy)


def gc_percent(seq: str) -> float:
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def bounded_levenshtein(a: str, b: str, limit: int) -> int:
    """Levenshtein distance with an early exit above limit."""
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for i, cb in enumerate(b, 1):
        current = [i]
        row_min = i
        for j, ca in enumerate(a, 1):
            value = min(
                current[j - 1] + 1,
                previous[j] + 1,
                previous[j - 1] + (ca != cb),
            )
            current.append(value)
            row_min = min(row_min, value)
        if row_min > limit:
            return limit + 1
        previous = current
    return previous[-1]


@dataclass
class ReferenceSequence:
    sequence: str
    ids: List[str]
    length: int
    gc: float


@dataclass
class ReferenceLibrary:
    unique: List[ReferenceSequence]
    id_to_sequence: Dict[str, str]
    sequence_to_index: Dict[str, int]
    duplicate_sequences: Dict[str, List[str]]
    source_rows: int


def read_reference(
    path: Path,
    id_column: str,
    sequence_column: str,
    left_primer: str,
    right_primer: str,
    reference_mode: str,
) -> ReferenceLibrary:
    rows: List[Tuple[str, str]] = []
    lower_name = path.name.lower()
    if lower_name.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz")):
        with open_text(path) as handle:
            current_id: Optional[str] = None
            parts: List[str] = []
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        rows.append((current_id, "".join(parts)))
                    current_id = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line)
            if current_id is not None:
                rows.append((current_id, "".join(parts)))
    else:
        with open_text(path) as handle:
            sample = handle.read(8192)
            handle.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            except csv.Error:
                dialect = csv.excel_tab if "\t" in sample else csv.excel
            reader = csv.DictReader(handle, dialect=dialect)
            if not reader.fieldnames:
                raise ValueError("reference table has no header")
            normalized_headers = {x.strip().lower(): x for x in reader.fieldnames}
            actual_id = normalized_headers.get(id_column.strip().lower())
            actual_seq = normalized_headers.get(sequence_column.strip().lower())
            if actual_id is None or actual_seq is None:
                raise ValueError(
                    "reference columns not found. "
                    f"Available columns: {', '.join(reader.fieldnames)}"
                )
            for line_number, row in enumerate(reader, 2):
                variant_id = (row.get(actual_id) or "").strip()
                raw_seq = (row.get(actual_seq) or "").strip()
                if not variant_id and not raw_seq:
                    continue
                if not variant_id or not raw_seq:
                    raise ValueError(
                        f"incomplete reference row {line_number}: ID or sequence missing"
                    )
                rows.append((variant_id, raw_seq))

    if not rows:
        raise ValueError("no reference sequences found")

    right_forward = reverse_complement(right_primer)
    left_binding_arm = left_primer[-15:]
    right_binding_arm = right_forward[:15]
    normalized: List[Tuple[str, str]] = []
    seen_ids = set()
    for variant_id, raw_seq in rows:
        if variant_id in seen_ids:
            raise ValueError(f"duplicate Variant_ID: {variant_id}")
        seen_ids.add(variant_id)
        seq = normalize_dna(raw_seq, allow_n=False)
        full_primer_flanked = seq.startswith(left_primer) and seq.endswith(right_forward)
        binding_arm_flanked = seq.startswith(left_binding_arm) and seq.endswith(
            right_binding_arm
        )
        trim_amplicon = reference_mode == "amplicon" or (
            reference_mode == "auto"
            and (full_primer_flanked or binding_arm_flanked)
        )
        if trim_amplicon:
            if full_primer_flanked:
                seq = seq[len(left_primer) : -len(right_forward)]
            elif binding_arm_flanked:
                seq = seq[len(left_binding_arm) : -len(right_binding_arm)]
            else:
                raise ValueError(
                    f"could not extract target between primer sequences or 15-bp "
                    f"binding arms for reference {variant_id}"
                )
            if not seq:
                raise ValueError(f"empty target after flank removal: {variant_id}")
        normalized.append((variant_id, seq))

    sequence_ids: Dict[str, List[str]] = defaultdict(list)
    id_to_sequence: Dict[str, str] = {}
    for variant_id, seq in normalized:
        sequence_ids[seq].append(variant_id)
        id_to_sequence[variant_id] = seq

    unique = [
        ReferenceSequence(
            sequence=seq,
            ids=ids,
            length=len(seq),
            gc=gc_percent(seq),
        )
        for seq, ids in sequence_ids.items()
    ]
    unique.sort(key=lambda x: (x.ids[0], x.sequence))
    sequence_to_index = {item.sequence: i for i, item in enumerate(unique)}
    duplicates = {seq: ids for seq, ids in sequence_ids.items() if len(ids) > 1}
    return ReferenceLibrary(
        unique=unique,
        id_to_sequence=id_to_sequence,
        sequence_to_index=sequence_to_index,
        duplicate_sequences=duplicates,
        source_rows=len(normalized),
    )


@dataclass
class FastqFile:
    path: Path
    sample: str
    read: str
    lane: str
    chunk: str
    key: str


@dataclass
class FastqGroup:
    sample: str
    chunks: List[Tuple[Optional[Path], Optional[Path]]] = field(default_factory=list)
    is_undetermined: bool = False


def parse_fastq_name(path: Path) -> Optional[FastqFile]:
    stem = FASTQ_SUFFIX_RE.sub("", path.name)
    match = READ_TOKEN_RE.match(stem)
    if not match:
        # Common fallback: locate _R1_/_R2_/_I1_/_I2_ anywhere.
        fallback = re.search(r"(?P<sep>[_\.])(?P<read>[RI][12])(?P<tail>[_\.]\d+)?$", stem)
        if not fallback:
            return None
        prefix = stem[: fallback.start()]
        read = fallback.group("read").upper()
        lane = ""
        chunk = (fallback.group("tail") or "").strip("_.")
    else:
        prefix = match.group("prefix")
        read = match.group("read").upper()
        if read in {"1", "2"}:
            read = "R" + read
        lane = match.group("lane") or ""
        chunk = match.group("chunk") or ""

    sample = re.sub(r"_S\d+$", "", prefix)
    key = f"{sample}|{lane}|{chunk}"
    return FastqFile(path=path, sample=sample, read=read, lane=lane, chunk=chunk, key=key)


def discover_fastqs(input_dir: Path) -> Tuple[List[FastqGroup], List[FastqFile], List[str]]:
    parsed: List[FastqFile] = []
    warnings: List[str] = []
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or not FASTQ_SUFFIX_RE.search(path.name):
            continue
        item = parse_fastq_name(path)
        if item is None:
            warnings.append(f"Unrecognized FASTQ filename: {path}")
            continue
        parsed.append(item)

    by_key: Dict[str, Dict[str, Path]] = defaultdict(dict)
    sample_for_key: Dict[str, str] = {}
    for item in parsed:
        if item.read in {"R1", "R2"}:
            if item.read in by_key[item.key]:
                warnings.append(f"Duplicate {item.read} for chunk {item.key}")
            by_key[item.key][item.read] = item.path
            sample_for_key[item.key] = item.sample

    groups_map: Dict[str, FastqGroup] = {}
    for key in sorted(by_key):
        sample = sample_for_key[key]
        group = groups_map.setdefault(
            sample,
            FastqGroup(
                sample=sample,
                is_undetermined=sample.lower().startswith("undetermined"),
            ),
        )
        group.chunks.append((by_key[key].get("R1"), by_key[key].get("R2")))

    groups = sorted(groups_map.values(), key=lambda x: (x.is_undetermined, x.sample))
    return groups, parsed, warnings


def iter_fastq(path: Path) -> Iterator[Tuple[str, str, str]]:
    with open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            line_number += 4
            if not quality:
                raise ValueError(f"truncated FASTQ: {path} near line {line_number}")
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n").upper()
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ structure: {path} near line {line_number}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"sequence/quality length mismatch: {path} near line {line_number}"
                )
            yield header, sequence, quality


def read_single_fasta(path: Path) -> str:
    parts: List[str] = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    if not parts:
        raise ValueError(f"no sequence found in FASTA: {path}")
    return normalize_dna("".join(parts), allow_n=False)


class PhiXDetector:
    """Fast PhiX174 classifier using exact circular-genome k-mers."""

    def __init__(
        self,
        reference_path: Path,
        kmer_length: int = 27,
        min_kmer_hits: int = 1,
        scan_stride: int = 7,
    ):
        self.reference_path = reference_path
        self.kmer_length = kmer_length
        self.min_kmer_hits = min_kmer_hits
        self.scan_stride = scan_stride
        genome = read_single_fasta(reference_path)
        if len(genome) < kmer_length:
            raise ValueError(
                f"PhiX reference is shorter than k-mer length {kmer_length}"
            )
        circular = genome + genome[: kmer_length - 1]
        reverse = reverse_complement(genome)
        circular_reverse = reverse + reverse[: kmer_length - 1]
        self.kmers = {
            circular[i : i + kmer_length]
            for i in range(len(genome))
        }
        self.kmers.update(
            circular_reverse[i : i + kmer_length]
            for i in range(len(genome))
        )
        self.genome_length = len(genome)

    def is_phix(self, *sequences: Optional[str]) -> bool:
        hits = 0
        k = self.kmer_length
        for sequence in sequences:
            if not sequence or len(sequence) < k:
                continue
            last = len(sequence) - k
            positions = range(0, last + 1, self.scan_stride)
            for pos in positions:
                if sequence[pos : pos + k] in self.kmers:
                    hits += 1
                    if hits >= self.min_kmer_hits:
                        return True
            if last % self.scan_stride and sequence[last : last + k] in self.kmers:
                hits += 1
                if hits >= self.min_kmer_hits:
                    return True
        return False


def phred_values(quality: str) -> List[int]:
    return [max(0, ord(char) - 33) for char in quality]


@dataclass
class CycleMetric:
    bases: Counter = field(default_factory=Counter)
    total: int = 0
    q20: int = 0
    q30: int = 0
    qsum: int = 0


@dataclass
class FastqMetrics:
    label: str
    path: str
    sample: str
    read: str
    total_reads: int = 0
    total_bases: int = 0
    q20_bases: int = 0
    q30_bases: int = 0
    qsum: int = 0
    n_bases: int = 0
    reads_with_n: int = 0
    min_length: Optional[int] = None
    max_length: int = 0
    cycles: List[CycleMetric] = field(default_factory=list)

    def update(self, sequence: str, quality: str) -> None:
        qs = phred_values(quality)
        length = len(sequence)
        self.total_reads += 1
        self.total_bases += length
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = max(self.max_length, length)
        n_count = sequence.count("N")
        self.n_bases += n_count
        self.reads_with_n += int(n_count > 0)
        self.q20_bases += sum(q >= 20 for q in qs)
        self.q30_bases += sum(q >= 30 for q in qs)
        self.qsum += sum(qs)
        while len(self.cycles) < length:
            self.cycles.append(CycleMetric())
        for i, (base, q) in enumerate(zip(sequence, qs)):
            cycle = self.cycles[i]
            cycle.bases[base if base in "ACGTN" else "N"] += 1
            cycle.total += 1
            cycle.q20 += int(q >= 20)
            cycle.q30 += int(q >= 30)
            cycle.qsum += q

    def summary(self) -> Dict[str, object]:
        return {
            "label": self.label,
            "sample": self.sample,
            "read": self.read,
            "path": self.path,
            "total_reads": self.total_reads,
            "total_bases": self.total_bases,
            "mean_length": self.total_bases / self.total_reads if self.total_reads else None,
            "min_length": self.min_length,
            "max_length": self.max_length or None,
            "mean_q": self.qsum / self.total_bases if self.total_bases else None,
            "q20_base_percent": pct(self.q20_bases, self.total_bases),
            "q30_base_percent": pct(self.q30_bases, self.total_bases),
            "n_base_percent": pct(self.n_bases, self.total_bases),
            "reads_with_n_percent": pct(self.reads_with_n, self.total_reads),
        }


def extract_header_index(header: str) -> Optional[str]:
    parts = header.split()
    if len(parts) < 2:
        return None
    candidate = parts[-1].split(":")[-1].upper()
    if 6 <= len(candidate.replace("+", "")) <= 40 and DNA_RE.match(candidate):
        return candidate.replace("-", "+").replace("_", "+")
    return None


class CappedCounter:
    """Approximate high-frequency counter with bounded memory."""

    def __init__(self, max_items: int = 200_000):
        self.max_items = max_items
        self.counter: Counter = Counter()
        self.total = 0
        self.pruned = False

    def update(self, key: str, count: int = 1) -> None:
        self.counter[key] += count
        self.total += count
        if len(self.counter) > self.max_items:
            self.counter = Counter(dict(self.counter.most_common(self.max_items // 2)))
            self.pruned = True


@dataclass
class ReadEvidence:
    prefix: Optional[str]
    prefix_q: Optional[List[int]]
    suffix: Optional[str]
    suffix_q: Optional[List[int]]
    full_sequences: List[Tuple[str, List[int]]]
    left_anchor: bool
    right_anchor: bool
    orientation: str


def find_anchor(
    sequence: str,
    anchor: str,
    search_start: int = 0,
    search_end: int = 70,
    max_mismatches: int = 1,
) -> Optional[Tuple[int, int]]:
    end_limit = min(len(sequence) - len(anchor), search_end)
    if end_limit < search_start:
        return None
    exact = sequence.find(anchor, search_start, end_limit + len(anchor) + 1)
    if exact >= 0 and exact <= end_limit:
        return exact, 0
    best: Optional[Tuple[int, int]] = None
    for start in range(search_start, end_limit + 1):
        distance = hamming(sequence[start : start + len(anchor)], anchor)
        if distance <= max_mismatches and (best is None or distance < best[1]):
            best = (start, distance)
            if distance == 0:
                break
    return best


def extract_forward(
    sequence: str,
    quality: str,
    primer: str,
    opposite_forward: str,
    anchor_mismatches: int,
) -> Tuple[Optional[str], Optional[List[int]], Optional[Tuple[str, List[int]]], bool]:
    found = find_anchor(sequence, primer, max_mismatches=anchor_mismatches)
    if found is None:
        return None, None, None, False
    start = found[0] + len(primer)
    body = sequence[start:]
    body_q = phred_values(quality[start:])
    opposite = find_anchor(
        body,
        opposite_forward,
        search_start=0,
        search_end=len(body),
        max_mismatches=anchor_mismatches,
    )
    full = None
    if opposite is not None:
        full_seq = body[: opposite[0]]
        full_q = body_q[: opposite[0]]
        if full_seq:
            full = (full_seq, full_q)
    return body, body_q, full, True


def extract_reverse(
    sequence: str,
    quality: str,
    reverse_primer: str,
    opposite_reverse: str,
    anchor_mismatches: int,
) -> Tuple[Optional[str], Optional[List[int]], Optional[Tuple[str, List[int]]], bool]:
    found = find_anchor(sequence, reverse_primer, max_mismatches=anchor_mismatches)
    if found is None:
        return None, None, None, False
    start = found[0] + len(reverse_primer)
    body = sequence[start:]
    body_q = phred_values(quality[start:])
    opposite = find_anchor(
        body,
        opposite_reverse,
        search_start=0,
        search_end=len(body),
        max_mismatches=anchor_mismatches,
    )
    full = None
    if opposite is not None:
        trimmed = body[: opposite[0]]
        trimmed_q = body_q[: opposite[0]]
        if trimmed:
            full = (reverse_complement(trimmed), list(reversed(trimmed_q)))
    suffix = reverse_complement(body)
    suffix_q = list(reversed(body_q))
    return suffix, suffix_q, full, True


def build_evidence(
    r1_seq: str,
    r1_qual: str,
    r2_seq: Optional[str],
    r2_qual: Optional[str],
    left_primer: str,
    right_primer: str,
    anchor_mismatches: int,
) -> ReadEvidence:
    right_forward = reverse_complement(right_primer)
    left_reverse = reverse_complement(left_primer)

    p1, p1q, full1, left1 = extract_forward(
        r1_seq, r1_qual, left_primer, right_forward, anchor_mismatches
    )
    s2 = s2q = None
    full2 = None
    right2 = False
    if r2_seq is not None and r2_qual is not None:
        s2, s2q, full2, right2 = extract_reverse(
            r2_seq, r2_qual, right_primer, left_reverse, anchor_mismatches
        )
    score_a = int(left1) + int(right2)

    p2 = p2q = s1 = s1q = None
    full3 = full4 = None
    left2 = right1 = False
    if r2_seq is not None and r2_qual is not None:
        p2, p2q, full3, left2 = extract_forward(
            r2_seq, r2_qual, left_primer, right_forward, anchor_mismatches
        )
        s1, s1q, full4, right1 = extract_reverse(
            r1_seq, r1_qual, right_primer, left_reverse, anchor_mismatches
        )
    score_b = int(left2) + int(right1)

    if score_b > score_a:
        fulls = [x for x in (full3, full4) if x is not None]
        return ReadEvidence(
            prefix=p2,
            prefix_q=p2q,
            suffix=s1,
            suffix_q=s1q,
            full_sequences=fulls,
            left_anchor=left2,
            right_anchor=right1,
            orientation="R2_forward",
        )
    fulls = [x for x in (full1, full2) if x is not None]
    return ReadEvidence(
        prefix=p1,
        prefix_q=p1q,
        suffix=s2,
        suffix_q=s2q,
        full_sequences=fulls,
        left_anchor=left1,
        right_anchor=right2,
        orientation="R1_forward",
    )


class ReferenceMatcher:
    def __init__(
        self,
        reference: ReferenceLibrary,
        signature_length: int,
        max_mismatches: int,
    ):
        self.reference = reference
        self.signature_length = signature_length
        self.max_mismatches = max_mismatches
        self.prefix_indexes: Dict[int, Dict[str, set]] = {}
        self.suffix_indexes: Dict[int, Dict[str, set]] = {}
        lengths = sorted(set([signature_length, min(12, signature_length)]), reverse=True)
        for k in lengths:
            prefix: Dict[str, set] = defaultdict(set)
            suffix: Dict[str, set] = defaultdict(set)
            for i, item in enumerate(reference.unique):
                if item.length >= k:
                    prefix[item.sequence[:k]].add(i)
                    suffix[item.sequence[-k:]].add(i)
            self.prefix_indexes[k] = prefix
            self.suffix_indexes[k] = suffix

    @staticmethod
    def _neighbor_hits(key: str, index: Dict[str, set]) -> set:
        hits = set(index.get(key, ()))
        bases = "ACGT"
        for pos, original in enumerate(key):
            for base in bases:
                if base != original:
                    neighbor = key[:pos] + base + key[pos + 1 :]
                    found = index.get(neighbor)
                    if found:
                        hits.update(found)
        return hits

    def _candidate_ids(self, evidence: ReadEvidence) -> set:
        for k in self.prefix_indexes:
            p_hits: Optional[set] = None
            s_hits: Optional[set] = None
            if evidence.prefix and len(evidence.prefix) >= k:
                p_hits = self._neighbor_hits(
                    evidence.prefix[:k], self.prefix_indexes[k]
                )
            if evidence.suffix and len(evidence.suffix) >= k:
                s_hits = self._neighbor_hits(
                    evidence.suffix[-k:], self.suffix_indexes[k]
                )
            if p_hits is not None and s_hits is not None:
                intersection = p_hits & s_hits
                if intersection:
                    return intersection
            one_sided = (p_hits or set()) | (s_hits or set())
            if one_sided:
                return one_sided
        return set()

    @staticmethod
    def _assemble(
        reference_length: int,
        prefix: Optional[str],
        prefix_q: Optional[List[int]],
        suffix: Optional[str],
        suffix_q: Optional[List[int]],
    ) -> Tuple[str, List[bool], int]:
        bases: List[Optional[str]] = [None] * reference_length
        quals = [-1] * reference_length
        disagreements = 0
        if prefix:
            n = min(reference_length, len(prefix))
            q_values = prefix_q or [0] * len(prefix)
            for pos in range(n):
                bases[pos] = prefix[pos]
                quals[pos] = q_values[pos] if pos < len(q_values) else 0
        if suffix:
            n = min(reference_length, len(suffix))
            start = reference_length - n
            seq_part = suffix[-n:]
            q_values = (suffix_q or [0] * len(suffix))[-n:]
            for offset, base in enumerate(seq_part):
                pos = start + offset
                q = q_values[offset] if offset < len(q_values) else 0
                if bases[pos] is None:
                    bases[pos] = base
                    quals[pos] = q
                elif bases[pos] != base:
                    disagreements += 1
                    if q > quals[pos]:
                        bases[pos] = base
                        quals[pos] = q
        covered = [x is not None for x in bases]
        assembled = "".join(x if x is not None else "N" for x in bases)
        return assembled, covered, disagreements

    def match(
        self,
        evidence: ReadEvidence,
    ) -> Tuple[str, Optional[int], Optional[str], Optional[int], int]:
        """
        Returns:
          match_type, reference_index, observed_target, distance, pair_disagreements
        """
        full_candidates: List[Tuple[str, List[int]]] = evidence.full_sequences
        for full_seq, _ in full_candidates:
            index = self.reference.sequence_to_index.get(full_seq)
            if index is not None:
                return "exact", index, full_seq, 0, 0

        candidate_ids = self._candidate_ids(evidence)
        if not candidate_ids:
            if evidence.left_anchor or evidence.right_anchor:
                return "no_candidate", None, None, None, 0
            return "anchor_missing", None, None, None, 0

        scored = []
        for i in candidate_ids:
            item = self.reference.unique[i]
            assembled, covered, disagreements = self._assemble(
                item.length,
                evidence.prefix,
                evidence.prefix_q,
                evidence.suffix,
                evidence.suffix_q,
            )
            coverage = sum(covered) / item.length
            if coverage < 0.95:
                continue
            distance = sum(
                covered[pos] and assembled[pos] != item.sequence[pos]
                for pos in range(item.length)
            )
            scored.append((distance, disagreements, i, assembled))

        if not scored:
            return "insufficient_evidence", None, None, None, 0
        scored.sort(key=lambda x: (x[0], x[1], x[2]))
        best_distance = scored[0][0]
        best_disagreement = scored[0][1]
        best = [row for row in scored if row[0] == best_distance]
        if len(best) != 1:
            return "ambiguous", None, scored[0][3], best_distance, best_disagreement
        _, disagreements, index, assembled = best[0]
        if best_distance == 0:
            return "exact", index, assembled, 0, disagreements
        if best_distance <= self.max_mismatches:
            return "near", index, assembled, best_distance, disagreements
        return "too_many_mismatches", None, assembled, best_distance, disagreements


class ObservedSequenceStore:
    def __init__(self, output_path: Path):
        self.output_path = output_path
        self.db_path = output_path.with_suffix(output_path.suffix + ".sqlite")
        if self.db_path.exists():
            self.db_path.unlink()
        self.connection = sqlite3.connect(self.db_path)
        self.connection.execute(
            """
            CREATE TABLE observations (
                sample TEXT NOT NULL,
                design_key TEXT NOT NULL,
                observed_sequence TEXT NOT NULL,
                match_type TEXT NOT NULL,
                count INTEGER NOT NULL,
                PRIMARY KEY(sample, design_key, observed_sequence, match_type)
            )
            """
        )
        self.buffer: Counter = Counter()
        self.flush_size = 50_000

    def add(
        self,
        sample: str,
        design_key: str,
        observed_sequence: str,
        match_type: str,
    ) -> None:
        self.buffer[(sample, design_key, observed_sequence, match_type)] += 1
        if len(self.buffer) >= self.flush_size:
            self.flush()

    def flush(self) -> None:
        if not self.buffer:
            return
        with self.connection:
            self.connection.executemany(
                """
                INSERT INTO observations
                    (sample, design_key, observed_sequence, match_type, count)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(sample, design_key, observed_sequence, match_type)
                DO UPDATE SET count = count + excluded.count
                """,
                [(*key, count) for key, count in self.buffer.items()],
            )
        self.buffer.clear()

    def close_and_export(self) -> None:
        self.flush()
        with gzip.open(self.output_path, "wt", encoding="utf-8") as out:
            out.write(
                "sample\tdesign_key\tobserved_target_sequence\tmatch_type\tcount\n"
            )
            query = """
                SELECT sample, design_key, observed_sequence, match_type, count
                FROM observations
                ORDER BY sample, count DESC, design_key, observed_sequence
            """
            for row in self.connection.execute(query):
                out.write("\t".join(map(str, row)) + "\n")
        self.connection.close()
        self.db_path.unlink(missing_ok=True)


@dataclass
class SampleResult:
    sample: str
    is_undetermined: bool
    total_units: int
    counts_exact: List[int]
    counts_near: List[int]
    categories: Counter
    anchor_left: int
    anchor_right: int
    anchors_both: int
    pair_header_mismatch: int
    orientation_counts: Counter
    distance_counts: Counter
    pair_disagreements: int
    phix_reads: int

    @property
    def counts_total(self) -> List[int]:
        return [a + b for a, b in zip(self.counts_exact, self.counts_near)]


def canonical_read_id(header: str) -> str:
    token = header.split()[0]
    return re.sub(r"/[12]$", "", token)


def process_group(
    group: FastqGroup,
    matcher: ReferenceMatcher,
    left_primer: str,
    right_primer: str,
    anchor_mismatches: int,
    max_reads: Optional[int],
    file_metrics: Dict[str, FastqMetrics],
    header_indexes: CappedCounter,
    observed_store: Optional[ObservedSequenceStore],
    phix_detector: Optional[PhiXDetector],
) -> SampleResult:
    n_refs = len(matcher.reference.unique)
    exact = [0] * n_refs
    near = [0] * n_refs
    categories: Counter = Counter()
    orientations: Counter = Counter()
    distances: Counter = Counter()
    total_units = 0
    left_count = right_count = both_count = 0
    header_mismatch = 0
    pair_disagreements = 0
    phix_reads = 0

    for r1_path, r2_path in group.chunks:
        if r1_path is None and r2_path is None:
            continue
        if r1_path is None:
            r1_path, r2_path = r2_path, None
        assert r1_path is not None
        metric_r1 = file_metrics[str(r1_path)]
        iterator_r1 = iter_fastq(r1_path)
        iterator_r2 = iter_fastq(r2_path) if r2_path is not None else None
        metric_r2 = file_metrics[str(r2_path)] if r2_path is not None else None

        if iterator_r2 is None:
            iterator: Iterable[
                Tuple[Tuple[str, str, str], Optional[Tuple[str, str, str]]]
            ] = ((record, None) for record in iterator_r1)
        else:
            iterator = zip_longest(iterator_r1, iterator_r2)

        for record1, record2 in iterator:
            if record1 is None or (iterator_r2 is not None and record2 is None):
                raise ValueError(
                    f"R1/R2 read-count mismatch in sample {group.sample}: "
                    f"{r1_path} vs {r2_path}"
                )
            h1, s1, q1 = record1
            metric_r1.update(s1, q1)
            index1 = extract_header_index(h1)
            if index1:
                header_indexes.update(index1)

            h2 = s2 = q2 = None
            if record2 is not None:
                h2, s2, q2 = record2
                assert metric_r2 is not None
                metric_r2.update(s2, q2)
                index2 = extract_header_index(h2)
                if index2 and index2 != index1:
                    header_indexes.update(index2)
                if canonical_read_id(h1) != canonical_read_id(h2):
                    header_mismatch += 1

            total_units += 1
            if phix_detector is not None and phix_detector.is_phix(s1, s2):
                phix_reads += 1
                categories["phix"] += 1
            else:
                evidence = build_evidence(
                    s1,
                    q1,
                    s2,
                    q2,
                    left_primer,
                    right_primer,
                    anchor_mismatches,
                )
                left_count += int(evidence.left_anchor)
                right_count += int(evidence.right_anchor)
                both_count += int(evidence.left_anchor and evidence.right_anchor)
                orientations[evidence.orientation] += 1

                match_type, ref_index, observed, distance, disagreements = matcher.match(
                    evidence
                )
                categories[match_type] += 1
                pair_disagreements += disagreements
                if distance is not None:
                    distances[distance] += 1
                if ref_index is not None:
                    if match_type == "exact":
                        exact[ref_index] += 1
                    elif match_type == "near":
                        near[ref_index] += 1
                    design_key = "|".join(matcher.reference.unique[ref_index].ids)
                else:
                    design_key = ""
                if observed_store is not None and observed and "N" not in observed:
                    observed_store.add(group.sample, design_key, observed, match_type)

            if max_reads is not None and total_units >= max_reads:
                break
        if max_reads is not None and total_units >= max_reads:
            break

    return SampleResult(
        sample=group.sample,
        is_undetermined=group.is_undetermined,
        total_units=total_units,
        counts_exact=exact,
        counts_near=near,
        categories=categories,
        anchor_left=left_count,
        anchor_right=right_count,
        anchors_both=both_count,
        pair_header_mismatch=header_mismatch,
        orientation_counts=orientations,
        distance_counts=distances,
        pair_disagreements=pair_disagreements,
        phix_reads=phix_reads,
    )


def process_unpaired_qc_file(
    item: FastqFile,
    metric: FastqMetrics,
    max_reads: Optional[int],
    header_indexes: CappedCounter,
) -> None:
    for number, (header, sequence, quality) in enumerate(iter_fastq(item.path), 1):
        metric.update(sequence, quality)
        if item.read not in {"I1", "I2"}:
            index = extract_header_index(header)
            if index:
                header_indexes.update(index)
        if max_reads is not None and number >= max_reads:
            break


def process_index_fastqs(
    parsed_files: Sequence[FastqFile],
    file_metrics: Dict[str, FastqMetrics],
    max_reads: Optional[int],
    raw_index_reads: CappedCounter,
) -> set:
    """QC I1/I2 FASTQs and count observed index sequences."""
    by_key: Dict[str, Dict[str, FastqFile]] = defaultdict(dict)
    for item in parsed_files:
        if item.read in {"I1", "I2"}:
            by_key[item.key][item.read] = item

    consumed = set()
    for key in sorted(by_key):
        i1_item = by_key[key].get("I1")
        i2_item = by_key[key].get("I2")
        if i1_item is None and i2_item is None:
            continue
        if i1_item is None:
            i1_item, i2_item = i2_item, None
        assert i1_item is not None
        consumed.add(str(i1_item.path))
        metric_i1 = file_metrics[str(i1_item.path)]
        iterator_i1 = iter_fastq(i1_item.path)
        if i2_item is not None:
            consumed.add(str(i2_item.path))
            metric_i2 = file_metrics[str(i2_item.path)]
            iterator_i2 = iter_fastq(i2_item.path)
            iterator = zip_longest(iterator_i1, iterator_i2)
        else:
            metric_i2 = None
            iterator = ((record, None) for record in iterator_i1)

        for number, (record1, record2) in enumerate(iterator, 1):
            if record1 is None or (i2_item is not None and record2 is None):
                raise ValueError(f"I1/I2 read-count mismatch for index chunk {key}")
            _, sequence1, quality1 = record1
            metric_i1.update(sequence1, quality1)
            observed = sequence1
            if record2 is not None:
                _, sequence2, quality2 = record2
                assert metric_i2 is not None
                metric_i2.update(sequence2, quality2)
                observed = f"{sequence1}+{sequence2}"
            raw_index_reads.update(observed)
            if max_reads is not None and number >= max_reads:
                break
    return consumed


def calculate_library_metrics(
    counts: Sequence[int],
    exact_counts: Sequence[int],
    near_counts: Sequence[int],
    total_units: int,
    references: Sequence[ReferenceSequence],
) -> Dict[str, object]:
    n = len(counts)
    total_mapped = sum(counts)
    sorted_counts = sorted(counts)
    nonzero = [x for x in counts if x > 0]
    p10 = percentile(sorted_counts, 0.10)
    p90 = percentile(sorted_counts, 0.90)
    p95 = percentile(sorted_counts, 0.95)
    p05 = percentile(sorted_counts, 0.05)
    lengths = [item.length for item in references]
    gcs = [item.gc for item in references]
    logcounts = [math.log10(x + 1) for x in counts]
    sorted_desc = sorted(counts, reverse=True)
    top1_n = max(1, math.ceil(n * 0.01))
    top10_n = max(1, math.ceil(n * 0.10))
    mean_count = statistics.fmean(counts) if counts else None
    median_count = statistics.median(counts) if counts else None
    stdev = statistics.pstdev(counts) if len(counts) > 1 else 0.0
    return {
        "designed_ids": sum(len(item.ids) for item in references),
        "unique_design_sequences": n,
        "total_read_pairs_or_reads": total_units,
        "mapped_reads": total_mapped,
        "exact_mapped_reads": sum(exact_counts),
        "near_mapped_reads": sum(near_counts),
        "mapped_fraction_percent": pct(total_mapped, total_units),
        "detected_1x": sum(x >= 1 for x in counts),
        "detected_1x_percent": pct(sum(x >= 1 for x in counts), n),
        "detected_10x": sum(x >= 10 for x in counts),
        "detected_10x_percent": pct(sum(x >= 10 for x in counts), n),
        "detected_100x": sum(x >= 100 for x in counts),
        "detected_100x_percent": pct(sum(x >= 100 for x in counts), n),
        "detected_500x": sum(x >= 500 for x in counts),
        "detected_500x_percent": pct(sum(x >= 500 for x in counts), n),
        "dropout_sequences": sum(x == 0 for x in counts),
        "mean_coverage": mean_count,
        "median_coverage": median_count,
        "min_nonzero_coverage": min(nonzero) if nonzero else None,
        "max_coverage": max(counts) if counts else None,
        "coverage_cv": stdev / mean_count if mean_count else None,
        "p05_coverage": p05,
        "p10_coverage": p10,
        "p90_coverage": p90,
        "p95_coverage": p95,
        "p90_over_p10": p90 / p10 if p10 and p90 is not None else None,
        "p95_over_p05": p95 / p05 if p05 and p95 is not None else None,
        "gini": gini(counts),
        "top_1_percent_read_share_percent": pct(sum(sorted_desc[:top1_n]), total_mapped),
        "top_10_percent_read_share_percent": pct(
            sum(sorted_desc[:top10_n]), total_mapped
        ),
        "effective_library_size": shannon_effective_size(counts),
        "spearman_length_logcount": spearman(lengths, logcounts),
        "spearman_gc_logcount": spearman(gcs, logcounts),
    }


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        "NA"
                        if value is None
                        else f"{value:.8g}"
                        if isinstance(value, float)
                        else value
                    )
                    for key, value in row.items()
                }
            )


def write_sample_outputs(
    outdir: Path,
    result: SampleResult,
    reference: ReferenceLibrary,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    sample_dir = outdir / "samples" / safe_name(result.sample)
    sample_dir.mkdir(parents=True, exist_ok=True)
    counts = result.counts_total
    metrics = calculate_library_metrics(
        counts,
        result.counts_exact,
        result.counts_near,
        result.total_units,
        reference.unique,
    )
    non_phix_units = max(0, result.total_units - result.phix_reads)
    mapped_reads = sum(counts)
    unresolved_non_phix = max(0, non_phix_units - mapped_reads)
    metrics.update(
        {
            "sample": result.sample,
            "is_undetermined": result.is_undetermined,
            "phix_reads": result.phix_reads,
            "phix_fraction_percent": pct(result.phix_reads, result.total_units),
            "non_phix_read_pairs_or_reads": non_phix_units,
            "mapped_fraction_non_phix_percent": pct(mapped_reads, non_phix_units),
            "unresolved_non_phix_reads": unresolved_non_phix,
            "unresolved_non_phix_percent": pct(unresolved_non_phix, non_phix_units),
            "library_to_phix_read_ratio": (
                non_phix_units / result.phix_reads if result.phix_reads else None
            ),
            "left_anchor_percent": pct(result.anchor_left, non_phix_units),
            "right_anchor_percent": pct(result.anchor_right, non_phix_units),
            "both_anchors_percent": pct(result.anchors_both, non_phix_units),
            "pair_header_mismatches": result.pair_header_mismatch,
            "pair_base_disagreements": result.pair_disagreements,
            "ambiguous_reads": result.categories.get("ambiguous", 0),
            "anchor_missing_reads": result.categories.get("anchor_missing", 0),
            "unassigned_reads": unresolved_non_phix,
        }
    )
    metric_order = [
        "sample",
        "is_undetermined",
        "designed_ids",
        "unique_design_sequences",
        "total_read_pairs_or_reads",
        "phix_reads",
        "phix_fraction_percent",
        "non_phix_read_pairs_or_reads",
        "library_to_phix_read_ratio",
        "mapped_reads",
        "exact_mapped_reads",
        "near_mapped_reads",
        "mapped_fraction_percent",
        "mapped_fraction_non_phix_percent",
        "unresolved_non_phix_reads",
        "unresolved_non_phix_percent",
        "detected_1x",
        "detected_1x_percent",
        "detected_10x",
        "detected_10x_percent",
        "detected_100x",
        "detected_100x_percent",
        "detected_500x",
        "detected_500x_percent",
        "dropout_sequences",
        "mean_coverage",
        "median_coverage",
        "min_nonzero_coverage",
        "max_coverage",
        "coverage_cv",
        "p05_coverage",
        "p10_coverage",
        "p90_coverage",
        "p95_coverage",
        "p90_over_p10",
        "p95_over_p05",
        "gini",
        "top_1_percent_read_share_percent",
        "top_10_percent_read_share_percent",
        "effective_library_size",
        "spearman_length_logcount",
        "spearman_gc_logcount",
        "left_anchor_percent",
        "right_anchor_percent",
        "both_anchors_percent",
        "pair_header_mismatches",
        "pair_base_disagreements",
        "ambiguous_reads",
        "anchor_missing_reads",
        "unassigned_reads",
    ]
    write_tsv(
        sample_dir / f"{safe_name(result.sample)}_library_metrics.tsv",
        ["metric", "value"],
        [{"metric": key, "value": metrics.get(key)} for key in metric_order],
    )

    total_mapped = sum(counts)
    variant_rows: List[Dict[str, object]] = []
    for i, item in enumerate(reference.unique):
        row = {
            "sequence_key": "|".join(item.ids),
            "variant_ids": "|".join(item.ids),
            "duplicate_id_count": len(item.ids),
            "target_sequence": item.sequence,
            "length": item.length,
            "gc_percent": item.gc,
            "exact_count": result.counts_exact[i],
            "near_count": result.counts_near[i],
            "total_count": counts[i],
            "frequency_percent": pct(counts[i], total_mapped),
            "reads_per_million": 1_000_000 * counts[i] / total_mapped
            if total_mapped
            else None,
            "detected": int(counts[i] > 0),
        }
        variant_rows.append(row)
    fields = [
        "sequence_key",
        "variant_ids",
        "duplicate_id_count",
        "target_sequence",
        "length",
        "gc_percent",
        "exact_count",
        "near_count",
        "total_count",
        "frequency_percent",
        "reads_per_million",
        "detected",
    ]
    write_tsv(
        sample_dir / f"{safe_name(result.sample)}_variant_counts.tsv",
        fields,
        variant_rows,
    )
    write_tsv(
        sample_dir / f"{safe_name(result.sample)}_dropout_variants.tsv",
        fields,
        (row for row in variant_rows if row["total_count"] == 0),
    )
    write_tsv(
        sample_dir / f"{safe_name(result.sample)}_top_variants.tsv",
        fields,
        sorted(
            variant_rows,
            key=lambda x: (-int(x["total_count"]), str(x["sequence_key"])),
        )[:100],
    )
    write_tsv(
        sample_dir / f"{safe_name(result.sample)}_classification.tsv",
        ["category", "count", "percent"],
        [
            {
                "category": key,
                "count": value,
                "percent": pct(value, result.total_units),
            }
            for key, value in sorted(result.categories.items())
        ],
    )
    return metrics, variant_rows


def aggregate_results(
    name: str,
    results: Sequence[SampleResult],
    reference: ReferenceLibrary,
) -> SampleResult:
    n = len(reference.unique)
    exact = [0] * n
    near = [0] * n
    categories: Counter = Counter()
    orientations: Counter = Counter()
    distances: Counter = Counter()
    for result in results:
        exact = [x + y for x, y in zip(exact, result.counts_exact)]
        near = [x + y for x, y in zip(near, result.counts_near)]
        categories.update(result.categories)
        orientations.update(result.orientation_counts)
        distances.update(result.distance_counts)
    return SampleResult(
        sample=name,
        is_undetermined=False,
        total_units=sum(x.total_units for x in results),
        counts_exact=exact,
        counts_near=near,
        categories=categories,
        anchor_left=sum(x.anchor_left for x in results),
        anchor_right=sum(x.anchor_right for x in results),
        anchors_both=sum(x.anchors_both for x in results),
        pair_header_mismatch=sum(x.pair_header_mismatch for x in results),
        orientation_counts=orientations,
        distance_counts=distances,
        pair_disagreements=sum(x.pair_disagreements for x in results),
        phix_reads=sum(x.phix_reads for x in results),
    )


def write_cross_sample_matrices(
    outdir: Path,
    reference: ReferenceLibrary,
    results: Sequence[SampleResult],
) -> None:
    """Write count and RPM matrices for downstream statistics and plotting."""
    combined_dir = outdir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    count_fields = [
        "sequence_key",
        "variant_ids",
        "target_sequence",
        "length",
        "gc_percent",
        *[safe_name(result.sample) for result in results],
    ]
    count_rows: List[Dict[str, object]] = []
    rpm_rows: List[Dict[str, object]] = []
    totals = {result.sample: sum(result.counts_total) for result in results}
    for index, item in enumerate(reference.unique):
        metadata: Dict[str, object] = {
            "sequence_key": "|".join(item.ids),
            "variant_ids": "|".join(item.ids),
            "target_sequence": item.sequence,
            "length": item.length,
            "gc_percent": item.gc,
        }
        count_row = dict(metadata)
        rpm_row = dict(metadata)
        for result in results:
            column = safe_name(result.sample)
            value = result.counts_total[index]
            count_row[column] = value
            rpm_row[column] = (
                1_000_000 * value / totals[result.sample]
                if totals[result.sample]
                else None
            )
        count_rows.append(count_row)
        rpm_rows.append(rpm_row)
    write_tsv(
        combined_dir / "variant_count_matrix.tsv",
        count_fields,
        count_rows,
    )
    write_tsv(
        combined_dir / "variant_rpm_matrix.tsv",
        count_fields,
        rpm_rows,
    )


def parse_sample_sheet(path: Optional[Path]) -> List[Dict[str, str]]:
    if path is None:
        return []
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().lower() == "[data]":
            data_start = i + 1
            break
    data_lines = [line for line in lines[data_start:] if line.strip()]
    if not data_lines:
        return []
    reader = csv.DictReader(data_lines)
    rows = []
    for row in reader:
        normalized = {
            (key or "").strip().lower(): (value or "").strip().upper()
            for key, value in row.items()
        }
        sample = (
            normalized.get("sample_id")
            or normalized.get("sampleid")
            or normalized.get("sample_name")
            or normalized.get("samplename")
            or ""
        )
        i7 = normalized.get("index") or normalized.get("index1") or ""
        i5 = normalized.get("index2") or ""
        i7 = re.sub(r"[^ACGTN]", "", i7)
        i5 = re.sub(r"[^ACGTN]", "", i5)
        if sample and i7:
            rows.append({"sample": sample, "i7": i7, "i5": i5})
    return rows


def parse_top_unknown(path: Optional[Path]) -> Counter:
    output: Counter = Counter()
    if path is None:
        return output
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            return output
        field_lookup = {name.lower().strip(): name for name in reader.fieldnames}
        sequence_fields = [
            original
            for normalized, original in field_lookup.items()
            if "barcode" in normalized or "index" in normalized
        ]
        count_fields = [
            original
            for normalized, original in field_lookup.items()
            if "count" in normalized or "read" in normalized or "cluster" in normalized
        ]
        for row in reader:
            sequence = None
            count = None
            ordered_sequence_keys = sequence_fields + [
                key for key in reader.fieldnames if key not in sequence_fields
            ]
            for key in ordered_sequence_keys:
                value = row.get(key)
                raw = (value or "").strip().upper()
                normalized = raw.replace("-", "+").replace("_", "+")
                if sequence is None and DNA_RE.match(normalized) and len(
                    normalized.replace("+", "")
                ) >= 6:
                    sequence = normalized
            ordered_count_keys = count_fields + [
                key for key in reader.fieldnames if key not in count_fields
            ]
            for key in ordered_count_keys:
                raw = (row.get(key) or "").strip()
                try:
                    numeric = int(float(raw.replace(",", "")))
                except ValueError:
                    continue
                if numeric >= 0:
                    count = numeric
                    break
            if sequence is not None:
                output[sequence] += count if count is not None else 1
    return output


def index_distance(observed: str, expected: str) -> int:
    observed_parts = observed.split("+")
    expected_parts = expected.split("+")
    if len(observed_parts) != len(expected_parts):
        return 999
    total = 0
    for observed_part, expected_part in zip(observed_parts, expected_parts):
        if len(observed_part) != len(expected_part):
            return 999
        total += hamming(observed_part, expected_part)
    return total


def run_index_qc(
    outdir: Path,
    sample_sheet_rows: List[Dict[str, str]],
    raw_index_reads: CappedCounter,
    header_indexes: CappedCounter,
    top_unknown: Counter,
    max_distance: int,
    assigned_read_units: int,
    undetermined_read_units: int,
    assigned_phix_read_units: int,
    undetermined_phix_read_units: int,
) -> Dict[str, object]:
    index_dir = outdir / "index_qc"
    index_dir.mkdir(parents=True, exist_ok=True)
    expected_records = []
    for row in sample_sheet_rows:
        as_given = row["i7"] + (f"+{row['i5']}" if row["i5"] else "")
        rc_i5 = row["i7"] + (
            f"+{reverse_complement(row['i5'])}" if row["i5"] else ""
        )
        expected_records.append((row["sample"], as_given, rc_i5))

    source = []
    for barcode, count in raw_index_reads.counter.most_common():
        source.append(
            {"barcode": barcode, "count": count, "source": "I1_I2_FASTQ"}
        )
    for barcode, count in header_indexes.counter.most_common():
        source.append(
            {"barcode": barcode, "count": count, "source": "FASTQ_header"}
        )
    for barcode, count in top_unknown.most_common():
        source.append(
            {"barcode": barcode, "count": count, "source": "Top_Unknown_Barcodes"}
        )
    write_tsv(
        index_dir / "top_observed_barcodes.tsv",
        ["barcode", "count", "source"],
        sorted(source, key=lambda x: (-int(x["count"]), str(x["barcode"])))[:5000],
    )

    if raw_index_reads.counter:
        analysis_observed = Counter(raw_index_reads.counter)
        analysis_source = "I1_I2_FASTQ"
    elif header_indexes.counter:
        analysis_observed = Counter(header_indexes.counter)
        analysis_source = "FASTQ_header"
    else:
        analysis_observed = Counter(top_unknown)
        analysis_source = "Top_Unknown_Barcodes"

    nearest_rows = []
    substitution_counts: Counter = Counter()
    distance_counts: Counter = Counter()
    near_total = 0
    exact_expected_total = 0
    g_expected_total = 0
    g_to_n = 0
    g_to_other = 0
    g_correct = 0
    orientation_votes = Counter()

    if expected_records:
        for barcode, count in analysis_observed.items():
            candidates = []
            for sample, as_given, rc_i5 in expected_records:
                if "+" not in barcode:
                    expected_i7 = as_given.split("+")[0]
                    d_as = index_distance(barcode, expected_i7)
                    candidates.append((d_as, sample, expected_i7, "I1_only"))
                else:
                    d_as = index_distance(barcode, as_given)
                    candidates.append((d_as, sample, as_given, "as_given"))
                    if rc_i5 != as_given:
                        d_rc = index_distance(barcode, rc_i5)
                        candidates.append(
                            (d_rc, sample, rc_i5, "i5_reverse_complement")
                        )
            candidates.sort(key=lambda x: (x[0], x[1], x[3]))
            best_distance = candidates[0][0]
            tied = [x for x in candidates if x[0] == best_distance]
            best = candidates[0]
            distance_counts[best_distance] += count
            orientation_votes[best[3]] += count
            if best_distance == 0:
                exact_expected_total += count
            if best_distance <= max_distance:
                near_total += count
                observed_parts = barcode.split("+")
                expected_parts = best[2].split("+")
                for part_number, (observed_part, expected_part) in enumerate(
                    zip(observed_parts, expected_parts), 1
                ):
                    for pos, (expected_base, observed_base) in enumerate(
                        zip(expected_part, observed_part), 1
                    ):
                        substitution_counts[
                            (
                                f"I{part_number}",
                                pos,
                                expected_base,
                                observed_base,
                            )
                        ] += count
                        if expected_base == "G":
                            g_expected_total += count
                            if observed_base == "G":
                                g_correct += count
                            elif observed_base == "N":
                                g_to_n += count
                            else:
                                g_to_other += count
            nearest_rows.append(
                {
                    "observed_barcode": barcode,
                    "count": count,
                    "nearest_sample": best[1],
                    "nearest_expected": best[2],
                    "distance": best_distance if best_distance < 999 else "NA",
                    "orientation": best[3],
                    "tied_best_candidates": len(tied),
                }
            )

    write_tsv(
        index_dir / "nearest_expected_barcodes.tsv",
        [
            "observed_barcode",
            "count",
            "nearest_sample",
            "nearest_expected",
            "distance",
            "orientation",
            "tied_best_candidates",
        ],
        sorted(nearest_rows, key=lambda x: (-int(x["count"]), str(x["observed_barcode"]))),
    )
    write_tsv(
        index_dir / "index_position_substitutions.tsv",
        ["index_read", "position", "expected_base", "observed_base", "weighted_reads"],
        [
            {
                "index_read": key[0],
                "position": key[1],
                "expected_base": key[2],
                "observed_base": key[3],
                "weighted_reads": value,
            }
            for key, value in sorted(substitution_counts.items())
        ],
    )

    total_sample_units = assigned_read_units + undetermined_read_units
    assigned_non_phix = max(0, assigned_read_units - assigned_phix_read_units)
    undetermined_non_phix = max(
        0, undetermined_read_units - undetermined_phix_read_units
    )
    total_non_phix = assigned_non_phix + undetermined_non_phix
    if analysis_source == "Top_Unknown_Barcodes":
        note = (
            "Substitution percentages describe only failed/unknown index reads "
            f"within the {max_distance}-base nearest-index threshold."
        )
    else:
        note = (
            f"Substitution percentages use {analysis_source} barcodes within the "
            f"{max_distance}-base nearest-index threshold. Barcodes farther from "
            "every expected index are excluded from the substitution denominator."
        )
    summary = {
        "expected_index_count": len(expected_records),
        "assigned_read_pairs_or_reads": assigned_read_units,
        "undetermined_read_pairs_or_reads": undetermined_read_units,
        "assigned_read_fraction_percent": pct(assigned_read_units, total_sample_units),
        "undetermined_read_fraction_percent": pct(
            undetermined_read_units, total_sample_units
        ),
        "assigned_phix_reads": assigned_phix_read_units,
        "undetermined_phix_reads": undetermined_phix_read_units,
        "assigned_non_phix_reads": assigned_non_phix,
        "undetermined_non_phix_reads": undetermined_non_phix,
        "undetermined_non_phix_fraction_percent": pct(
            undetermined_non_phix, total_non_phix
        ),
        "index_substitution_analysis_source": analysis_source,
        "raw_I1_I2_index_reads_counted": raw_index_reads.total,
        "raw_I1_I2_unique_retained": len(raw_index_reads.counter),
        "raw_I1_I2_counter_pruned": raw_index_reads.pruned,
        "header_index_reads_counted": header_indexes.total,
        "header_index_unique_retained": len(header_indexes.counter),
        "header_index_counter_pruned": header_indexes.pruned,
        "top_unknown_reads_counted": sum(top_unknown.values()),
        "observed_index_reads_analyzed": sum(analysis_observed.values()),
        "exact_expected_index_reads": exact_expected_total,
        f"within_{max_distance}_bp_expected_reads": near_total,
        "expected_G_position_weighted_reads_within_threshold": g_expected_total,
        "expected_G_to_G": g_correct,
        "expected_G_to_N": g_to_n,
        "expected_G_to_other": g_to_other,
        "expected_G_to_N_percent": pct(g_to_n, g_expected_total),
        "expected_G_to_other_percent": pct(g_to_other, g_expected_total),
        "preferred_i5_orientation": (
            orientation_votes.most_common(1)[0][0] if orientation_votes else "NA"
        ),
        "note": note,
    }
    write_tsv(
        index_dir / "index_summary.tsv",
        ["metric", "value"],
        [{"metric": key, "value": value} for key, value in summary.items()],
    )
    return summary


def write_qc_outputs(
    outdir: Path,
    metrics: Sequence[FastqMetrics],
) -> List[Dict[str, object]]:
    summary_rows = [metric.summary() for metric in metrics]
    fields = [
        "label",
        "sample",
        "read",
        "path",
        "total_reads",
        "total_bases",
        "mean_length",
        "min_length",
        "max_length",
        "mean_q",
        "q20_base_percent",
        "q30_base_percent",
        "n_base_percent",
        "reads_with_n_percent",
    ]
    write_tsv(outdir / "run_qc_summary.tsv", fields, summary_rows)
    cycle_rows = []
    for metric in metrics:
        for i, cycle in enumerate(metric.cycles, 1):
            row = {
                "label": metric.label,
                "sample": metric.sample,
                "read": metric.read,
                "cycle": i,
                "total_bases": cycle.total,
                "mean_q": cycle.qsum / cycle.total if cycle.total else None,
                "q20_percent": pct(cycle.q20, cycle.total),
                "q30_percent": pct(cycle.q30, cycle.total),
            }
            for base in "ACGTN":
                row[f"{base}_percent"] = pct(cycle.bases.get(base, 0), cycle.total)
            cycle_rows.append(row)
    write_tsv(
        outdir / "per_cycle_qc.tsv",
        [
            "label",
            "sample",
            "read",
            "cycle",
            "total_bases",
            "mean_q",
            "q20_percent",
            "q30_percent",
            "A_percent",
            "C_percent",
            "G_percent",
            "T_percent",
            "N_percent",
        ],
        cycle_rows,
    )
    return summary_rows


def rank_svg(counts: Sequence[int], width: int = 760, height: int = 260) -> str:
    if not counts or max(counts) == 0:
        return "<p>No mapped reads.</p>"
    ordered = sorted(counts, reverse=True)
    values = [math.log10(x + 1) for x in ordered]
    max_y = max(values) or 1
    left, top, right, bottom = 55, 15, 15, 35
    plot_w = width - left - right
    plot_h = height - top - bottom
    points = []
    for i, value in enumerate(values):
        x = left + plot_w * i / max(1, len(values) - 1)
        y = top + plot_h * (1 - value / max_y)
        points.append(f"{x:.1f},{y:.1f}")
    return (
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Rank abundance plot">'
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" '
        'stroke="#666"/>'
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
        f'y2="{top + plot_h}" stroke="#666"/>'
        f'<polyline fill="none" stroke="#2d6cdf" stroke-width="2" '
        f'points="{" ".join(points)}"/>'
        f'<text x="{left + plot_w / 2:.1f}" y="{height - 5}" '
        'text-anchor="middle" font-size="12">Variant rank</text>'
        f'<text x="14" y="{top + plot_h / 2:.1f}" text-anchor="middle" '
        'font-size="12" transform="rotate(-90 14 '
        f'{top + plot_h / 2:.1f})">log10(read count + 1)</text>'
        "</svg>"
    )


def html_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    head = "".join(f"<th>{html.escape(str(x))}</th>" for x in headers)
    body = []
    for row in rows:
        body.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(x))}</td>" for x in row)
            + "</tr>"
        )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def create_html_report(
    outdir: Path,
    reference: ReferenceLibrary,
    qc_rows: List[Dict[str, object]],
    sample_metrics: Dict[str, Dict[str, object]],
    sample_results: Dict[str, SampleResult],
    index_summary: Dict[str, object],
    warnings: Sequence[str],
    command: str,
    target_name: str,
) -> None:
    target_label = html.escape(target_name)
    qc_table_rows = [
        [
            row["sample"],
            row["read"],
            fmt_int(row["total_reads"]),
            fmt_float(row["mean_q"]),
            fmt_float(row["q30_base_percent"], 2),
            fmt_float(row["n_base_percent"], 4),
        ]
        for row in qc_rows
    ]
    sample_table_rows = []
    sections = []
    for sample, metrics in sample_metrics.items():
        sample_table_rows.append(
            [
                sample,
                fmt_int(metrics.get("total_read_pairs_or_reads")),
                fmt_float(metrics.get("phix_fraction_percent"), 2),
                fmt_float(metrics.get("mapped_fraction_percent"), 2),
                fmt_float(metrics.get("mapped_fraction_non_phix_percent"), 2),
                fmt_float(metrics.get("detected_1x_percent"), 2),
                fmt_int(metrics.get("dropout_sequences")),
                fmt_float(metrics.get("mean_coverage"), 1),
                fmt_float(metrics.get("gini"), 3),
                fmt_float(metrics.get("top_10_percent_read_share_percent"), 2),
            ]
        )
        result = sample_results[sample]
        sections.append(
            f"<h3>{html.escape(sample)}</h3>"
            + rank_svg(result.counts_total)
        )

    index_rows = [[key, value] for key, value in index_summary.items()]
    warning_html = (
        "<ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in warnings) + "</ul>"
        if warnings
        else "<p>No pipeline warnings.</p>"
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{target_label} amplicon NGS QC report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       margin: 2rem auto; max-width: 1100px; padding: 0 1rem; color: #172033; }}
h1, h2, h3 {{ color: #153b73; }}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0 2rem; }}
th, td {{ border-bottom: 1px solid #d8dee9; padding: .45rem .6rem; text-align: left; }}
th {{ background: #eef4fb; position: sticky; top: 0; }}
.note {{ background: #fff8df; border-left: 4px solid #e0a800; padding: .8rem 1rem; }}
.meta {{ color: #52606d; font-size: .9rem; }}
svg {{ width: 100%; max-width: 760px; background: #fafcff; }}
code {{ word-break: break-all; }}
</style>
</head>
<body>
<h1>{target_label} Amplicon NGS QC Report</h1>
<p class="meta">Generated {time.strftime("%Y-%m-%d %H:%M:%S")}</p>
<p>Designed IDs: {reference.source_rows:,}; unique sequences:
{len(reference.unique):,}; duplicated design sequences:
{len(reference.duplicate_sequences):,}.</p>

<h2>Run QC</h2>
{html_table(["Sample", "Read", "Reads", "Mean Q", "Q30 base %", "N base %"], qc_table_rows)}

<h2>Library and PhiX overview</h2>
{html_table(["Sample", "Read pairs/reads", "PhiX %", "Mapped % (all)",
             "Mapped % (non-PhiX)", "Detected %", "Dropouts", "Mean coverage",
             "Gini", "Top 10% read share %"], sample_table_rows)}

<div class="note">ALL_WITH_UNDETERMINED is a diagnostic reference. Do not use it
as the official sample count if another library sharing the same design sequences
was multiplexed in the run.</div>

<h2>Rank-abundance curves</h2>
{''.join(sections)}

<h2>Index QC</h2>
{html_table(["Metric", "Value"], index_rows)}

<h2>Warnings</h2>
{warning_html}

<h2>Reproducibility</h2>
<p><code>{html.escape(command)}</code></p>
</body>
</html>
"""
    (outdir / "report.html").write_text(document, encoding="utf-8")


def create_share_summary(
    outdir: Path,
    reference: ReferenceLibrary,
    qc_rows: List[Dict[str, object]],
    sample_metrics: Dict[str, Dict[str, object]],
    index_summary: Dict[str, object],
    sample_results: Dict[str, SampleResult],
    warnings: Sequence[str],
    args: argparse.Namespace,
) -> None:
    lines = [
        f"{args.target_name.upper()} AMPLICON NGS ANALYSIS SUMMARY",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "[RUN SETTINGS]",
        f"Reference IDs: {reference.source_rows}",
        f"Unique reference sequences: {len(reference.unique)}",
        f"Duplicate design sequences: {len(reference.duplicate_sequences)}",
        f"Forward primer: {args.left_primer}",
        f"Reverse primer: {args.right_primer}",
        f"Target label: {args.target_name}",
        f"Max sequence mismatches: {args.max_mismatches}",
        f"Max reads per sample: {args.max_reads or 'ALL'}",
        f"PhiX QC enabled: {args.phix_enabled}",
        f"PhiX reference: {args.phix_reference if args.phix_enabled else 'DISABLED'}",
        f"Expected PhiX percent: "
        f"{args.expected_phix_percent if args.expected_phix_percent is not None else 'NA'}",
        "",
        "[FASTQ QC]",
    ]
    for row in qc_rows:
        lines.append(
            f"{row['sample']} {row['read']}: reads={fmt_int(row['total_reads'])}; "
            f"meanQ={fmt_float(row['mean_q'], 2)}; "
            f"Q30={fmt_float(row['q30_base_percent'], 2)}%; "
            f"N={fmt_float(row['n_base_percent'], 4)}%"
        )

    lines.extend(["", "[INDEX QC]"])
    for key in [
        "expected_index_count",
        "assigned_read_pairs_or_reads",
        "undetermined_read_pairs_or_reads",
        "assigned_read_fraction_percent",
        "undetermined_read_fraction_percent",
        "assigned_phix_reads",
        "undetermined_phix_reads",
        "assigned_non_phix_reads",
        "undetermined_non_phix_reads",
        "undetermined_non_phix_fraction_percent",
        "index_substitution_analysis_source",
        "raw_I1_I2_index_reads_counted",
        "header_index_reads_counted",
        "top_unknown_reads_counted",
        "observed_index_reads_analyzed",
        "exact_expected_index_reads",
        f"within_{args.max_index_distance}_bp_expected_reads",
        "expected_G_position_weighted_reads_within_threshold",
        "expected_G_to_G",
        "expected_G_to_N",
        "expected_G_to_other",
        "expected_G_to_N_percent",
        "expected_G_to_other_percent",
        "preferred_i5_orientation",
    ]:
        if key in index_summary:
            value = index_summary[key]
            if isinstance(value, float):
                value = fmt_float(value, 3)
            lines.append(f"{key}: {value}")
    lines.append(str(index_summary.get("note", "")))

    lines.extend(["", "[PHIX QC]"])
    if args.phix_enabled:
        overall = sample_metrics.get("ALL_WITH_UNDETERMINED", {})
        observed = overall.get("phix_fraction_percent")
        lines.extend(
            [
                f"All FASTQs PhiX reads: {fmt_int(overall.get('phix_reads'))}",
                f"All FASTQs total reads: "
                f"{fmt_int(overall.get('total_read_pairs_or_reads'))}",
                f"Observed PhiX fraction: {fmt_float(observed, 4)}%",
                f"Observed non-PhiX:PhiX ratio: "
                f"{fmt_float(overall.get('library_to_phix_read_ratio'), 4)}:1",
            ]
        )
        if args.expected_phix_percent is not None and observed is not None:
            lines.append(
                f"Observed minus expected PhiX: "
                f"{fmt_float(observed - args.expected_phix_percent, 4)} percentage points"
            )
        if args.max_reads is not None:
            lines.append(
                "NOTE: max_reads was used per sample; the combined PhiX fraction "
                "does not represent the full run mixture."
            )
    else:
        lines.append("PhiX classification disabled.")

    lines.extend(["", "[LIBRARY QC]"])
    for sample, metrics in sample_metrics.items():
        result = sample_results[sample]
        lines.extend(
            [
                f"Sample: {sample}",
                f"  Read pairs/reads: {fmt_int(metrics.get('total_read_pairs_or_reads'))}",
                f"  PhiX reads: {fmt_int(metrics.get('phix_reads'))} "
                f"({fmt_float(metrics.get('phix_fraction_percent'), 4)}%)",
                f"  Exact mapped: {fmt_int(metrics.get('exact_mapped_reads'))}",
                f"  Near mapped: {fmt_int(metrics.get('near_mapped_reads'))}",
                f"  Mapped fraction (all reads): "
                f"{fmt_float(metrics.get('mapped_fraction_percent'), 3)}%",
                f"  Mapped fraction (non-PhiX): "
                f"{fmt_float(metrics.get('mapped_fraction_non_phix_percent'), 3)}%",
                f"  Unresolved non-PhiX: "
                f"{fmt_int(metrics.get('unresolved_non_phix_reads'))} "
                f"({fmt_float(metrics.get('unresolved_non_phix_percent'), 3)}%)",
                f"  Detected >=1x: {fmt_int(metrics.get('detected_1x'))}/"
                f"{fmt_int(metrics.get('unique_design_sequences'))} "
                f"({fmt_float(metrics.get('detected_1x_percent'), 3)}%)",
                f"  Detected >=10x: {fmt_int(metrics.get('detected_10x'))} "
                f"({fmt_float(metrics.get('detected_10x_percent'), 3)}%)",
                f"  Detected >=100x: {fmt_int(metrics.get('detected_100x'))} "
                f"({fmt_float(metrics.get('detected_100x_percent'), 3)}%)",
                f"  Detected >=500x: {fmt_int(metrics.get('detected_500x'))} "
                f"({fmt_float(metrics.get('detected_500x_percent'), 3)}%)",
                f"  Dropouts: {fmt_int(metrics.get('dropout_sequences'))}",
                f"  Mean/median coverage: "
                f"{fmt_float(metrics.get('mean_coverage'), 2)}/"
                f"{fmt_float(metrics.get('median_coverage'), 2)}",
                f"  Min nonzero/max coverage: "
                f"{fmt_int(metrics.get('min_nonzero_coverage'))}/"
                f"{fmt_int(metrics.get('max_coverage'))}",
                f"  Gini: {fmt_float(metrics.get('gini'), 4)}",
                f"  P90/P10: {fmt_float(metrics.get('p90_over_p10'), 3)}",
                f"  Top 10% read share: "
                f"{fmt_float(metrics.get('top_10_percent_read_share_percent'), 3)}%",
                f"  Effective library size: "
                f"{fmt_float(metrics.get('effective_library_size'), 1)}",
                f"  Spearman length vs logcount: "
                f"{fmt_float(metrics.get('spearman_length_logcount'), 4)}",
                f"  Spearman GC vs logcount: "
                f"{fmt_float(metrics.get('spearman_gc_logcount'), 4)}",
                f"  Both primer anchors: "
                f"{fmt_float(metrics.get('both_anchors_percent'), 3)}%",
                "  Classification: "
                + ", ".join(
                    f"{key}={value:,}" for key, value in sorted(result.categories.items())
                ),
            ]
        )

    lines.extend(["", "[WARNINGS]"])
    lines.extend(warnings or ["None"])
    lines.extend(
        [
            "",
            "[FILES TO CHECK]",
            "report.html",
            "run_qc_summary.tsv",
            "per_cycle_qc.tsv",
            "index_qc/index_summary.tsv",
            "index_qc/index_position_substitutions.tsv",
            "combined/ALL_ASSIGNED_variant_counts.tsv",
            "combined/ALL_WITH_UNDETERMINED_variant_counts.tsv",
            "combined/all_sample_metrics.tsv",
            "combined/variant_count_matrix.tsv",
            "combined/variant_rpm_matrix.tsv",
        ]
    )
    (outdir / "RESULTS_TO_SHARE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _config_path(
    value: str,
    config_dir: Path,
) -> Optional[Path]:
    value = value.strip()
    if not value:
        return None
    path = Path(value).expanduser()
    return path if path.is_absolute() else config_dir / path


def load_config(path: Path) -> Dict[str, object]:
    """Load an INI file and return argparse-compatible defaults."""
    if not path.is_file():
        raise ValueError(f"config file not found: {path}")
    parser = configparser.ConfigParser(interpolation=None)
    with path.open("r", encoding="utf-8") as handle:
        parser.read_file(handle)
    base = path.resolve().parent
    values: Dict[str, object] = {}

    def text(section: str, key: str, default: Optional[str] = None) -> Optional[str]:
        if not parser.has_option(section, key):
            return default
        value = parser.get(section, key).strip()
        return value if value else default

    def integer(section: str, key: str, default: Optional[int] = None) -> Optional[int]:
        if not parser.has_option(section, key):
            return default
        value = parser.get(section, key).strip()
        return int(value) if value else default

    def boolean(section: str, key: str, default: bool = False) -> bool:
        return parser.getboolean(section, key) if parser.has_option(section, key) else default

    values["target_name"] = text("project", "target_name", "target")
    for key in ("input_dir", "reference", "outdir", "sample_sheet", "top_unknown"):
        raw = text("input", key)
        values[key] = _config_path(raw, base) if raw else None
    values["id_column"] = text("reference", "id_column", "Variant_ID")
    values["sequence_column"] = text("reference", "sequence_column", "target_sequence")
    values["reference_mode"] = text("reference", "reference_mode", "auto")
    values["left_primer"] = text("amplicon", "left_primer")
    values["right_primer"] = text("amplicon", "right_primer")
    values["anchor_mismatches"] = integer("matching", "anchor_mismatches", 1)
    values["signature_length"] = integer("matching", "signature_length", 20)
    values["max_mismatches"] = integer("matching", "max_mismatches", 2)
    values["max_index_distance"] = integer("index_qc", "max_index_distance", 2)
    raw_phix_reference = text("control_qc", "phix_reference")
    values["phix_reference"] = (
        _config_path(raw_phix_reference, base) if raw_phix_reference else None
    )
    values["phix_enabled"] = boolean("control_qc", "phix_enabled", True)
    values["expected_phix_percent"] = (
        parser.getfloat("control_qc", "expected_phix_percent")
        if parser.has_option("control_qc", "expected_phix_percent")
        and parser.get("control_qc", "expected_phix_percent").strip()
        else None
    )
    values["phix_kmer_length"] = integer("control_qc", "phix_kmer_length", 27)
    values["phix_min_kmer_hits"] = integer(
        "control_qc", "phix_min_kmer_hits", 1
    )
    values["max_reads"] = integer("runtime", "max_reads")
    values["write_observed_sequences"] = boolean(
        "runtime", "write_observed_sequences", True
    )
    values["overwrite"] = boolean("runtime", "overwrite", False)
    return values


def build_parser(defaults: Optional[Dict[str, object]] = None) -> argparse.ArgumentParser:
    defaults = defaults or {}
    parser = argparse.ArgumentParser(
        description=(
            "QC and quantify a reference-defined amplicon library from "
            "single-end or paired-end FASTQ files."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="INI configuration file")
    parser.add_argument("--target-name", default="target")
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--id-column", default="Variant_ID")
    parser.add_argument("--sequence-column", default="target_sequence")
    parser.add_argument(
        "--reference-mode",
        choices=["auto", "target", "amplicon"],
        default="auto",
    )
    parser.add_argument("--left-primer")
    parser.add_argument("--right-primer")
    parser.add_argument("--anchor-mismatches", type=int, default=1)
    parser.add_argument("--signature-length", type=int, default=20)
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument("--sample-sheet", type=Path)
    parser.add_argument("--top-unknown", type=Path)
    parser.add_argument("--max-index-distance", type=int, default=2)
    parser.add_argument(
        "--phix-reference",
        type=Path,
        help="PhiX174 FASTA; bundled NC_001422.1 is used when omitted",
    )
    parser.add_argument("--expected-phix-percent", type=float)
    parser.add_argument("--phix-kmer-length", type=int, default=27)
    parser.add_argument("--phix-min-kmer-hits", type=int, default=1)
    parser.add_argument(
        "--phix-qc",
        dest="phix_enabled",
        action="store_true",
        default=True,
    )
    parser.add_argument("--no-phix-qc", dest="phix_enabled", action="store_false")
    parser.add_argument("--max-reads", type=int)
    parser.add_argument(
        "--write-observed-sequences",
        dest="write_observed_sequences",
        action="store_true",
    )
    parser.add_argument(
        "--no-write-observed-sequences",
        dest="write_observed_sequences",
        action="store_false",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.set_defaults(**defaults)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    for name in ("input_dir", "reference", "outdir", "left_primer", "right_primer"):
        if getattr(args, name, None) in (None, ""):
            raise ValueError(
                f"{name.replace('_', ' ')} is required in the config file "
                f"or as --{name.replace('_', '-')}"
            )
    if not args.input_dir.is_dir():
        raise ValueError(f"input directory not found: {args.input_dir}")
    if not args.reference.is_file():
        raise ValueError(f"reference file not found: {args.reference}")
    for optional in ("sample_sheet", "top_unknown"):
        value = getattr(args, optional)
        if value is not None and not value.is_file():
            raise ValueError(f"{optional.replace('_', ' ')} not found: {value}")
    if args.signature_length < 8:
        raise ValueError("--signature-length must be at least 8")
    if args.max_mismatches < 0 or args.max_mismatches > 5:
        raise ValueError("--max-mismatches must be between 0 and 5")
    if args.anchor_mismatches < 0 or args.anchor_mismatches > 3:
        raise ValueError("--anchor-mismatches must be between 0 and 3")
    if args.max_reads is not None and args.max_reads <= 0:
        raise ValueError("--max-reads must be positive")
    if args.expected_phix_percent is not None and not (
        0 <= args.expected_phix_percent <= 100
    ):
        raise ValueError("--expected-phix-percent must be between 0 and 100")
    if args.phix_kmer_length < 17 or args.phix_kmer_length > 51:
        raise ValueError("--phix-kmer-length must be between 17 and 51")
    if args.phix_min_kmer_hits < 1:
        raise ValueError("--phix-min-kmer-hits must be at least 1")
    if args.phix_enabled:
        if args.phix_reference is None:
            args.phix_reference = (
                Path(__file__).resolve().parent
                / "references"
                / "NC_001422.1_phiX174.fasta"
            )
        if not args.phix_reference.is_file():
            raise ValueError(f"PhiX reference not found: {args.phix_reference}")
    args.target_name = re.sub(r"\s+", " ", str(args.target_name).strip()) or "target"
    args.left_primer = normalize_dna(args.left_primer)
    args.right_primer = normalize_dna(args.right_primer)


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(argv) if argv is not None else sys.argv[1:]
    config_probe = argparse.ArgumentParser(add_help=False)
    config_probe.add_argument("--config", type=Path)
    config_args, _ = config_probe.parse_known_args(argv)
    defaults: Dict[str, object] = {}
    if config_args.config is not None:
        try:
            defaults = load_config(config_args.config)
        except (ValueError, configparser.Error) as error:
            config_probe.error(str(error))
    parser = build_parser(defaults)
    args = parser.parse_args(argv)
    try:
        validate_args(args)
    except ValueError as error:
        parser.error(str(error))

    if args.outdir.exists():
        if not args.overwrite:
            parser.error(
                f"output directory already exists: {args.outdir}. "
                "Use a new path or --overwrite."
            )
        shutil.rmtree(args.outdir)
    args.outdir.mkdir(parents=True)
    started = time.time()
    warnings: List[str] = []

    print("[1/7] Loading designed reference library", flush=True)
    reference = read_reference(
        args.reference,
        args.id_column,
        args.sequence_column,
        args.left_primer,
        args.right_primer,
        args.reference_mode,
    )
    if reference.duplicate_sequences:
        warnings.append(
            f"{len(reference.duplicate_sequences)} duplicated design sequence(s) "
            "cannot be distinguished by sequencing."
        )
    duplicate_rows = []
    for sequence, ids in reference.duplicate_sequences.items():
        duplicate_rows.append(
            {
                "target_sequence": sequence,
                "variant_ids": "|".join(ids),
                "id_count": len(ids),
            }
        )
    write_tsv(
        args.outdir / "design_sequence_duplicates.tsv",
        ["target_sequence", "variant_ids", "id_count"],
        duplicate_rows,
    )

    print("[2/7] Discovering FASTQ files", flush=True)
    groups, parsed_files, discovery_warnings = discover_fastqs(args.input_dir)
    warnings.extend(discovery_warnings)
    if not groups:
        parser.error("no R1/R2 FASTQ files found under --input-dir")

    file_metrics: Dict[str, FastqMetrics] = {}
    for item in parsed_files:
        label = f"{item.sample}:{item.read}:{item.lane or '-'}:{item.chunk or '-'}"
        file_metrics[str(item.path)] = FastqMetrics(
            label=label,
            path=str(item.path),
            sample=item.sample,
            read=item.read,
        )

    matcher = ReferenceMatcher(
        reference=reference,
        signature_length=args.signature_length,
        max_mismatches=args.max_mismatches,
    )
    phix_detector = (
        PhiXDetector(
            args.phix_reference,
            kmer_length=args.phix_kmer_length,
            min_kmer_hits=args.phix_min_kmer_hits,
        )
        if args.phix_enabled
        else None
    )
    observed_store = (
        ObservedSequenceStore(args.outdir / "observed_target_sequences.tsv.gz")
        if args.write_observed_sequences
        else None
    )
    header_indexes = CappedCounter()
    raw_index_reads = CappedCounter()
    sample_results_list: List[SampleResult] = []

    print("[3/7] Reading FASTQs and matching target sequences", flush=True)
    for number, group in enumerate(groups, 1):
        print(f"  [{number}/{len(groups)}] {group.sample}", flush=True)
        result = process_group(
            group,
            matcher,
            args.left_primer,
            args.right_primer,
            args.anchor_mismatches,
            args.max_reads,
            file_metrics,
            header_indexes,
            observed_store,
            phix_detector,
        )
        sample_results_list.append(result)

    # I1/I2 files are not consumed by paired sample processing.
    paired_paths = {
        str(path)
        for group in groups
        for chunk in group.chunks
        for path in chunk
        if path is not None
    }
    index_paths = process_index_fastqs(
        parsed_files,
        file_metrics,
        args.max_reads,
        raw_index_reads,
    )
    for item in parsed_files:
        if str(item.path) not in paired_paths and str(item.path) not in index_paths:
            process_unpaired_qc_file(
                item,
                file_metrics[str(item.path)],
                args.max_reads,
                header_indexes,
            )

    if observed_store is not None:
        print("[4/7] Exporting observed target sequence table", flush=True)
        observed_store.close_and_export()
    else:
        print("[4/7] Observed-sequence export skipped", flush=True)

    print("[5/7] Writing sample and combined library metrics", flush=True)
    sample_metrics: Dict[str, Dict[str, object]] = {}
    sample_results: Dict[str, SampleResult] = {}
    for result in sample_results_list:
        metrics, _ = write_sample_outputs(args.outdir, result, reference)
        sample_metrics[result.sample] = metrics
        sample_results[result.sample] = result

    assigned = [x for x in sample_results_list if not x.is_undetermined]
    all_assigned = aggregate_results("ALL_ASSIGNED", assigned, reference)
    all_with_undetermined = aggregate_results(
        "ALL_WITH_UNDETERMINED", sample_results_list, reference
    )
    combined_dir = args.outdir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_metrics_rows = []
    for combined in (all_assigned, all_with_undetermined):
        metrics, rows = write_sample_outputs(args.outdir, combined, reference)
        sample_metrics[combined.sample] = metrics
        sample_results[combined.sample] = combined
        combined_metrics_rows.append(metrics)
        source_dir = args.outdir / "samples" / safe_name(combined.sample)
        for suffix in (
            "library_metrics.tsv",
            "variant_counts.tsv",
            "dropout_variants.tsv",
            "top_variants.tsv",
            "classification.tsv",
        ):
            source = source_dir / f"{safe_name(combined.sample)}_{suffix}"
            shutil.copy2(source, combined_dir / source.name)
    write_tsv(
        combined_dir / "combined_metrics.tsv",
        sorted({key for row in combined_metrics_rows for key in row}),
        combined_metrics_rows,
    )
    all_metric_rows = list(sample_metrics.values())
    write_tsv(
        combined_dir / "all_sample_metrics.tsv",
        sorted({key for row in all_metric_rows for key in row}),
        all_metric_rows,
    )
    write_cross_sample_matrices(
        args.outdir,
        reference,
        sample_results_list,
    )

    qc_rows = write_qc_outputs(
        args.outdir,
        sorted(file_metrics.values(), key=lambda x: (x.sample, x.read, x.path)),
    )
    sample_sheet_rows = parse_sample_sheet(args.sample_sheet)
    top_unknown = parse_top_unknown(args.top_unknown)

    print("[6/7] Running index QC", flush=True)
    index_summary = run_index_qc(
        args.outdir,
        sample_sheet_rows,
        raw_index_reads,
        header_indexes,
        top_unknown,
        args.max_index_distance,
        assigned_read_units=sum(x.total_units for x in assigned),
        undetermined_read_units=sum(
            x.total_units for x in sample_results_list if x.is_undetermined
        ),
        assigned_phix_read_units=sum(
            x.phix_reads for x in assigned
        ),
        undetermined_phix_read_units=sum(
            x.phix_reads for x in sample_results_list if x.is_undetermined
        ),
    )
    if not sample_sheet_rows:
        warnings.append(
            "No usable SampleSheet indexes were provided; expected-index "
            "substitution analysis is unavailable."
        )
    if not raw_index_reads.counter and not header_indexes.counter and not top_unknown:
        warnings.append(
            "No raw index strings were found in FASTQ headers and no Top Unknown "
            "Barcode file was provided."
        )

    print("[7/7] Creating HTML and shareable summary", flush=True)
    command = " ".join([Path(sys.argv[0]).name, *sys.argv[1:]])
    create_html_report(
        args.outdir,
        reference,
        qc_rows,
        sample_metrics,
        sample_results,
        index_summary,
        warnings,
        command,
        args.target_name,
    )
    create_share_summary(
        args.outdir,
        reference,
        qc_rows,
        sample_metrics,
        index_summary,
        sample_results,
        warnings,
        args,
    )

    manifest = {
        "pipeline": "amplicon_library_qc",
        "version": PIPELINE_VERSION,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "runtime_seconds": round(time.time() - started, 3),
        "python": sys.version,
        "reference_path": str(args.reference),
        "reference_sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        "input_dir": str(args.input_dir),
        "sample_names": [x.sample for x in sample_results_list],
        "arguments": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "warnings": warnings,
    }
    (args.outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        f"Done in {time.time() - started:.1f} seconds.\n"
        f"Share: {args.outdir / 'RESULTS_TO_SHARE.txt'}\n"
        f"Report: {args.outdir / 'report.html'}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
