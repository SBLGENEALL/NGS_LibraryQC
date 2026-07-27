#!/usr/bin/env python3
"""
Multiprocess launcher for amplicon_qc.py.

The original pipeline and output format are preserved. Read pairs are batched
in the parent process, classified by worker processes, and merged before the
original reporting steps continue.
"""

from __future__ import annotations

import argparse
import multiprocessing
import sys
import time
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import amplicon_qc as core


WORKER_STATE: Dict[str, object] = {}
# Use half of the visible logical CPUs, capped at 128.  This gives the user's
# 256-thread server 128 workers without oversubscribing smaller machines.
DEFAULT_WORKERS = max(1, min(128, multiprocessing.cpu_count() // 2))
DEFAULT_BATCH_SIZE = 5_000


def initialize_worker(
    matcher,
    left_primer: str,
    right_primer: str,
    anchor_mismatches: int,
    phix_detector,
    metric_metadata: Dict[str, Tuple[str, str, str, str]],
    collect_observations: bool,
) -> None:
    WORKER_STATE.clear()
    WORKER_STATE.update(
        {
            "matcher": matcher,
            "left_primer": left_primer,
            "right_primer": right_primer,
            "anchor_mismatches": anchor_mismatches,
            "phix_detector": phix_detector,
            "metric_metadata": metric_metadata,
            "collect_observations": collect_observations,
        }
    )


def new_metric(path: str):
    label, stored_path, sample, read = WORKER_STATE["metric_metadata"][path]
    return core.FastqMetrics(
        label=label,
        path=stored_path,
        sample=sample,
        read=read,
    )


def process_batch(payload):
    sample, r1_path, r2_path, records = payload
    matcher = WORKER_STATE["matcher"]
    left_primer = WORKER_STATE["left_primer"]
    right_primer = WORKER_STATE["right_primer"]
    anchor_mismatches = WORKER_STATE["anchor_mismatches"]
    phix_detector = WORKER_STATE["phix_detector"]

    n_refs = len(matcher.reference.unique)
    exact = [0] * n_refs
    near = [0] * n_refs
    categories = Counter()
    orientations = Counter()
    distances = Counter()
    header_indexes = Counter()
    observations = Counter()
    left_count = 0
    right_count = 0
    both_count = 0
    header_mismatch = 0
    pair_disagreements = 0
    phix_reads = 0

    metric_r1 = new_metric(r1_path)
    metric_r2 = new_metric(r2_path) if r2_path is not None else None

    for record1, record2 in records:
        h1, s1, q1 = record1
        metric_r1.update(s1, q1)
        index1 = core.extract_header_index(h1)
        if index1:
            header_indexes[index1] += 1

        h2 = s2 = q2 = None
        if record2 is not None:
            h2, s2, q2 = record2
            metric_r2.update(s2, q2)
            index2 = core.extract_header_index(h2)
            if index2 and index2 != index1:
                header_indexes[index2] += 1
            if core.canonical_read_id(h1) != core.canonical_read_id(h2):
                header_mismatch += 1

        if phix_detector is not None and phix_detector.is_phix(s1, s2):
            phix_reads += 1
            categories["phix"] += 1
            continue

        evidence = core.build_evidence(
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
        if WORKER_STATE["collect_observations"] and observed and "N" not in observed:
            observations[(sample, design_key, observed, match_type)] += 1

    result = core.SampleResult(
        sample=sample,
        is_undetermined=sample.lower().startswith("undetermined"),
        total_units=len(records),
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
    return result, r1_path, r2_path, metric_r1, metric_r2, header_indexes, observations


def merge_metric(destination, source) -> None:
    destination.total_reads += source.total_reads
    destination.total_bases += source.total_bases
    destination.q20_bases += source.q20_bases
    destination.q30_bases += source.q30_bases
    destination.qsum += source.qsum
    destination.n_bases += source.n_bases
    destination.reads_with_n += source.reads_with_n
    if source.min_length is not None:
        destination.min_length = (
            source.min_length
            if destination.min_length is None
            else min(destination.min_length, source.min_length)
        )
    destination.max_length = max(destination.max_length, source.max_length)

    while len(destination.cycles) < len(source.cycles):
        destination.cycles.append(core.CycleMetric())
    for dest_cycle, source_cycle in zip(destination.cycles, source.cycles):
        dest_cycle.bases.update(source_cycle.bases)
        dest_cycle.total += source_cycle.total
        dest_cycle.q20 += source_cycle.q20
        dest_cycle.q30 += source_cycle.q30
        dest_cycle.qsum += source_cycle.qsum


def empty_result(group, n_refs: int):
    return core.SampleResult(
        sample=group.sample,
        is_undetermined=group.is_undetermined,
        total_units=0,
        counts_exact=[0] * n_refs,
        counts_near=[0] * n_refs,
        categories=Counter(),
        anchor_left=0,
        anchor_right=0,
        anchors_both=0,
        pair_header_mismatch=0,
        orientation_counts=Counter(),
        distance_counts=Counter(),
        pair_disagreements=0,
        phix_reads=0,
    )


def merge_result(destination, source) -> None:
    destination.total_units += source.total_units
    destination.counts_exact = [
        a + b for a, b in zip(destination.counts_exact, source.counts_exact)
    ]
    destination.counts_near = [
        a + b for a, b in zip(destination.counts_near, source.counts_near)
    ]
    destination.categories.update(source.categories)
    destination.anchor_left += source.anchor_left
    destination.anchor_right += source.anchor_right
    destination.anchors_both += source.anchors_both
    destination.pair_header_mismatch += source.pair_header_mismatch
    destination.orientation_counts.update(source.orientation_counts)
    destination.distance_counts.update(source.distance_counts)
    destination.pair_disagreements += source.pair_disagreements
    destination.phix_reads += source.phix_reads


def iter_batches(
    group,
    batch_size: int,
    max_reads: Optional[int],
) -> Iterator[Tuple[str, str, Optional[str], list]]:
    submitted = 0
    for r1_path, r2_path in group.chunks:
        if r1_path is None and r2_path is None:
            continue
        if r1_path is None:
            r1_path, r2_path = r2_path, None
        assert r1_path is not None

        iterator_r1 = core.iter_fastq(r1_path)
        iterator_r2 = core.iter_fastq(r2_path) if r2_path is not None else None
        if iterator_r2 is None:
            iterator: Iterable = ((record, None) for record in iterator_r1)
        else:
            iterator = zip_longest(iterator_r1, iterator_r2)

        batch = []
        for record1, record2 in iterator:
            if record1 is None or (iterator_r2 is not None and record2 is None):
                raise ValueError(
                    f"R1/R2 read-count mismatch in sample {group.sample}: "
                    f"{r1_path} vs {r2_path}"
                )
            batch.append((record1, record2))
            submitted += 1
            if len(batch) >= batch_size:
                yield group.sample, str(r1_path), (
                    str(r2_path) if r2_path is not None else None
                ), batch
                batch = []
            if max_reads is not None and submitted >= max_reads:
                break
        if batch:
            yield group.sample, str(r1_path), (
                str(r2_path) if r2_path is not None else None
            ), batch
        if max_reads is not None and submitted >= max_reads:
            break


def make_parallel_process_group(workers: int, batch_size: int):
    def process_group_parallel(
        group,
        matcher,
        left_primer,
        right_primer,
        anchor_mismatches,
        max_reads,
        file_metrics,
        header_indexes,
        observed_store,
        phix_detector,
    ):
        n_refs = len(matcher.reference.unique)
        combined = empty_result(group, n_refs)
        metric_metadata = {
            path: (metric.label, metric.path, metric.sample, metric.read)
            for path, metric in file_metrics.items()
        }

        try:
            context = multiprocessing.get_context("fork")
        except ValueError:
            context = multiprocessing.get_context()

        started = time.time()
        last_report = started
        with context.Pool(
            processes=workers,
            initializer=initialize_worker,
            initargs=(
                matcher,
                left_primer,
                right_primer,
                anchor_mismatches,
                phix_detector,
                metric_metadata,
                observed_store is not None,
            ),
        ) as pool:
            batches = iter_batches(group, batch_size, max_reads)
            for output in pool.imap_unordered(process_batch, batches, chunksize=1):
                (
                    partial,
                    r1_path,
                    r2_path,
                    metric_r1,
                    metric_r2,
                    batch_indexes,
                    observations,
                ) = output
                merge_result(combined, partial)
                merge_metric(file_metrics[r1_path], metric_r1)
                if r2_path is not None and metric_r2 is not None:
                    merge_metric(file_metrics[r2_path], metric_r2)
                for key, count in batch_indexes.items():
                    header_indexes.update(key, count)
                if observed_store is not None and observations:
                    observed_store.buffer.update(observations)
                    if len(observed_store.buffer) >= observed_store.flush_size:
                        observed_store.flush()

                now = time.time()
                if now - last_report >= 30:
                    rate = combined.total_units / max(now - started, 0.001)
                    print(
                        f"    {combined.total_units:,} read pairs processed "
                        f"({rate:,.0f}/s; {workers} workers)",
                        flush=True,
                    )
                    last_report = now

        return combined

    return process_group_parallel


def parse_wrapper_args(argv: Sequence[str]):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    options, remaining = parser.parse_known_args(argv)
    if options.workers < 1:
        parser.error("--workers must be at least 1")
    if options.batch_size < 100:
        parser.error("--batch-size must be at least 100")
    return options, remaining


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(argv) if argv is not None else sys.argv[1:]
    options, core_argv = parse_wrapper_args(argv)
    core.process_group = make_parallel_process_group(
        workers=options.workers,
        batch_size=options.batch_size,
    )
    print(
        f"Parallel mode: {options.workers} workers, "
        f"batch size {options.batch_size:,}",
        flush=True,
    )
    return core.main(core_argv)


if __name__ == "__main__":
    raise SystemExit(main())
