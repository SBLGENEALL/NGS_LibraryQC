#!/usr/bin/env python3
"""Estimate full-run resources from a completed 1% run."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def tree_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def human_bytes(value: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} TiB"


def human_time(seconds: float) -> str:
    hours, remainder = divmod(int(round(seconds)), 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours} h {minutes} min {secs} sec"


def main() -> int:
    result = Path(sys.argv[1])
    subset = Path(sys.argv[2])
    full = Path(sys.argv[3])
    manifest = json.loads((result / "manifest.json").read_text(encoding="utf-8"))
    runtime = float(manifest["runtime_seconds"])
    subset_size = tree_size(subset)
    full_size = tree_size(full)
    result_size = tree_size(result)
    if subset_size <= 0 or full_size <= 0:
        raise RuntimeError("FASTQ size is zero; cannot estimate the full run.")
    scale = full_size / subset_size

    peak_rss = "not available"
    resource_file = result / "resource_usage.txt"
    if resource_file.is_file():
        match = re.search(
            r"Maximum resident set size \(kbytes\):\s*(\d+)",
            resource_file.read_text(encoding="utf-8", errors="replace"),
        )
        if match:
            peak_rss = human_bytes(int(match.group(1)) * 1024)

    text = (
        "Full-run estimate from the completed 1% run\n"
        f"Subset FASTQ size: {human_bytes(subset_size)}\n"
        f"Full FASTQ size: {human_bytes(full_size)}\n"
        f"Observed size ratio: {scale:.2f}x\n"
        f"Observed 1% runtime: {human_time(runtime)}\n"
        f"Estimated full runtime: {human_time(runtime * scale)}\n"
        f"Observed peak memory: {peak_rss}\n"
        f"Observed 1% result size: {human_bytes(result_size)}\n"
        f"Estimated full result size: {human_bytes(result_size * scale)}\n"
        "\nThe time and storage estimates assume approximately linear scaling. "
        "Peak memory is reported from the 1% run and should not be multiplied "
        "by the read-count ratio.\n"
    )
    output = result / "full_run_estimate.txt"
    output.write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, KeyError, IndexError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
