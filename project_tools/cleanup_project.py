#!/usr/bin/env python3
"""Safely consolidate the user's 5UTR NGS workspace.

The script is intentionally interactive: it discovers candidates and asks only
yes/no questions. Nothing is permanently deleted. Unselected and obsolete
items are moved to a timestamped quarantine beside 03_NGS.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_ROOT = Path("/data/user/MCET03/03_NGS")
PROJECT_NAME = "01_5UTR_Plasmid"
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
REFERENCE_SUFFIXES = {".csv", ".tsv", ".txt", ".fa", ".fasta", ".fna"}


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} TiB"


def path_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    total = 0
    for item in path.rglob("*"):
        if item.is_file():
            try:
                total += item.stat().st_size
            except FileNotFoundError:
                pass
    return total


def is_fastq(path: Path) -> bool:
    lower = path.name.lower()
    return any(lower.endswith(suffix) for suffix in FASTQ_SUFFIXES)


def is_inside(path: Path, parents: Iterable[Path]) -> bool:
    resolved = path.resolve()
    for parent in parents:
        try:
            resolved.relative_to(parent.resolve())
            return True
        except ValueError:
            continue
    return False


def is_pipeline_material(path: Path) -> bool:
    lowered = [part.lower() for part in path.parts]
    return any(
        "ngs_libraryqc" in part
        or "amplicon_library_qc_pipeline" in part
        or part in {"00_tools", "04_scripts", "tests", "test"}
        for part in lowered
    )


def discover_result_dirs(root: Path, package_dir: Path) -> List[Path]:
    results = set()
    for marker in root.rglob("ALL_ASSIGNED_variant_counts.tsv"):
        if is_inside(marker, [package_dir]) or is_pipeline_material(marker):
            continue
        if marker.parent.name == "combined":
            results.add(marker.parent.parent)
    return sorted(results)


def discover_fastq_groups(
    root: Path, package_dir: Path, result_dirs: Sequence[Path]
) -> List[Tuple[Path, List[Path]]]:
    groups: Dict[Path, List[Path]] = defaultdict(list)
    for path in root.rglob("*"):
        if not path.is_file() or not is_fastq(path):
            continue
        if is_inside(path, [package_dir, *result_dirs]) or is_pipeline_material(path):
            continue
        groups[path.parent].append(path)
    return sorted((parent, sorted(files)) for parent, files in groups.items())


def discover_references(
    root: Path, package_dir: Path, result_dirs: Sequence[Path]
) -> List[Path]:
    candidates = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in REFERENCE_SUFFIXES:
            continue
        lowered = path.name.lower()
        if "reference" not in lowered and "design" not in lowered:
            continue
        if is_inside(path, [package_dir, *result_dirs]) or is_pipeline_material(path):
            continue
        candidates.append(path)
    return sorted(candidates)


def discover_metadata(
    root: Path, package_dir: Path, result_dirs: Sequence[Path]
) -> List[Path]:
    candidates = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        compact = re.sub(r"[^a-z0-9]", "", path.name.lower())
        if "samplesheet" not in compact and "topunknownbarcodes" not in compact:
            continue
        if is_inside(path, [package_dir, *result_dirs]) or is_pipeline_material(path):
            continue
        candidates.append(path)
    return sorted(candidates)


def discover_figures(
    root: Path, package_dir: Path, result_dirs: Sequence[Path]
) -> List[Path]:
    candidates = []
    allowed = {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".html", ".r", ".rmd"}
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in allowed:
            continue
        lowered = {part.lower() for part in path.parts}
        if not ({"03_figures", "04_reports"} & lowered):
            continue
        if is_inside(path, [package_dir, *result_dirs]) or is_pipeline_material(path):
            continue
        candidates.append(path)
    return sorted(candidates)


def guess_fastq_role(path: Path) -> str:
    text = str(path).lower()
    if re.search(r"(?:1pct|1_percent|1percent|1%|subsample|subset)", text):
        return "subset"
    return "raw"


def ask(question: str) -> bool:
    while True:
        answer = input(f"{question} [yes/no]: ").strip().lower()
        if answer in {"yes", "y"}:
            return True
        if answer in {"no", "n"}:
            return False
        print("Please answer yes or no.")


def unique_destination(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    if path.name.lower().endswith(".fastq.gz"):
        stem = path.name[:-9]
        suffix = ".fastq.gz"
    elif path.name.lower().endswith(".fq.gz"):
        stem = path.name[:-6]
        suffix = ".fq.gz"
    number = 2
    while True:
        candidate = path.with_name(f"{stem}_{number}{suffix}")
        if not candidate.exists():
            return candidate
        number += 1


def normalized_date(path: Path) -> str:
    text = str(path)
    match = re.search(r"(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)", text)
    if match:
        return "".join(match.groups())
    match = re.search(r"(?<!\d)(\d{2})([01]\d)([0-3]\d)(?!\d)", text)
    if match:
        return "20" + "".join(match.groups())
    return "original_run"


def result_name(path: Path) -> str:
    date = normalized_date(path)
    text = str(path).lower()
    kind = "_1pct" if re.search(r"1pct|1_percent|1percent|1%", text) else ""
    version = ""
    manifest = path / "manifest.json"
    if manifest.is_file():
        try:
            raw_version = str(json.loads(manifest.read_text()).get("version", "")).strip()
            if raw_version:
                parts = raw_version.lstrip("v").split(".")
                version = "_v" + ".".join(parts[:2])
        except (OSError, ValueError, TypeError):
            pass
    if not version:
        match = re.search(r"v(\d+\.\d+)", text)
        if match:
            version = "_v" + match.group(1)
    return f"{date}{kind}{version}" if date != "original_run" else path.name


def print_inventory(
    fastq_groups: Sequence[Tuple[Path, Sequence[Path]]],
    references: Sequence[Path],
    metadata: Sequence[Path],
    results: Sequence[Path],
    figures: Sequence[Path],
) -> None:
    print("\nDetected FASTQ groups")
    if not fastq_groups:
        print("  NONE")
    for parent, files in fastq_groups:
        total = sum(path_size(path) for path in files)
        print(
            f"  {parent}\n"
            f"    {len(files)} FASTQ file(s), {human_size(total)}, "
            f"guess={guess_fastq_role(parent)}"
        )

    for title, items in (
        ("Reference candidates", references),
        ("SampleSheet / index metadata candidates", metadata),
        ("Essential analysis result candidates", results),
        ("Figure / report candidates", figures),
    ):
        print(f"\n{title}")
        if not items:
            print("  NONE")
        for item in items:
            print(f"  {item} ({human_size(path_size(item))})")


def move_file(source: Path, destination_dir: Path, records: List[dict]) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = unique_destination(destination_dir / source.name)
    shutil.move(str(source), str(destination))
    records.append({"source": str(source), "destination": str(destination)})


def ensure_not_running() -> None:
    result = subprocess.run(
        ["pgrep", "-f", "[a]mplicon_qc.py"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if result.returncode == 0:
        raise RuntimeError("amplicon_qc.py is running. Stop it before cleanup.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Discover and safely consolidate the 5UTR NGS project."
    )
    parser.add_argument("action", choices=["scan", "apply"], nargs="?", default="scan")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()

    root = args.root.resolve()
    package_dir = Path(__file__).resolve().parents[1]
    if not root.is_dir():
        parser.error(f"NGS root not found: {root}")
    if package_dir.parent.resolve() != root:
        parser.error(
            f"Extract the GitHub ZIP directly under {root} before running cleanup. "
            f"Current package: {package_dir}"
        )

    result_dirs = discover_result_dirs(root, package_dir)
    fastq_groups = discover_fastq_groups(root, package_dir, result_dirs)
    references = discover_references(root, package_dir, result_dirs)
    metadata = discover_metadata(root, package_dir, result_dirs)
    figures = discover_figures(root, package_dir, result_dirs)
    print_inventory(fastq_groups, references, metadata, result_dirs, figures)

    if args.action == "scan":
        print("\nScan only: no files were changed.")
        print(f"Next: bash {package_dir.name}/cleanup_project.sh apply")
        return 0

    ensure_not_running()
    selected_raw: List[Path] = []
    selected_subset: List[Path] = []
    for parent, files in fastq_groups:
        print(f"\nFASTQ group: {parent}")
        print(f"  {len(files)} file(s), {human_size(sum(path_size(x) for x in files))}")
        guess = guess_fastq_role(parent)
        if guess == "subset":
            if ask("Preserve this group as the extracted 1% dataset?"):
                selected_subset.extend(files)
            elif ask("Preserve this group as original raw data instead?"):
                selected_raw.extend(files)
        else:
            if ask("Preserve this group as original raw data?"):
                selected_raw.extend(files)
            elif ask("Preserve this group as the extracted 1% dataset instead?"):
                selected_subset.extend(files)

    selected_references = [
        path for path in references if ask(f"Preserve this reference?\n  {path}")
    ]
    selected_metadata = [
        path for path in metadata if ask(f"Preserve this metadata file?\n  {path}")
    ]
    selected_results = [
        path for path in result_dirs if ask(f"Preserve this analysis result?\n  {path}")
    ]
    selected_figures = [
        path for path in figures if ask(f"Preserve this figure/report?\n  {path}")
    ]

    if not selected_raw:
        print("\nNo original raw FASTQ group was selected. No changes were made.")
        return 2
    if not selected_references:
        print("\nNo reference was selected. No changes were made.")
        return 2
    if not ask(
        "Create the clean project and move every unselected item to recoverable quarantine?"
    ):
        print("Cancelled. No changes were made.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    staging_root = root.parent / f"{root.name}_staging_{stamp}"
    quarantine = root.parent / f"{root.name}_cleanup_quarantine_{stamp}"
    if staging_root.exists() or quarantine.exists():
        raise RuntimeError("Timestamped staging or quarantine path already exists.")

    project = staging_root / PROJECT_NAME
    for directory in (
        project / "raw_data",
        project / "reference",
        project / "subset_1pct",
        project / "results",
        project / "config",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    quarantine.mkdir(parents=True)
    moves: List[dict] = []

    # Move whole analysis roots before moving individual files.
    selected_result_roots = []
    for source in selected_results:
        destination = unique_destination(project / "results" / result_name(source))
        shutil.move(str(source), str(destination))
        moves.append({"source": str(source), "destination": str(destination)})
        selected_result_roots.append(source)

    for source in selected_raw:
        if source.exists():
            move_file(source, project / "raw_data" / normalized_date(source), moves)
    for source in selected_subset:
        if source.exists():
            move_file(source, project / "subset_1pct", moves)
    for source in selected_references:
        if source.exists():
            move_file(source, project / "reference", moves)
    for source in selected_metadata:
        if source.exists():
            move_file(source, project / "config", moves)
    for source in selected_figures:
        if source.exists():
            move_file(source, project / "results" / "figures", moves)

    quarantined = []
    for child in sorted(root.iterdir()):
        if child.resolve() == package_dir.resolve():
            continue
        destination = unique_destination(quarantine / child.name)
        shutil.move(str(child), str(destination))
        quarantined.append({"source": str(child), "destination": str(destination)})

    pipeline_destination = project / "pipeline"
    shutil.move(str(package_dir), str(pipeline_destination))
    moves.append({"source": str(package_dir), "destination": str(pipeline_destination)})
    final_project = root / PROJECT_NAME
    shutil.move(str(project), str(final_project))
    staging_root.rmdir()

    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project": str(final_project),
        "quarantine": str(quarantine),
        "moves": moves,
        "quarantined": quarantined,
    }
    manifest_path = final_project / "config" / "cleanup_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print("\nCleanup staging completed.")
    print(f"Clean project: {final_project}")
    print(f"Recoverable quarantine: {quarantine}")
    print("Nothing in quarantine has been deleted.")
    print("\nNext:")
    print(f"  bash {final_project}/pipeline/run_project.sh 1pct")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, shutil.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
