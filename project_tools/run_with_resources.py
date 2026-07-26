#!/usr/bin/env python3
"""Run one command while streaming logs and recording Linux resource usage."""

from __future__ import annotations

import argparse
import resource
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resource-file", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    started = time.monotonic()
    with args.log_file.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        status = process.wait()

    elapsed = time.monotonic() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    args.resource_file.write_text(
        "Resource measurement: Python resource.getrusage\n"
        f"Elapsed wall clock time (seconds): {elapsed:.3f}\n"
        f"User time (seconds): {usage.ru_utime:.3f}\n"
        f"System time (seconds): {usage.ru_stime:.3f}\n"
        f"Maximum resident set size (kbytes): {usage.ru_maxrss}\n"
        f"Exit status: {status}\n",
        encoding="utf-8",
    )
    return status


if __name__ == "__main__":
    raise SystemExit(main())
