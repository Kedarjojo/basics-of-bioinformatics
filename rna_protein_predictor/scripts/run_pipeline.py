#!/usr/bin/env python3
"""Restartable one-command orchestrator for the conference analysis."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = ROOT / "outputs" / "logs"
STATE_DIR = ROOT / "outputs" / "state"

STAGES = [
    ("download", [sys.executable, "scripts/download_data.py"]),
    ("rna_summary", [sys.executable, "-m", "src.build_rna_table"]),
    ("rna_symbols", [sys.executable, "-m", "src.add_gene_symbols_to_rna"]),
    ("gene_features", [sys.executable, "-m", "src.build_gene_features"]),
    ("hpa_target", [sys.executable, "-m", "src.build_master_table"]),
    ("augment", [sys.executable, "-m", "src.augment_master_table"]),
    ("localization", [sys.executable, "-m", "src.build_subcellular_features"]),
    ("audit", [sys.executable, "-m", "src.validate_dataset"]),
    ("loco_logistic", [sys.executable, "-m", "src.evaluate_loco", "--model", "logistic", "--bootstrap", "2000", "--seed", "7"]),
    ("loco_annotated", [sys.executable, "-m", "src.evaluate_loco", "--model", "logistic", "--bootstrap", "2000", "--seed", "7", "--annotated-only"]),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-stage", choices=[name for name, _ in STAGES])
    parser.add_argument("--until-stage", choices=[name for name, _ in STAGES])
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--resume", action="store_true",
        help="Skip stages with a successful completion marker from an earlier run.",
    )
    parser.add_argument("--bootstrap", type=int, help="Override LOCO bootstrap replicates")
    args = parser.parse_args()

    names = [name for name, _ in STAGES]
    start = names.index(args.from_stage) if args.from_stage else 0
    end = names.index(args.until_stage) + 1 if args.until_stage else len(STAGES)
    selected = STAGES[start:end]
    if args.skip_download:
        selected = [(name, cmd) for name, cmd in selected if name != "download"]

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest = {
        "run_id": run_id,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "stages": [],
    }
    for name, original in selected:
        marker = STATE_DIR / f"{name}.done.json"
        if args.resume and marker.exists():
            print(f"Skipping completed stage {name}: {marker.relative_to(ROOT)}")
            continue
        command = list(original)
        if args.bootstrap is not None and "--bootstrap" in command:
            command[command.index("--bootstrap") + 1] = str(args.bootstrap)
        print(f"\n=== {name}: {shlex.join(command)} ===", flush=True)
        if args.dry_run:
            continue
        log_path = LOG_DIR / f"{run_id}_{name}.log"
        started = datetime.now(timezone.utc)
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True
            )
        record = {
            "stage": name,
            "command": command,
            "started_utc": started.isoformat(),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "return_code": completed.returncode,
            "log": str(log_path.relative_to(ROOT)),
        }
        manifest["stages"].append(record)
        manifest_path = LOG_DIR / f"{run_id}_run.json"
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        if completed.returncode:
            print(f"Stage {name} failed. See {log_path}", file=sys.stderr)
            raise SystemExit(completed.returncode)
        marker.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print(f"Completed {name}; log: {log_path.relative_to(ROOT)}")
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    if not args.dry_run:
        (LOG_DIR / f"{run_id}_run.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )


if __name__ == "__main__":
    main()
