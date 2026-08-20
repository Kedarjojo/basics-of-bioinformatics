#!/usr/bin/env python3
"""Download, extract and inventory every public input in config/downloads.tsv."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "config" / "downloads.tsv"


def sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    command = [
        "curl", "--fail", "--location", "--retry", "5",
        "--retry-delay", "5", "--continue-at", "-",
        "--output", str(partial), url,
    ]
    subprocess.run(command, check=True)
    partial.replace(destination)


def extract(archive: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".part")
    if archive.suffix == ".zip":
        with zipfile.ZipFile(archive) as bundle:
            members = [name for name in bundle.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise ValueError(f"Expected one file in {archive}, found {members}")
            with bundle.open(members[0]) as source, temporary.open("wb") as target:
                shutil.copyfileobj(source, target)
    elif archive.suffix == ".gz" and output.suffix != ".gz":
        with gzip.open(archive, "rb") as source, temporary.open("wb") as target:
            shutil.copyfileobj(source, target)
    else:
        if archive.resolve() == output.resolve():
            return
        shutil.copy2(archive, temporary)
    temporary.replace(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8", newline="") as handle:
        resources = list(csv.DictReader(handle, delimiter="\t"))
    inventory = []
    for item in resources:
        archive = ROOT / item["archive_path"]
        output = ROOT / item["output_path"]
        print(f"[{item['name']}] {item['url']} -> {output.relative_to(ROOT)}")
        if args.dry_run:
            continue
        if args.force or not archive.exists():
            download(item["url"], archive)
        else:
            print(f"  using existing archive: {archive.relative_to(ROOT)}")
        if args.force or not output.exists():
            extract(archive, output)
        else:
            print(f"  using existing output: {output.relative_to(ROOT)}")
        inventory.append({
            **item,
            "bytes": output.stat().st_size,
            "sha256": sha256(output),
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
        })

    if not args.dry_run:
        inventory_path = ROOT / "data" / "raw" / "download_inventory.tsv"
        inventory_path.parent.mkdir(parents=True, exist_ok=True)
        with inventory_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=inventory[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(inventory)
        print(f"Wrote {inventory_path.relative_to(ROOT)}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

