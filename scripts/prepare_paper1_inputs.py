#!/usr/bin/env python3
"""Verify pinned Paper 1 inputs and materialize only the handoff allowlist.

Uses the Python standard library. Never writes to Paper 1 or overwrites a
mismatching destination. Text source checks allow only CRLF/LF checkout changes;
archive and data checks are exact byte-level SHA-256 checks.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], stderr=subprocess.PIPE
    ).decode("utf-8").strip()


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def inside(root: Path, relative: str) -> Path:
    rel = PurePosixPath(relative)
    if rel.is_absolute() or ".." in rel.parts or ":" in relative or "\\" in relative:
        raise ValueError(f"Unsafe relative path: {relative}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes root: {relative}")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper1-root", type=Path, default=ROOT.parents[1])
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    lock_path = ROOT / "paper1_sources.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    source = args.paper1_root.resolve()
    repo = lock["paper1_repository"]
    if Path(git(source, "rev-parse", "--show-toplevel")).resolve() != source:
        raise ValueError("--paper1-root must be the Paper 1 Git root")
    commit = repo["commit"]
    if git(source, "rev-parse", f"{commit}^{{commit}}") != commit:
        raise ValueError("Pinned Paper 1 commit is unavailable")
    for relative, expected in repo["source_files"].items():
        path = inside(source, relative)
        if not path.is_file():
            raise ValueError(f"Missing source file: {relative}")
        pinned = git(source, "rev-parse", f"{commit}:{relative}")
        if pinned != expected["git_blob_oid"]:
            raise ValueError(f"Locked Git blob differs: {relative}")
        blob = subprocess.check_output(["git", "-C", str(source), "cat-file", "blob", pinned])
        if hashlib.sha256(blob).hexdigest() != expected["sha256_git_blob"]:
            raise ValueError(f"Source SHA-256 mismatch: {relative}")
        if path.read_bytes().replace(b"\r\n", b"\n") != blob.replace(b"\r\n", b"\n"):
            raise ValueError(f"Source differs from locked baseline: {relative}")
    release = lock["release"]
    archive = inside(source, release["archive_path"])
    if not archive.is_file() or digest(archive) != release["archive_sha256"]:
        raise ValueError("Archive missing, changed, or still an LFS pointer; retrieve the frozen LFS release")
    destination = inside(ROOT, f"data/external/paper1/{commit[:12]}/handoff")
    # Verify all data before materializing any of them. No extractall and no links.
    with tarfile.open(archive, "r:gz") as bundle:
        members = {member.name: member for member in bundle.getmembers()}
        checked = []
        for name, expected in release["handoff_files"].items():
            target = inside(destination, name)
            member_name = release["handoff_prefix"] + name
            member = members.get(member_name)
            if member is None or not member.isfile():
                raise ValueError(f"Missing or non-regular archive member: {member_name}")
            with bundle.extractfile(member) as handle:
                actual = hashlib.file_digest(handle, "sha256").hexdigest()
            if actual != expected:
                raise ValueError(f"Handoff checksum mismatch: {name}")
            if target.exists() and (not target.is_file() or digest(target) != expected):
                raise ValueError(f"Existing destination differs; refusing overwrite: {target}")
            checked.append((name, member, target, expected))
        if not args.verify_only:
            for name, member, target, expected in checked:
                if target.exists():
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                temporary = None
                try:
                    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as output:
                        temporary = Path(output.name)
                        with bundle.extractfile(member) as handle:
                            for block in iter(lambda: handle.read(1 << 20), b""):
                                output.write(block)
                    if digest(temporary) != expected:
                        raise ValueError(f"Written checksum mismatch: {name}")
                    # On Windows rename refuses an existing target. The exclusive
                    # creation below provides the same behavior on every platform.
                    with target.open("xb") as output, temporary.open("rb") as handle:
                        for block in iter(lambda: handle.read(1 << 20), b""):
                            output.write(block)
                finally:
                    if temporary is not None:
                        temporary.unlink(missing_ok=True)
            provenance = destination.parent / "provenance.json"
            evidence = {
                "prepared_utc": datetime.now(timezone.utc).isoformat(),
                "paper1_commit": commit,
                "source_root": str(source),
                "lock_sha256": digest(lock_path),
                "archive_sha256": release["archive_sha256"],
                "handoff_files": release["handoff_files"],
            }
            if not provenance.exists():
                with provenance.open("x", encoding="utf-8", newline="\n") as handle:
                    json.dump(evidence, handle, indent=2, ensure_ascii=False)
                    handle.write("\n")
    print(f"PASS: pinned commit, {len(repo['source_files'])} source files, archive and {len(checked)} handoff files")
    print("Verified only" if args.verify_only else f"Prepared: {destination}")
    print("Phase 0 simulation-field/transport suitability remains unaudited.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError, tarfile.TarError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
