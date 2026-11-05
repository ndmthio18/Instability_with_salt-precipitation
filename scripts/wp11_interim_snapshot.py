"""Create an immutable, read-only snapshot of an active WP11 run.

The script never writes below ``run_dir``.  It freezes byte prefixes of the
three append-only diagnostic CSVs, selects their largest common consecutive
accepted step, validates a small set of stable VTU files already registered in
the PVD, builds the snapshot on node-local temporary storage, and only then
copies one archive plus one record to durable storage.
"""
from __future__ import print_function

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import tempfile
import time
import traceback
import xml.etree.ElementTree as ET


TARGET_T = (2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 12.0, 15.0, 18.0)
FRAME_STRIDE = 50
STEPS_PER_T = 1200.0


def utc_now():
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixed_complete_prefix(path):
    size = os.stat(path).st_size
    with open(path, "rb") as handle:
        data = handle.read(size)
    if len(data) != size:
        raise RuntimeError("short fixed-prefix read: %s" % path)
    end = data.rfind(b"\n")
    if end < 0:
        raise RuntimeError("no complete line: %s" % path)
    return data[:end + 1], {"source_size_at_freeze": size, "complete_prefix_bytes": end + 1}


def csv_rows(data):
    text = data.decode("utf-8")
    reader = csv.DictReader(io.StringIO(text))
    rows = list(reader)
    if not reader.fieldnames:
        raise RuntimeError("missing CSV header")
    return reader.fieldnames, rows


def maximum_consecutive_single(rows):
    mapping = {}
    for row in rows:
        step = int(row["step"])
        if step in mapping:
            raise RuntimeError("duplicate step %d" % step)
        mapping[step] = row
    step = 0
    while step in mapping:
        step += 1
    return step - 1, mapping


def maximum_consecutive_budget(rows, expected_regions):
    groups = {}
    for row in rows:
        groups.setdefault(int(row["step"]), []).append(row)
    step = 0
    expected = set(expected_regions)
    while step in groups:
        regions = [row["region"] for row in groups[step]]
        if len(regions) != len(expected_regions) or set(regions) != expected:
            break
        step += 1
    return step - 1, groups


def write_rows(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_verified(source, destination):
    shutil.copyfile(source, destination)
    source_hash = sha256(source)
    destination_hash = sha256(destination)
    if source_hash != destination_hash:
        raise RuntimeError("copy hash mismatch: %s" % source)
    return {"bytes": os.stat(destination).st_size, "sha256": destination_hash}


def tail_complete_lines(path, count=200, maximum_bytes=1024 * 1024):
    size = os.stat(path).st_size
    start = max(0, size - maximum_bytes)
    with open(path, "rb") as handle:
        handle.seek(start)
        data = handle.read(size - start)
    if start:
        newline = data.find(b"\n")
        data = data[newline + 1:] if newline >= 0 else b""
    end = data.rfind(b"\n")
    if end >= 0:
        data = data[:end + 1]
    return b"".join(data.splitlines(True)[-count:])


def pvd_records(data):
    pattern = re.compile(r'<DataSet\s+timestep="([^"]+)".*?file="([^"]+\.vtu)"')
    records = []
    for match in pattern.finditer(data.decode("utf-8", "replace")):
        filename = match.group(2)
        frame_match = re.search(r"-(\d+)\.vtu$", filename)
        if not frame_match:
            continue
        frame = int(frame_match.group(1))
        records.append({
            "frame": frame,
            "step": frame * FRAME_STRIDE,
            "pvd_time_s": float(match.group(1)),
            "pvd_time_text": match.group(1),
            "filename": filename,
        })
    return records


def select_field_records(records, snapshot_step):
    available = [record for record in records if record["step"] <= snapshot_step]
    if not available:
        raise RuntimeError("no registered VTU at or before snapshot step")
    selected = {}
    for target in TARGET_T:
        record = min(available, key=lambda item: abs(item["step"] / STEPS_PER_T - target))
        selected[record["filename"]] = dict(record, requested_time_over_T=target)
    last = max(available, key=lambda item: item["step"])
    selected[last["filename"]] = dict(last, requested_time_over_T="latest_stable")
    return [selected[name] for name in sorted(selected, key=lambda item: selected[item]["frame"])]


def validate_stable_fields(run_dir, records):
    before = {}
    for record in records:
        path = os.path.join(run_dir, record["filename"])
        stat = os.stat(path)
        before[path] = (stat.st_size, stat.st_mtime_ns)
    time.sleep(3)
    for record in records:
        path = os.path.join(run_dir, record["filename"])
        stat = os.stat(path)
        if before[path] != (stat.st_size, stat.st_mtime_ns):
            raise RuntimeError("VTU changed during stability window: %s" % path)
        if stat.st_size <= 0:
            raise RuntimeError("empty VTU: %s" % path)
        ET.parse(path)


def make_partial_pvd(path, records):
    with open(path, "w") as handle:
        handle.write('<?xml version="1.0"?>\n')
        handle.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        handle.write('<Collection>\n')
        for record in records:
            handle.write('  <DataSet timestep="%.17g" group="" part="0" name="" file="%s"/>\n' %
                         (record["actual_time_s"], record["filename"]))
        handle.write('</Collection>\n</VTKFile>\n')


def file_manifest(root):
    entries = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if os.path.isfile(path) and name != "manifest.json":
            entries.append({"path": name, "bytes": os.stat(path).st_size, "sha256": sha256(path)})
    return entries


def capture(args):
    run_dir = os.path.abspath(args.run_dir)
    temporary_parent = args.temporary_parent or os.path.dirname(run_dir)
    temporary_root = tempfile.mkdtemp(prefix="wp11_interim_", dir=temporary_parent)
    snapshot_dir = os.path.join(temporary_root, args.label)
    os.mkdir(snapshot_dir)
    try:
        metadata_path = os.path.join(run_dir, "transport_metadata.json")
        with open(metadata_path) as handle:
            metadata = json.load(handle)
        expected_regions = [item["name"] for item in metadata["cuts"]]

        frozen = {}
        parsed = {}
        for name in ("surface_events.csv", "weighted_online.csv", "transport_budget.csv"):
            data, audit = fixed_complete_prefix(os.path.join(run_dir, name))
            fields, rows = csv_rows(data)
            frozen[name] = audit
            parsed[name] = (fields, rows)

        surface_max, surface_map = maximum_consecutive_single(parsed["surface_events.csv"][1])
        weighted_max, weighted_map = maximum_consecutive_single(parsed["weighted_online.csv"][1])
        budget_max, budget_groups = maximum_consecutive_budget(parsed["transport_budget.csv"][1], expected_regions)
        snapshot_step = min(surface_max, weighted_max, budget_max)
        if snapshot_step < 0:
            raise RuntimeError("no common complete accepted step")
        snapshot_time = float(surface_map[snapshot_step]["time_s"])

        maximum_time_error = 0.0
        for step in range(snapshot_step + 1):
            reference = float(surface_map[step]["time_s"])
            maximum_time_error = max(maximum_time_error,
                                     abs(float(weighted_map[step]["time_s"]) - reference))
            for row in budget_groups[step]:
                maximum_time_error = max(maximum_time_error, abs(float(row["time_s"]) - reference))
        if maximum_time_error > 1e-8:
            raise RuntimeError("diagnostic timestamp mismatch %.17g s" % maximum_time_error)

        filtered = {
            "surface_events.csv": [surface_map[step] for step in range(snapshot_step + 1)],
            "weighted_online.csv": [weighted_map[step] for step in range(snapshot_step + 1)],
            "transport_budget.csv": [row for step in range(snapshot_step + 1) for row in budget_groups[step]],
        }
        for name, rows in filtered.items():
            write_rows(os.path.join(snapshot_dir, name), parsed[name][0], rows)

        copied = {}
        for name in ("transport_metadata.json", "params.input", "eigenprofile.dat"):
            copied[name] = copy_verified(os.path.join(run_dir, name), os.path.join(snapshot_dir, name))
        solver_log = os.path.join(run_dir, "solver.log")
        if os.path.exists(solver_log):
            with open(os.path.join(snapshot_dir, "solver_tail.log"), "wb") as handle:
                handle.write(tail_complete_lines(solver_log))

        pvd_name = args.run_name + ".pvd"
        pvd_data, pvd_audit = fixed_complete_prefix(os.path.join(run_dir, pvd_name))
        registered = pvd_records(pvd_data)
        selected = select_field_records(registered, snapshot_step)
        validate_stable_fields(run_dir, selected)
        for record in selected:
            step = record["step"]
            record["actual_time_s"] = float(surface_map[step]["time_s"])
            record["actual_time_over_T"] = step / STEPS_PER_T
            record["pvd_time_error_s"] = record["pvd_time_s"] - record["actual_time_s"]
            record["source"] = os.path.join(run_dir, record["filename"])
            record.update(copy_verified(record["source"], os.path.join(snapshot_dir, record["filename"])))
        make_partial_pvd(os.path.join(snapshot_dir, args.label + ".pvd"), selected)

        status = {
            "schema_version": 1,
            "label": args.label,
            "captured_utc": utc_now(),
            "production_job_id": args.production_job_id,
            "production_state_at_capture": "RUNNING",
            "production_solver_untouched": True,
            "source_run_dir": run_dir,
            "snapshot_end_step": snapshot_step,
            "snapshot_end_time_s": snapshot_time,
            "snapshot_end_time_over_T": snapshot_step / STEPS_PER_T,
            "maximum_common_timestamp_error_s": maximum_time_error,
            "diagnostic_source_max_steps": {
                "surface_events.csv": surface_max,
                "weighted_online.csv": weighted_max,
                "transport_budget.csv": budget_max,
            },
            "diagnostic_prefix_audit": frozen,
            "pvd_prefix_audit": pvd_audit,
            "selected_fields": selected,
            "provenance": {
                "executable_sha256": args.executable_sha256,
                "input_sha256": copied["params.input"]["sha256"],
                "profile_sha256": copied["eigenprofile.dat"]["sha256"],
            },
            "scope": "read-only interim snapshot; not a terminal simulation state",
        }
        with open(os.path.join(snapshot_dir, "source_status.json"), "w") as handle:
            json.dump(status, handle, indent=2, sort_keys=True)

        manifest = dict(status)
        manifest["files"] = file_manifest(snapshot_dir)
        manifest_path = os.path.join(snapshot_dir, "manifest.json")
        with open(manifest_path, "w") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)

        archive_path = os.path.join(temporary_root, args.label + ".tar.gz")
        with tarfile.open(archive_path, "w:gz") as archive:
            archive.add(snapshot_dir, arcname=args.label)
        record = {
            "schema_version": 1,
            "label": args.label,
            "created_utc": utc_now(),
            "production_job_id": args.production_job_id,
            "snapshot_end_step": snapshot_step,
            "snapshot_end_time_s": snapshot_time,
            "snapshot_end_time_over_T": snapshot_step / STEPS_PER_T,
            "archive_bytes": os.stat(archive_path).st_size,
            "archive_sha256": sha256(archive_path),
            "manifest_sha256": sha256(manifest_path),
            "complete": True,
        }
        with open(args.durable_archive, "wb") as destination, open(archive_path, "rb") as source:
            shutil.copyfileobj(source, destination, 1024 * 1024)
            destination.flush()
            os.fsync(destination.fileno())
        if sha256(args.durable_archive) != record["archive_sha256"]:
            raise RuntimeError("durable archive hash mismatch")
        with open(args.durable_record, "w") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        print(json.dumps(record, indent=2, sort_keys=True))
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    parser.add_argument("durable_archive")
    parser.add_argument("durable_record")
    parser.add_argument("--label", default="WP11-R96-interim-18T")
    parser.add_argument("--run-name", default="WP11-R96-M1-long25T")
    parser.add_argument("--production-job-id", required=True)
    parser.add_argument("--executable-sha256", required=True)
    parser.add_argument("--temporary-parent")
    args = parser.parse_args()
    try:
        capture(args)
    except Exception:
        error = {"schema_version": 1, "label": args.label, "complete": False,
                 "failed_utc": utc_now(), "error": traceback.format_exc()}
        with open(args.durable_record, "w") as handle:
            json.dump(error, handle, indent=2, sort_keys=True)
        raise


if __name__ == "__main__":
    main()
