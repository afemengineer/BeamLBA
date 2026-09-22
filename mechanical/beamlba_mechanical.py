# -*- coding: utf-8 -*-
"""Run in Mechanical's Scripting panel. Python 2.7/3-compatible launcher.

Default: read-only inventory. To run numerical assessment, provide a reviewed
matrix bundle and explicitly set RUN_ASSESSMENT=True. This alpha does NOT derive
solver-equation transformations from Named Selections or modify/solve the model.
See docs/ANSYS_INTEGRATION.md before using this on a project.
"""
from __future__ import print_function

import hashlib
import io
import json
import math
import os
import subprocess
import time
import uuid

STATIC_ANALYSIS_NAME = "Static Structural"
LBA_ANALYSIS_NAME = "Eigenvalue Buckling"
MEMBER_PREFIX = "MEMBER_"
OUTPUT_ROOT = os.path.join(os.path.expanduser("~"), "BeamLBA")
RUN_ASSESSMENT = False
PYTHON_EXE = r"C:\BeamLBA\.venv\Scripts\python.exe"
BUNDLE_DIRECTORY = None
SEARCH_CONFIG = None
MAX_WORKER_SECONDS = 300
MAX_WORKER_MEMORY_MB = 4096


def _text(value):
    try:
        return unicode(value)
    except NameError:
        return str(value)


def _read(obj, name, default=None):
    try:
        return getattr(obj, name)
    except Exception:
        return default


def _walk(obj):
    for child in list(_read(obj, "Children", []) or []):
        yield child
        for nested in _walk(child):
            yield nested


def select_analysis(analyses, name):
    matches = [item for item in analyses if _text(item.Name) == name]
    if len(matches) != 1:
        available = ", ".join(_text(item.Name) for item in analyses)
        raise ValueError("Choose one exact, unique analysis name. Available: " + available)
    return matches[0]


def member_elements(selection, mesh, geo_data):
    location = selection.Location
    kind = _text(location.SelectionType).split(".")[-1].lower()
    if kind == "meshelements":
        result = set(int(v) for v in location.Ids)
    elif kind == "geometryentities":
        result = set()
        for entity_id in location.Ids:
            entity = geo_data.GeoEntityById(int(entity_id))
            body_type = _text(_read(entity, "BodyType", "")).lower()
            if "wire" not in body_type and "line" not in body_type:
                raise ValueError("Use whole line-body selections or explicit mesh elements; not faces, edges, surfaces or solids")
            result.update(int(v) for v in mesh.MeshRegionById(int(entity_id)).ElementIds)
    else:
        raise ValueError("Unsupported scoping: " + kind + ". Nodal selections cannot define physical member ownership")
    if not result:
        raise ValueError("Empty or unmeshed member: " + _text(selection.Name))
    return sorted(result)


def collect(ext_api, static_name, lba_name, prefix):
    model = ext_api.DataModel.Project.Model
    analyses = list(model.Analyses)
    static = select_analysis(analyses, static_name)
    lba = select_analysis(analyses, lba_name)
    if static is lba:
        raise ValueError("Static and buckling references must differ")
    mesh = ext_api.DataModel.MeshDataByName("Global")
    members, ownership, all_nodes, element_rows = [], {}, set(), {}
    for selection in _walk(model.NamedSelections):
        if not _text(_read(selection, "Name", "")).startswith(prefix):
            continue
        if _read(selection, "Location") is None:
            continue
        name = _text(selection.Name)
        if any(row["name"] == name for row in members):
            raise ValueError("Duplicate member name: " + name)
        ids = member_elements(selection, mesh, ext_api.DataModel.GeoData)
        for element_id in ids:
            if element_id in ownership:
                raise ValueError("Element %s belongs to both %s and %s" % (element_id, ownership[element_id], name))
            ownership[element_id] = name
            element = mesh.ElementById(element_id)
            nodes = [int(v) for v in element.NodeIds]
            all_nodes.update(nodes)
            element_rows[element_id] = {"id": element_id, "node_ids": nodes,
                                        "mesh_type": _text(_read(element, "Type", "unavailable"))}
        members.append({"name": name, "named_selection_id": int(selection.ObjectId),
                        "element_ids": ids, "physical_member_identity_reviewed": False})
    if not members:
        raise ValueError("No nonempty member selections matching " + prefix)
    nodes = []
    for node_id in sorted(all_nodes):
        node = mesh.NodeById(node_id)
        nodes.append({"id": node_id, "coordinates": [float(node.X), float(node.Y), float(node.Z)]})
    working = _read(lba, "WorkingDir", _read(lba.Solution, "WorkingDir", ""))
    files = []
    if working and os.path.isdir(_text(working)):
        for filename in sorted(os.listdir(_text(working))):
            if filename.lower().endswith((".full", ".emat", ".mode", ".rst")):
                path = os.path.join(_text(working), filename)
                stat = os.stat(path)
                files.append({"path": path, "size": stat.st_size, "mtime": stat.st_mtime})
    references = {"static_name": _text(static.Name), "static_id": int(static.ObjectId),
                  "lba_name": _text(lba.Name), "lba_id": int(lba.ObjectId),
                  "static_type": _text(_read(static, "AnalysisType", "unavailable")),
                  "lba_type": _text(_read(lba, "AnalysisType", "unavailable"))}
    snapshot = {"references": references, "members": sorted(members, key=lambda row: row["name"]),
                "elements": [element_rows[key] for key in sorted(element_rows)],
                "nodes": nodes, "solver_files": files}
    encoded = json.dumps(snapshot, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    factors = []
    for result in _walk(lba.Solution):
        if _read(result, "IsSolved", False) is not True:
            continue
        value = _read(result, "LoadMultiplier")
        mode = _read(result, "Mode", _read(result, "SetNumber"))
        if value is not None:
            value = float(value)
            if not math.isnan(value) and not math.isinf(value):
                factors.append({"result_name": _text(result.Name), "factor": value,
                                "mode_or_set": None if mode is None else _text(mode)})
    return {"schema_version": 1, "snapshot": snapshot,
            "snapshot_sha256": hashlib.sha256(encoded).hexdigest(),
            "existing_result_factors": factors,
            "status": "INVENTORY_ONLY", "solver_equation_mapping_available": False,
            "warnings": ["Physical member identity and the static/prestress link require review.",
                         "The snapshot is not a complete hash of loads, materials or restraints.",
                         "Total-deformation magnitudes are not signed solver eigenvectors.",
                         "No project changes, solves or formulation changes were made."]}


def _write_json(path, data):
    with io.open(path, "w", encoding="utf-8") as stream:
        stream.write(_text(json.dumps(data, indent=2, ensure_ascii=True, allow_nan=False)))


def worker_arguments(bundle, output, config=None):
    args = ["-m", "beamlba.cli", "run", os.path.abspath(bundle), "--out", os.path.abspath(output)]
    if config:
        args.extend(["--config", os.path.abspath(config)])
    return args


def run_worker(python_exe, args, run_directory, timeout_seconds, memory_mb):
    """A separate numerical process permits a hard timeout without killing ANSYS."""
    if not os.path.isabs(python_exe) or not os.path.isfile(python_exe):
        raise ValueError("Set PYTHON_EXE to the absolute path of the installed CPython executable")
    if timeout_seconds <= 0 or memory_mb <= 0:
        raise ValueError("Worker limits must be positive")
    from System.Diagnostics import Process, ProcessStartInfo
    info = ProcessStartInfo()
    info.FileName = python_exe
    info.Arguments = subprocess.list2cmdline(args)
    info.WorkingDirectory = run_directory
    info.UseShellExecute = False
    info.CreateNoWindow = True
    info.RedirectStandardOutput = True
    info.RedirectStandardError = True
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        info.EnvironmentVariables[name] = "1"
    process = Process()
    process.StartInfo = info
    process.Start()
    stdout = process.StandardOutput.ReadToEndAsync()
    stderr = process.StandardError.ReadToEndAsync()
    started, failure = time.time(), None
    try:
        while not process.WaitForExit(250):
            process.Refresh()
            if time.time() - started > timeout_seconds:
                failure = "Worker time budget exceeded; no approval issued"
                break
            if process.PrivateMemorySize64 > memory_mb * 1024 * 1024:
                failure = "Worker memory budget exceeded; no approval issued"
                break
    finally:
        if not process.HasExited:
            process.Kill()
        process.WaitForExit()
        with io.open(os.path.join(run_directory, "worker_stdout.log"), "w", encoding="utf-8") as stream:
            stream.write(_text(stdout.Result))
        with io.open(os.path.join(run_directory, "worker_stderr.log"), "w", encoding="utf-8") as stream:
            stream.write(_text(stderr.Result))
    if failure:
        raise RuntimeError(failure)
    if process.ExitCode != 0:
        raise RuntimeError("Numerical worker failed; inspect worker_stderr.log")


def main(ext_api):
    inventory = collect(ext_api, STATIC_ANALYSIS_NAME, LBA_ANALYSIS_NAME, MEMBER_PREFIX)
    run_directory = os.path.join(OUTPUT_ROOT, time.strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8])
    os.makedirs(run_directory)
    _write_json(os.path.join(run_directory, "inventory.json"), inventory)
    if not RUN_ASSESSMENT:
        ext_api.Log.WriteMessage("BeamLBA inventory saved: " + run_directory + ". No buckling assessment was run.")
        return run_directory
    try:
        if not BUNDLE_DIRECTORY:
            raise ValueError("A reviewed matrix bundle is required. See docs/ANSYS_INTEGRATION.md")
        with io.open(os.path.join(BUNDLE_DIRECTORY, "manifest.json"), "r", encoding="utf-8") as stream:
            manifest = json.load(stream)
        metadata = manifest.get("metadata", {})
        if metadata.get("source") != "ansys":
            raise ValueError("Mechanical assessments require an ANSYS-source bundle")
        if metadata.get("mechanical_snapshot_sha256") != inventory["snapshot_sha256"]:
            raise ValueError("Bundle snapshot differs from the current model; regenerate and review the export")
        result_directory = os.path.join(run_directory, "assessment")
        args = worker_arguments(BUNDLE_DIRECTORY, result_directory, SEARCH_CONFIG)
        run_worker(PYTHON_EXE, args, run_directory, MAX_WORKER_SECONDS, MAX_WORKER_MEMORY_MB)
        ext_api.Log.WriteMessage("BeamLBA candidate report: " + os.path.join(result_directory, "report.html"))
    except Exception as exc:
        _write_json(os.path.join(run_directory, "status.json"),
                    {"status": "NOT_ASSESSED", "reason": _text(exc), "design_check": "NOT_PERFORMED"})
        ext_api.Log.WriteError("BeamLBA: " + _text(exc))
        raise
    return run_directory


if "ExtAPI" in globals():
    main(ExtAPI)
