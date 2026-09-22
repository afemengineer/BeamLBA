"""Portable hash-checked equation-space bundles. Never unpickle solver data."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse
from scipy.io import mmread

from .numerics import Config, InvalidModel, Member, Pencil


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def safe_path(root: Path, value: str) -> Path:
    path = (root / value).resolve()
    if not path.is_relative_to(root.resolve()) or path == root.resolve():
        raise InvalidModel("Bundle path escapes its directory")
    return path


def write_bundle(directory: str | Path, k: Any, g: Any, members: list[Member],
                 metadata: dict, native: tuple[np.ndarray, np.ndarray] | None = None) -> Path:
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    sparse.save_npz(root / "K.npz", sparse.csc_matrix(k))
    sparse.save_npz(root / "G.npz", sparse.csc_matrix(g))
    arrays, rows = {}, []
    for i, member in enumerate(members):
        row = {"name": member.name, "dofs": member.dofs.tolist(), "elements": member.elements}
        for field in ("kg", "ke", "lateral", "twist"):
            value = getattr(member, field)
            if value is not None:
                key = f"m{i}_{field}"
                arrays[key], row[field] = value, key
        rows.append(row)
    if native is not None:
        arrays["native_values"], arrays["native_vectors"] = native
    np.savez_compressed(root / "operators.npz", **arrays)
    manifest = {"schema_version": 1, "equation": "K + lambda G = 0", "metadata": metadata,
                "members": rows, "sha256": {name: digest(root / name) for name in
                                             ("K.npz", "G.npz", "operators.npz")}}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False), encoding="utf-8")
    return root


def load_bundle(directory: str | Path, config: Config | None = None
                ) -> tuple[Pencil, dict, tuple[np.ndarray, np.ndarray] | None]:
    root = Path(directory)
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("equation") != "K + lambda G = 0":
        raise InvalidModel("Unsupported bundle schema or matrix equation")
    for name in ("K.npz", "G.npz", "operators.npz"):
        if manifest.get("sha256", {}).get(name) != digest(safe_path(root, name)):
            raise InvalidModel(f"Integrity check failed: {name}")
    metadata = manifest.get("metadata", {})
    for key in ("source", "load_case", "units", "load_scaling", "mapping"):
        if not metadata.get(key):
            raise InvalidModel(f"Missing provenance: {key}")
    if metadata["source"] not in ("synthetic", "ansys"):
        raise InvalidModel("Source must be synthetic or ansys")
    if metadata["load_scaling"] != "proportional":
        raise InvalidModel("v0.1 supports proportional loading only")
    members, native = [], None
    with np.load(root / "operators.npz", allow_pickle=False) as arrays:
        for row in manifest.get("members", []):
            kwargs = {field: arrays[row[field]] if field in row else None
                      for field in ("kg", "ke", "lateral", "twist")}
            members.append(Member(row["name"], np.asarray(row["dofs"]),
                                  elements=row.get("elements", []), **kwargs))
        has_values, has_vectors = "native_values" in arrays, "native_vectors" in arrays
        if has_values != has_vectors:
            raise InvalidModel("Native factors and vectors must be supplied together")
        if has_values:
            native = arrays["native_values"], arrays["native_vectors"]
    if not members:
        raise InvalidModel("Define at least one physical member")
    pencil = Pencil(sparse.load_npz(root / "K.npz"), sparse.load_npz(root / "G.npz"), members, config)
    return pencil, manifest, native


def read_matrix_market(path: str | Path, triangle: str = "none") -> sparse.csc_matrix:
    """Only explicitly declared triangular completion; never guess storage."""
    matrix = sparse.csc_matrix(mmread(path))
    if triangle not in ("none", "lower", "upper"):
        raise InvalidModel("Unknown triangular storage")
    if triangle != "none":
        opposite = sparse.triu(matrix, 1) if triangle == "lower" else sparse.tril(matrix, -1)
        if opposite.nnz:
            raise InvalidModel("Opposite triangle is populated; refusing to double entries")
        matrix = matrix + matrix.T - sparse.diags(matrix.diagonal())
    return matrix.tocsc()


def pack_matrix_market(k_path: str, g_path: str, specification: str, output: str,
                       g_sign: int = 1, triangle: str = "none") -> Path:
    """Import an explicit equation-space specification, never a guessed mapping."""
    if g_sign not in (-1, 1):
        raise InvalidModel("G sign must be +1 or -1")
    path = Path(specification)
    spec = json.loads(path.read_text(encoding="utf-8"))
    k, g = read_matrix_market(k_path, triangle), read_matrix_market(g_path, triangle) * g_sign
    arrays = {}
    if spec.get("operators"):
        with np.load(safe_path(path.parent, spec["operators"]), allow_pickle=False) as archive:
            arrays = {key: archive[key] for key in archive.files}
    members = []
    for row in spec["members"]:
        kwargs = {field: arrays[row[field]] if field in row else None
                  for field in ("kg", "ke", "lateral", "twist")}
        if kwargs["kg"] is not None:
            kwargs["kg"] = kwargs["kg"] * g_sign
        members.append(Member(row["name"], np.asarray(row["dofs"]), elements=row.get("elements", []), **kwargs))
    native = None
    if "native_values" in arrays and "native_vectors" in arrays:
        native = arrays["native_values"], arrays["native_vectors"]
    if "native_vectors_mtx" in spec:
        vectors = np.asarray(mmread(safe_path(path.parent, spec["native_vectors_mtx"])))
        native = np.asarray(spec["native_values"], dtype=float), vectors
    metadata = dict(spec["metadata"])
    metadata["import_g_sign"] = g_sign
    metadata["source_matrix_sha256"] = {"K": digest(Path(k_path)), "G": digest(Path(g_path))}
    pencil = Pencil(k, g, members)
    if native is not None:
        pencil.native_parity(*native)
    return write_bundle(output, k, g, members, metadata, native)
