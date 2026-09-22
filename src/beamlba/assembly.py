"""Assemble REVIEWED element contributions after coordinate/constraint reduction.

The ANSYS adapter must supply each transformation. This module intentionally
never infers a transform from a node number, element number or Named Selection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .numerics import InvalidModel, Member, _real_array


@dataclass
class ElementContribution:
    element_id: int
    equations: np.ndarray
    transform: np.ndarray  # element coordinates = transform @ solver support coordinates
    ke: np.ndarray
    kg: np.ndarray

    def reduced(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ids = np.asarray(self.equations)
        if (ids.ndim != 1 or not np.issubdtype(ids.dtype, np.integer) or not len(ids)
                or np.any(ids < 0) or len(set(ids)) != len(ids)):
            raise InvalidModel("Invalid element solver-equation support")
        t = _real_array(self.transform, "element transform")
        k, g = _real_array(self.ke, "element K"), _real_array(self.kg, "element G")
        if (t.ndim != 2 or t.shape[1] != len(ids) or k.shape != (t.shape[0], t.shape[0])
                or g.shape != k.shape):
            raise InvalidModel("Element matrices/transform dimensions differ")
        return ids, t.T @ k @ t, t.T @ g @ t


def assemble_members(selections: dict[str, list[int]], elements: list[ElementContribution]) -> list[Member]:
    index = {item.element_id: item for item in elements}
    if len(index) != len(elements):
        raise InvalidModel("Duplicate element contribution")
    used, result = set(), []
    for name, element_ids in selections.items():
        if not element_ids or len(set(element_ids)) != len(element_ids) or used.intersection(element_ids):
            raise InvalidModel("Empty or overlapping physical-member ownership")
        if set(element_ids) - set(index):
            raise InvalidModel("Selected member has missing element matrices")
        used.update(element_ids)
        blocks = [index[element_id].reduced() for element_id in element_ids]
        support = np.unique(np.concatenate([block[0] for block in blocks]))
        k, g = np.zeros((len(support), len(support))), np.zeros((len(support), len(support)))
        for ids, ek, eg in blocks:
            locations = np.searchsorted(support, ids)
            k[np.ix_(locations, locations)] += ek
            g[np.ix_(locations, locations)] += eg
        result.append(Member(name, support, kg=g, ke=k, elements=list(element_ids)))
    return result


def assembly_residuals(k: sparse.spmatrix, g: sparse.spmatrix,
                       elements: list[ElementContribution]) -> dict:
    """Compare COMPLETE reconstructed assembly against full matrices, scaled by K.

    Supply non-beam/support/load-stiffness contributions too where applicable.
    Low residuals do not independently certify physical element orientation.
    """
    k, g = sparse.csc_matrix(k), sparse.csc_matrix(g)
    if k.shape != g.shape or np.any(k.diagonal() <= 0):
        raise InvalidModel("Invalid assembly reference matrices")
    rows, columns, kvals, gvals = [], [], [], []
    seen = set()
    for element in elements:
        if element.element_id in seen:
            raise InvalidModel("Duplicate element in full assembly")
        seen.add(element.element_id)
        ids, ek, eg = element.reduced()
        if np.max(ids) >= k.shape[0]:
            raise InvalidModel("Element equation exceeds global dimension")
        rows.extend(np.repeat(ids, len(ids)))
        columns.extend(np.tile(ids, len(ids)))
        kvals.extend(ek.ravel())
        gvals.extend(eg.ravel())
    rebuilt_k = sparse.coo_matrix((kvals, (rows, columns)), shape=k.shape).tocsc()
    rebuilt_g = sparse.coo_matrix((gvals, (rows, columns)), shape=g.shape).tocsc()
    d = sparse.diags(1 / np.sqrt(k.diagonal()))
    def error(a, b):
        return float(sparse.linalg.norm(d @ (a-b) @ d) / max(sparse.linalg.norm(d @ a @ d), 1e-300))
    return {"K_relative_residual": error(k, rebuilt_k), "G_relative_residual": error(g, rebuilt_g)}
