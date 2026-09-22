"""Real symmetric K + lambda G = 0, with stable positive-definite K.

The elastic-restraint predictor changes the problem. Only full-matrix searches
produce original-system eigenpairs. Coordinates are scaled by congruence, never
regularized with artificial stiffness. Complex/unsymmetric problems are rejected.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigsh, splu


class InvalidModel(ValueError):
    """Inconsistent, unsupported, or numerically unsuitable input."""


@dataclass
class Member:
    name: str
    dofs: np.ndarray
    kg: np.ndarray | None = None
    ke: np.ndarray | None = None
    lateral: np.ndarray | None = None
    twist: np.ndarray | None = None
    elements: list[int] = field(default_factory=list)


@dataclass
class Config:
    modes_per_shift: int = 8
    max_shifts: int = 12
    max_seconds: float = 180.0  # Soft budget; cannot interrupt a native factorization.
    maxiter: int = 1500
    eigen_tol: float = 1e-9
    residual_tol: float = 1e-7
    cluster_ratio: float = 1.15
    participation_min: float = 0.15
    max_member_dofs: int = 160
    rhs_block: int = 16
    dense_limit: int = 768
    random_seed: int = 1729

    def validate(self) -> None:
        for name in ("modes_per_shift", "max_shifts", "maxiter", "max_member_dofs",
                     "rhs_block", "dense_limit"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise InvalidModel(f"{name} must be a positive integer")
        for name in ("max_seconds", "eigen_tol", "residual_tol", "cluster_ratio",
                     "participation_min"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise InvalidModel(f"{name} must be finite and positive")
        if self.cluster_ratio <= 1 or self.participation_min > 1:
            raise InvalidModel("Invalid clustering or participation threshold")
        if not isinstance(self.random_seed, int) or self.random_seed < 0:
            raise InvalidModel("random_seed must be a nonnegative integer")


def _real_array(value: Any, name: str) -> np.ndarray:
    if np.iscomplexobj(value):
        raise InvalidModel(f"{name}: complex data are unsupported")
    result = np.asarray(value, dtype=float)
    if not np.isfinite(result).all():
        raise InvalidModel(f"{name}: nonfinite data")
    return result


def symmetric_matrix(value: Any, name: str) -> sparse.csc_matrix:
    if np.iscomplexobj(value):
        raise InvalidModel(f"{name}: complex matrices are unsupported")
    a = sparse.csc_matrix(value, dtype=float)
    a.sum_duplicates()
    a.eliminate_zeros()
    if a.shape[0] != a.shape[1] or a.shape[0] < 2:
        raise InvalidModel(f"{name}: expected square matrix with at least two DOFs")
    if not np.isfinite(a.data).all():
        raise InvalidModel(f"{name}: nonfinite matrix entries")
    if sparse.linalg.norm(a - a.T) > 1e-10 * max(sparse.linalg.norm(a), 1e-300):
        raise InvalidModel(f"{name}: nonsymmetric; triangular storage must be declared at import")
    return a


def positive_roots(k: np.ndarray, g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Dense SPD-K reference: -G x = mu K x; lambda = 1/mu.

    This supports indefinite/singular G. It does NOT call eigh(K, -G), whose
    second matrix would incorrectly be required to be positive definite.
    Very small reciprocal roots are dropped relative to the spectral scale.
    """
    mu, vectors = linalg.eigh(-g, k, check_finite=True)
    scale = max(np.max(np.abs(mu)), np.finfo(float).tiny)
    indices = np.flatnonzero(mu > 1e-13 * scale)
    values = 1.0 / mu[indices]
    order = np.argsort(values)
    return values[order], vectors[:, indices[order]]


class Pencil:
    """Validated symmetric buckling pencil and member contributions.

    dofs are ZERO-BASED SOLVER-EQUATION indices, not node or element IDs.
    ke/kg are actual assembled member contributions in those coordinates.
    They must already include the SAME releases, transformations and constraint
    reductions as the global matrices; a global principal submatrix is not a
    substitute for the contribution of a selected set of elements.
    """

    def __init__(self, k: Any, g: Any, members: list[Member], config: Config | None = None):
        self.config = config or Config()
        self.config.validate()
        self.k, self.g = symmetric_matrix(k, "K"), symmetric_matrix(g, "G")
        if self.k.shape != self.g.shape:
            raise InvalidModel("K and G dimensions differ")
        self.n = self.k.shape[0]
        diagonal = self.k.diagonal()
        if np.any(diagonal <= 0):
            raise InvalidModel("K has nonpositive diagonal: review constraints/mechanisms")
        self.scale = 1.0 / np.sqrt(diagonal)
        d = sparse.diags(self.scale, format="csc")
        self.ks, self.gs = (d @ self.k @ d).tocsc(), (d @ self.g @ d).tocsc()
        try:
            if self.n <= self.config.dense_limit:
                linalg.cholesky(self.ks.toarray(), lower=True)
            else:
                start = np.random.default_rng(self.config.random_seed).normal(size=self.n)
                minimum = eigsh(self.ks, k=1, which="SA", v0=start,
                                maxiter=self.config.maxiter, tol=1e-9,
                                return_eigenvectors=False)[0]
                if minimum <= 1e-12:
                    raise InvalidModel("K is indefinite or numerically near singular")
        except (linalg.LinAlgError, ArpackNoConvergence) as exc:
            raise InvalidModel("Positive definiteness of K could not be established") from exc
        self.members = members
        self._elastic_factor = None
        names: set[str] = set()
        owned: set[int] = set()
        for member in members:
            if not member.name or member.name in names:
                raise InvalidModel("Member names must be nonempty and unique")
            names.add(member.name)
            raw = np.asarray(member.dofs)
            if raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer):
                raise InvalidModel(f"{member.name}: DOFs must be a 1D integer array")
            member.dofs = raw.astype(int)
            count = len(raw)
            if not count or len(set(raw)) != count or raw.min() < 0 or raw.max() >= self.n:
                raise InvalidModel(f"{member.name}: invalid or duplicate equation indices")
            if any(not isinstance(v, (int, np.integer)) or v <= 0 for v in member.elements):
                raise InvalidModel(f"{member.name}: element IDs must be positive integers")
            if len(set(member.elements)) != len(member.elements) or owned.intersection(member.elements):
                raise InvalidModel("Duplicate or overlapping physical-member element ownership")
            owned.update(member.elements)
            for name in ("kg", "ke"):
                value = getattr(member, name)
                if value is not None:
                    a = _real_array(value, f"{member.name}/{name}")
                    if a.shape != (count, count) or linalg.norm(a - a.T) > 1e-10 * max(linalg.norm(a), 1e-300):
                        raise InvalidModel(f"{member.name}: invalid or nonsymmetric {name}")
                    setattr(member, name, a)
                    if name == "ke":
                        minimum = linalg.eigvalsh(a, subset_by_index=[0, 0])[0]
                        if minimum < -1e-9 * max(linalg.norm(a), 1e-300):
                            raise InvalidModel(f"{member.name}: ke is not positive semidefinite")
            for name in ("lateral", "twist"):
                value = getattr(member, name)
                if value is not None:
                    a = _real_array(value, f"{member.name}/{name}")
                    if a.ndim != 2 or a.shape[1] != count or a.shape[0] == 0:
                        raise InvalidModel(f"{member.name}: invalid {name} observation operator")
                    setattr(member, name, a)

    def residual(self, value: float, x: np.ndarray) -> float:
        """Scaled relative residual, evaluated against the ORIGINAL full matrices."""
        x = np.asarray(x)
        if (np.iscomplexobj(x) or x.shape != (self.n,) or not np.isfinite(x).all()
                or not np.any(x) or not np.isfinite(value)):
            return float("inf")
        y = x / self.scale
        a, b = self.ks @ y, value * (self.gs @ y)
        return float(linalg.norm(a + b) / max(linalg.norm(a) + linalg.norm(b), 1e-300))

    def dense(self) -> tuple[np.ndarray, np.ndarray]:
        if self.n > self.config.dense_limit:
            raise InvalidModel("Dense solve exceeds configured size limit")
        values, vectors = positive_roots(self.ks.toarray(), self.gs.toarray())
        return values, self.scale[:, None] * vectors

    def predict(self, member: Member) -> tuple[float, np.ndarray]:
        """Elastic-restraint predictor, NOT a full-frame eigenpair.

        v0.1 retains all target support equations, not just end interfaces.
        K^-1 is never formed. A cached LU serves selected blocks of right-hand
        sides. Target contribution G_i is used; the remainder stays elastic.
        """
        ids = member.dofs
        if member.kg is None:
            raise InvalidModel("No assembled member geometric-stiffness contribution")
        if len(ids) > self.config.max_member_dofs:
            raise InvalidModel("Predictor DOF budget exceeded; provide a reviewed shift")
        if self._elastic_factor is None:
            self._elastic_factor = splu(self.ks)
        compliance = np.empty((len(ids), len(ids)))
        for first in range(0, len(ids), self.config.rhs_block):
            last = min(first + self.config.rhs_block, len(ids))
            rhs = np.zeros((self.n, last - first))
            rhs[ids[first:last], np.arange(last - first)] = 1.0
            compliance[:, first:last] = self._elastic_factor.solve(rhs)[ids, :]
        compliance = (compliance + compliance.T) * 0.5  # Roundoff only.
        effective = linalg.solve(compliance, np.eye(len(ids)), assume_a="pos")
        s = self.scale[ids]
        gi = s[:, None] * member.kg * s[None, :]
        values, vectors = positive_roots(effective, gi)
        if not len(values):
            raise InvalidModel("No positive predictor for this loading direction")
        start = np.zeros(self.n)
        start[ids] = vectors[:, 0]
        return float(values[0]), start

    def near(self, shift: float, start: np.ndarray | None = None
             ) -> tuple[np.ndarray, np.ndarray, dict]:
        if not np.isfinite(shift) or shift <= 0:
            raise InvalidModel("Buckling shift must be finite and strictly positive")
        cfg = self.config
        if self.n <= 3:
            values, vectors = self.dense()
            return values, vectors, {"solver": "dense-small", "partial": False}
        factor = splu(self.ks + shift * self.gs)
        inverse = LinearOperator(self.ks.shape, matvec=factor.solve, dtype=float)
        v0 = np.random.default_rng(cfg.random_seed).normal(size=self.n)
        if start is not None:
            start = _real_array(start, "starting vector")
            if start.shape != (self.n,):
                raise InvalidModel("Invalid starting-vector size")
            if linalg.norm(start) > 0:
                v0 = start / linalg.norm(start) + 1e-3 * v0 / linalg.norm(v0)
        count = min(cfg.modes_per_shift, self.n - 2)
        partial = False
        try:
            values, vectors = eigsh(self.ks, M=-self.gs, k=count, sigma=shift,
                                    which="LM", mode="buckling", OPinv=inverse,
                                    v0=v0, tol=cfg.eigen_tol, maxiter=cfg.maxiter,
                                    ncv=min(self.n, max(2 * count + 1, 20)))
        except ArpackNoConvergence as exc:
            values, vectors = exc.eigenvalues, exc.eigenvectors
            partial = True
        info = {"solver": "arpack-buckling", "partial": partial, "requested": count,
                "factor_nnz": int(factor.L.nnz + factor.U.nnz)}
        if values is None or vectors is None:
            return np.array([]), np.empty((self.n, 0)), info
        take = np.flatnonzero(np.isfinite(values) & (values > 0))
        take = take[np.argsort(values[take])]
        return values[take], self.scale[:, None] * vectors[:, take], info

    def native_parity(self, values: np.ndarray, vectors: np.ndarray) -> dict:
        values = _real_array(values, "native values")
        vectors = _real_array(vectors, "native vectors")
        if (values.ndim != 1 or not len(values) or vectors.shape != (self.n, len(values))
                or np.any(values == 0)):
            raise InvalidModel("Invalid native eigenpair arrays")
        errors = [self.residual(float(v), vectors[:, i]) for i, v in enumerate(values)]
        if max(errors) > self.config.residual_tol:
            raise InvalidModel("Native parity failed: review signs, mapping, and analysis stage")
        return {"algebraic_parity": True, "native_mode_count": len(values),
                "maximum_residual": max(errors), "physical_mapping_certified": False}
