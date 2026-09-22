"""Independent physical and synthetic verification fixtures, not ANSYS results."""
from __future__ import annotations

import numpy as np
from scipy import sparse

from .numerics import Member, Pencil


def uniform_mcr(e: float, iz: float, g: float, j: float, iw: float, length: float) -> float:
    """Prismatic doubly symmetric beam, uniform moment, fork ends, free warping."""
    values = np.array([e, iz, g, j, iw, length], dtype=float)
    if not np.isfinite(values).all() or np.any(values[[0, 1, 2, 3, 5]] <= 0) or iw < 0:
        raise ValueError("Invalid material, section properties or length")
    return float(np.pi / length * np.sqrt(e * iz * (g * j + np.pi**2 * e * iw / length**2)))


def _hermite(t: float, h: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = np.array([1-3*t*t+2*t**3, h*(t-2*t*t+t**3), 3*t*t-2*t**3, h*(-t*t+t**3)])
    d = np.array([-6*t+6*t*t, h*(1-4*t+3*t*t), 6*t-6*t*t, h*(-2*t+3*t*t)]) / h
    dd = np.array([-6+12*t, h*(-4+6*t), 6-12*t, h*(-2+6*t)]) / h**2
    return n, d, dd


def vlasov_beam(elements: int = 16, *, length: float = 6.0,
                reference_moment: float = 20000.0) -> tuple[Pencil, float]:
    """Two-field Hermite FE: v, v', phi, phi' at every station.

    Energy = .5 integral(EIz v''^2 + GJ phi'^2 + EIw phi''^2)
             + lambda integral(M0 phi v'').
    Excludes axial force, shear deformation and section distortion.
    SI units. v=phi=0 at both ends; slopes/warping are free.
    """
    if elements < 2 or length <= 0 or reference_moment <= 0:
        raise ValueError("Require >=2 elements and positive length/moment")
    e, iz, g, j, iw = 210e9, 8e-6, 210e9/2.6, 2e-7, 2e-7
    size, h = 4*(elements+1), length/elements
    k, kg = np.zeros((size, size)), np.zeros((size, size))
    points, weights = np.polynomial.legendre.leggauss(4)
    for element in range(elements):
        bending = np.array([4*element, 4*element+1, 4*element+4, 4*element+5])
        twisting = np.array([4*element+2, 4*element+3, 4*element+6, 4*element+7])
        for point, weight in zip(points, weights):
            n, d, dd = _hermite((point+1)/2, h)
            w = weight*h/2
            k[np.ix_(bending, bending)] += e*iz*np.outer(dd, dd)*w
            k[np.ix_(twisting, twisting)] += (g*j*np.outer(d, d)+e*iw*np.outer(dd, dd))*w
            coupling = reference_moment*np.outer(dd, n)*w
            kg[np.ix_(bending, twisting)] += coupling
            kg[np.ix_(twisting, bending)] += coupling.T
    free = np.setdiff1d(np.arange(size), [0, 2, 4*elements, 4*elements+2])
    k, kg = k[np.ix_(free, free)], kg[np.ix_(free, free)]
    identity = np.eye(len(free))
    member = Member("MEMBER_BEAM", np.arange(len(free)), kg, k,
                    identity[free % 4 == 0], identity[free % 4 == 2], list(range(1, elements+1)))
    return Pencil(k, kg, [member]), uniform_mcr(e, iz, g, j, iw, length)


def crowded_frame(count: int = 300, coupling: float = 0.002) -> Pencil:
    """Two-coordinate algebraic components; NOT a beam/frame FE validation.

    Positive roots span approximately 1..20; springs couple neighbors.
    G is symmetric indefinite. A component has both positive and negative roots.
    """
    if count < 3 or coupling < 0:
        raise ValueError("Require >=3 components and nonnegative coupling")
    targets = np.linspace(1, 20, count)
    k = sparse.eye(2*count, format="lil")
    kg = sparse.lil_matrix(k.shape)
    members = []
    for i, factor in enumerate(targets):
        ids = np.array([2*i, 2*i+1])
        gi = np.array([[0, -1/factor], [-1/factor, 0]])
        kg[np.ix_(ids, ids)] = gi
        members.append(Member(f"MEMBER_{i+1:04d}", ids, gi, np.eye(2),
                              np.array([[1., 0.]]), np.array([[0., 1.]]), [i+1]))
        if i:
            for offset in (0, 1):
                a, b = 2*i+offset, 2*(i-1)+offset
                k[a, a] += coupling
                k[b, b] += coupling
                k[a, b] -= coupling
                k[b, a] -= coupling
    return Pencil(k.tocsc(), kg.tocsc(), members)


def fixture_metadata(name: str) -> dict:
    return {"source": "synthetic", "load_case": name, "units": "fixture-defined; see benchmarks.py",
            "load_scaling": "proportional", "mapping": "fixture equation order, zero-based",
            "mapping_reviewed": False, "ansys_runtime_validated": False}
