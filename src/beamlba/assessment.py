"""Budgeted member-informed searches; unresolved never means safe."""
from __future__ import annotations

from time import monotonic

import numpy as np
from scipy import linalg

from .numerics import InvalidModel, Pencil


def assess(pencil: Pencil, shifts: dict[str, float] | None = None,
           initial_modes: tuple[np.ndarray, np.ndarray] | None = None) -> tuple[dict, np.ndarray]:
    started = monotonic()
    cfg = pencil.config
    overrides = shifts or {}
    by_name = {m.name: m for m in pencil.members}
    if set(overrides) - set(by_name):
        raise InvalidModel("Shift overrides refer to unknown members")
    for value in overrides.values():
        if not np.isfinite(value) or value <= 0:
            raise InvalidModel("User shifts must be finite and positive")
    rows = {m.name: {"member": m.name, "status": "UNRESOLVED", "predictor": None,
                    "candidate_modes": [], "notes": [], "coverage_certified": False,
                    "design_check": "NOT_PERFORMED", "Mcr": None} for m in pencil.members}
    modes: list[dict] = []
    vectors: list[np.ndarray] = []
    searches: list[dict] = []

    def retain(values: np.ndarray, xs: np.ndarray) -> None:
        for value, x in zip(values, xs.T):
            if not np.isfinite(value) or value <= 0:
                continue
            error = pencil.residual(float(value), x)
            if error > cfg.residual_tol:
                continue
            energy = float(x @ (pencil.k @ x))
            if energy <= 0 or not np.isfinite(energy):
                continue
            x = x / np.sqrt(energy)
            # Do not merge independent vectors in a repeated eigenspace.
            duplicate = any(abs(float(value) - m["factor"]) <= 1e-7 * max(float(value), 1.0)
                            and abs(float(x @ (pencil.k @ y))) > 1 - 1e-6
                            for m, y in zip(modes, vectors))
            if not duplicate:
                modes.append({"factor": float(value), "residual": error})
                vectors.append(x)

    if initial_modes is not None:
        pencil.native_parity(*initial_modes)
        retain(*initial_modes)
    initial_count = len(modes)

    def represented(name: str) -> bool:
        member = by_name[name]
        return member.ke is not None and any(
            float(x[member.dofs] @ member.ke @ x[member.dofs]) >= cfg.participation_min
            for x in vectors)

    candidates = []
    for member in pencil.members:
        if represented(member.name):
            rows[member.name]["notes"].append("Candidate present in supplied native modes")
            continue
        if monotonic() - started > cfg.max_seconds:
            rows[member.name]["notes"].append("Soft time budget exhausted before predictor")
            continue
        try:
            prediction, seed = ((float(overrides[member.name]), None) if member.name in overrides
                                else pencil.predict(member))
            rows[member.name]["predictor"] = prediction
            # Store only local support values, not one full-frame vector per member.
            local_seed = None if seed is None else seed[member.dofs].copy()
            candidates.append((prediction, member.name, local_seed))
        except (InvalidModel, linalg.LinAlgError, RuntimeError) as exc:
            rows[member.name]["notes"].append(str(exc))
    clusters: list[list] = []
    for item in sorted(candidates, key=lambda item: item[0]):
        if not clusters or item[0] / clusters[-1][0][0] > cfg.cluster_ratio:
            clusters.append([item])
        else:
            clusters[-1].append(item)

    def search(group: list, phase: str) -> None:
        if len(searches) >= cfg.max_shifts or monotonic() - started > cfg.max_seconds:
            for _, name, _ in group:
                rows[name]["notes"].append("Search budget exhausted")
            return
        # A tiny detuning avoids factoring exactly at a predicted pole. A large
        # relative offset would miss targets in a densely populated spectrum.
        shift = float(np.exp(np.mean(np.log([p[0] for p in group]))) * (1 - 1e-5))
        seed = np.zeros(pencil.n)
        for _, name, local_seed in group:
            if local_seed is not None:
                seed[by_name[name].dofs] += local_seed
        entry = {"shift": shift, "members": [p[1] for p in group], "phase": phase}
        searches.append(entry)
        try:
            values, xs, info = pencil.near(shift, seed)
            entry.update(info)
            retain(values, xs)
        except (RuntimeError, ValueError, linalg.LinAlgError) as exc:
            entry["error"] = str(exc)

    for group in clusters:
        if not all(represented(item[1]) for item in group):
            search(group, "cluster")
    # Shared windows can miss edge members. Spend remaining budget on them.
    for item in candidates:
        if not represented(item[1]) and not any(
                len(entry["members"]) == 1 and entry["members"][0] == item[1] for entry in searches):
            search([item], "refinement")
    order = sorted(range(len(modes)), key=lambda j: modes[j]["factor"])
    modes, vectors = [modes[j] for j in order], [vectors[j] for j in order]
    for j, (mode, x) in enumerate(zip(modes, vectors)):
        mode["id"], mode["members"] = j + 1, []
        for member in pencil.members:
            if member.ke is None:
                continue
            local = x[member.dofs]
            participation = float(local @ member.ke @ local)  # x^T K x = 1.
            if participation > 1 + 1e-5 or participation < -1e-8:
                raise InvalidModel(f"{member.name}: ke is inconsistent with global K")
            if participation < cfg.participation_min:
                continue
            lateral = None if member.lateral is None else float(linalg.norm(member.lateral @ local))
            twist = None if member.twist is None else float(linalg.norm(member.twist @ local))
            observed = lateral is not None and twist is not None and lateral > 0 and twist > 0
            record = {"mode": j + 1, "factor": mode["factor"], "participation": participation,
                      "classification": "MEMBER_MODE_CANDIDATE",
                      "nonzero_lateral_and_twist_observations": bool(observed),
                      "lateral_observation_norm": lateral, "twist_observation_norm": twist}
            rows[member.name]["candidate_modes"].append(record)
            rows[member.name]["status"] = "CANDIDATE_REQUIRES_REVIEW"
            mode["members"].append(member.name)
        mode["coupled_candidate"] = len(mode["members"]) > 1
    for member in pencil.members:
        row = rows[member.name]
        if member.ke is None:
            row["notes"].append("No member elastic energy operator; attribution disabled")
        elif not row["candidate_modes"]:
            row["notes"].append("No candidate exceeded the attribution threshold; not a stability pass")
    report = {"schema_version": 1, "equation": "K + lambda G = 0", "dofs": pencil.n,
              "elapsed_seconds": monotonic() - started, "initial_mode_count": initial_count,
              "searches": searches, "modes": modes, "members": list(rows.values()),
              "coverage_certified": False, "design_check": "NOT_PERFORMED",
              "warnings": ["A small residual does not establish lowest-mode coverage.",
                           "Observation norms and energy attribution do not validate LTB classification.",
                           "Repeated-eigenspace member attribution can be basis dependent in this alpha.",
                           "Only the selected proportional load pattern is assessed.",
                           "Soft time limits cannot interrupt a running factorization."]}
    return report, np.column_stack(vectors) if vectors else np.empty((pencil.n, 0))
