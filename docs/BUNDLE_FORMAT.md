# Bundle and equation-space contract

A bundle directory contains `K.npz`, `G.npz`, `operators.npz` and `manifest.json`. Matrices use SciPy sparse NPZ; operators use non-object NumPy arrays loaded with `allow_pickle=False`. The manifest records schema version 1, the literal equation `K + lambda G = 0`, SHA-256 digests and provenance. Hashes detect accidental file mismatch; they are not signed attestations of physical correctness.

## Input specification for Matrix Market packaging

```json
{
  "metadata": {
    "source": "ansys",
    "load_case": "Wind case 07 / exact prestress step",
    "units": "N, m, rad; include warping convention",
    "load_scaling": "proportional",
    "mapping": "reviewed zero-based independent solver equations",
    "mapping_reviewed": false,
    "element_formulation_reviewed": false,
    "ansys_release": "2025 R2",
    "mechanical_snapshot_sha256": "copy from the matching inventory"
  },
  "operators": "member_operators.npz",
  "members": [
    {
      "name": "MEMBER_RAFTER_001",
      "dofs": [0, 1, 4, 7],
      "elements": [101, 102],
      "ke": "rafter_elastic",
      "kg": "rafter_geometric",
      "lateral": "rafter_lateral_observations",
      "twist": "rafter_twist_observations"
    }
  ],
  "native_vectors_mtx": "bl_native.mtx",
  "native_values": [2.1, 2.4, 3.2]
}
```

All numbers above are illustrative, not a usable ANSYS fixture. DOFs are zero-based solver-equation indices, NOT node IDs or element IDs. The local matrix row order must exactly match `dofs`.

For a member with p support equations, ke/kg are p-by-p matrices in the original physical equation coordinates. Observation operators have p columns and map those coordinates to reviewed local observations. Shared equations between members are valid; duplicate ownership of an element is not. ke and kg are optional, but missing kg disables automatic prediction and missing ke disables member attribution. Missing values are never silently filled with zeros.

Native vectors have n rows and one column per signed multiplier, in exactly the same solver order as K/G. Alternatively, put `native_values` and `native_vectors` in the operator NPZ and omit the Matrix Market native-vector keys.

```powershell
beamlba pack --k bl_k.mtx --g bl_g.mtx --spec spec.json --out bundle
beamlba verify bundle
beamlba run bundle --out results/case_07
```

`--g-sign -1` explicitly negates the imported global and member geometric matrices when needed to match the declared equation. It does not flip supplied native multipliers. Establish the correct convention by native parity, not by trying signs until a convenient positive answer appears.

`--triangle lower` or `upper` is only for a GENERAL-format file that truly stores one triangle. A Matrix Market SYMMETRIC header is expanded by the reader already; use the default `none`. Opposite-triangle entries cause declared triangular import to fail rather than double them.

## Review gate

A normal ANSYS run requires at least three native pairs where the dimension permits, valid residuals, and mapping/formulation review flags. All supplied pairs must pass. Diagnostic mode allows absent review but labels all member rows DIAGNOSTIC_ONLY. A failed native residual remains a hard error even in diagnostic mode.

The flags are recorded human/integration review decisions, not independently verified physical certificates. The CLI does not infer section properties, translations/rotations/warping meaning, physical load height, connection restraint or the governing load case from raw matrices.

## Search configuration

`examples/search_config.json` contains conservative starting caps. An optional `shifts` object maps exact member names to positive reviewed candidate multipliers. Shifts are locations to search, not accepted critical loads.

## Output

`report.json` records provenance, search shifts, numerical residuals, candidate member participation and unresolved reasons. `modes.npz` contains original-equation-space signed mode vectors and their factors. `report.html` is a human-readable summary. No output automatically maps new vectors back to Mechanical result objects in this release.
