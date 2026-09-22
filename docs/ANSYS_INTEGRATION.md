# ANSYS integration: current procedure and remaining gate

## Read this first

The intended final experience is one Mechanical action: select a static case, its LBA and member Named Selections, then assess member candidates. The alpha implements the numerical backend, a read-only Mechanical inventory, a bounded external-worker launcher, matrix import, and element-assembly primitives. It does **not** yet automatically translate every ANSYS beam/connection formulation into correctly reduced member contributions. That is an explicit release gate, not a hidden fallback.

No ANSYS executable or license was available during implementation. Mechanical calls have been tested with mock objects only. The MAPDL export helper has been checked against documentation but not executed in ANSYS. In particular, Python/.NET process-launch behaviour must be checked on the installation being deployed.

## Supported numerical scope

Real symmetric matrices, stable positive-definite elastic K, proportional reference loading, and elastic bifurcation. Start integration with a simple prismatic BEAM188 model with reviewed warping behaviour, then add BEAM189, offsets, releases and connections through separate fixture tests. Presence of line bodies is not evidence that an arbitrary model is supported.

An unsymmetric pencil, nonlinear-base perturbation interpretation, follower-load case, contact-state change, mechanism, or missing warping behaviour must not be silently converted into this scope. A load case designated 'worst' is assessed only for its own load pattern.

## Phase 1: read-only inventory

Open `mechanical/beamlba_mechanical.py` in Mechanical's Scripting panel. Set the exact analysis names. The names must be unique. Set the member prefix if necessary; default `MEMBER_`.

The script accepts whole line-body selections and explicit mesh-element selections. It rejects nodal selections and surface/face scoping. A selection must represent one physical member; automated splitting of a family selection is not implemented. The analyst must review member identity and whether each selected mesh element is actually a supported beam.

Leave `RUN_ASSESSMENT=False`. The script records selected elements, mesh connectivity, coordinates, analysis identifiers, available solver-file paths/sizes/timestamps and available solved-result multipliers. It writes a new run directory under `~/BeamLBA`. It does not solve, change KEYOPTs, add constraints or modify the project.

The snapshot hash detects changes in that inventory. It is NOT a complete fingerprint of materials, loads, connection details or all model settings. A matching snapshot alone does not prove model equivalence.

## Phase 2: preserve and inspect the actual buckling data

Preserve matching `.full`, `.mode`, `.emat` (where needed), solver input/output and reference load-step information. Files must belong to the same completed buckling case/stage. Copy them into a separate scratch directory; never experiment in the live project solver directory.

The helper `mapdl/export_existing.mac` expects short filenames `file.full` and `file.mode` in the scratch directory. It imports the assembled matrices and converts internal-order MODE vectors to solver order using the documented NOD2SOLV multiplication. It exports Matrix Market files. It deliberately does not call SOLVE, WRFULL or RESUME.

Important distinctions:

- ANSYS documents MASS as the stress-stiffening matrix slot for an appropriate buckling FULL file. This does not make an arbitrary static/modal FULL file a valid buckling pair.
- WRFULL is not supported in linear perturbation. Do not insert it blindly into Mechanical's buckling system.
- FULL matrices use solver ordering. MODE/RST results use internal ordering. Node numbering, internal ordering and solver equation ordering are different.
- Constraints, CP/CE equations and MPC connections can eliminate or transform DOFs. Do not infer equation rows as `(node-1)*6` or assume uniform six/seven-DOF layouts.
- The helper is not a DMP file combiner. Distributed files need a validated gathering route before this workflow.
- Associate exported mode-vector columns with the actual signed native multipliers. Total-deformation magnitudes are not signed eigenvectors.

## Phase 3: construct member contributions with the same transformations

For each selected element, the adapter must provide the actual elastic/geometric matrices and a transformation T such that `u_element = T q_solver_support`. This includes coordinate rotations, offsets, releases and the same constraint elimination used in global assembly.

Use `beamlba.assembly.ElementContribution` and `assemble_members` to form `T.T @ ke @ T` and `T.T @ kg @ T`, then scatter them into member supports. A principal submatrix of global K or G is **not** the target member's contribution: it can include connected elements.

Use `assembly_residuals` on a complete set of element/support/load-stiffness contributions to compare the reconstructed K and G with the exported global matrices. Document any contributions intentionally outside the member sets. An adapter that cannot reproduce the assembly must not be marked reviewed.

This release supplies the assembly functions, but the automatic ANSYS element-reader/transformation bridge remains to be implemented and validated for the target release. The inventory alone is insufficient to construct it.

## Phase 4: native parity and physical audit

Package data using `beamlba pack` or `write_bundle`; see BUNDLE_FORMAT.md. Record units, reference load case, proportional scaling, mapping convention, solver/version and the Mechanical snapshot hash.

Before setting review flags:

1. Compare several native eigenvalues and signed vectors against the extracted original pencil. The default gate requires at least three native pairs, or the system dimension for very small problems. Every supplied pair must meet the residual tolerance.
2. Check physical mapping by reconstructing selected native displacements/rotations/warping values and comparing them with ANSYS, including nonconsecutive node numbers and constrained joints.
3. Compare complete element assembly and verify member element ownership.
4. Review active beam formulation, warping continuity/releases, offsets, load application and the static/prestress link.
5. Compare selected targeted modes to a deeper native reference solve on a small frame before using the workflow on a larger project.

Only after those checks may `mapping_reviewed` and `element_formulation_reviewed` be true. Those booleans record an engineering review; they are not independently generated certificates. A low residual cannot detect a consistently wrong physical model or certify lowest-member-mode coverage.

## Phase 5: launch a reviewed bundle from Mechanical

Install the package in a separate CPython environment. In the launcher set PYTHON_EXE, BUNDLE_DIRECTORY, optional SEARCH_CONFIG and RUN_ASSESSMENT=True. Store the current inventory's `snapshot_sha256` as `metadata.mechanical_snapshot_sha256` in the prepared bundle. A mismatch aborts the run.

The launcher starts a numerical worker without a shell, redirects both output streams, restricts BLAS threads, and checks time and private-memory budgets. It can terminate that worker, not ANSYS. The private-memory threshold is a sampled watchdog rather than a preallocated hard memory reservation. A native allocation could fail before the next sample; no approval is issued after failure. The numerical CLI itself only has soft between-operation time limits.

The launcher is synchronous and can occupy Mechanical's scripting UI while the worker runs. No background-service or progress-dialog integration is included. Do not run it from an after-solve callback until reentrancy and project lifecycle have been validated.

## Primary references

- [ANSYS 2025 R2 HBMAT: buckling matrix meaning and constraints](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_cmd/Hlp_C_HBMAT.html)
- [ANSYS DOF ordering and MODE-to-solver conversion](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_apdl/apdldofordering.html)
- [ANSYS WRFULL limitations](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_cmd/Hlp_C_WRFULL.html)
- [ANSYS SMAT imports](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_cmd/Hlp_C_SMAT.html)
- [ANSYS EXPORT formats and restrictions](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_cmd/Hlp_C_EXPORT.html)
