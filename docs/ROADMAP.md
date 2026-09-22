# Roadmap and acceptance criteria

## 1. Release-specific ANSYS extraction adapter: first deployment blocker

Deliver a small supported fixture set and automate the matching-stage FULL/MODE/EMAT extraction, element ownership and reduction transforms. Reproduce native eigenpairs and complete assembled K/G, then reconstruct physical displacements/rotations/warping values. Include beam offsets, releases, CP/CE/MPC and distributed-file handling only when individually verified. Preserve the original model. Do not substitute an guessed identity transform or global principal submatrix for missing element data.

Acceptance: repeatable read-only Mechanical entry point produces a reviewed bundle and matches native ANSYS for every advertised fixture. Today the inventory and numerical worker exist; this bridge is not complete.

## 2. Kinematic classification and coverage

Develop rigid-motion removal and member-local bending/twist/warping descriptors. Add basis-invariant participation for repeated eigenspaces. Distinguish flexural, coupled frame and LTB candidates. A low residual certifies an eigenpair, not that it is the lowest relevant mode.

Acceptance: labelled fixtures, ambiguous outputs preserved, lower relevant modes found against deep reference spectra. Add interval/inertia evidence or another explicit coverage argument before setting coverage_certified true. Do not tune thresholds only to the current synthetic examples.

## 3. Scalability

Replace all-member-support predictors with validated interface-only condensation and member interior reconstruction. Add reusable sparse factorizations across compatible cases, optional scalable solver backends and bounded/adaptive reduced-space enrichment. Benchmark fill-in, memory and time; soft budgets alone are insufficient.

Acceptance: measured performance on representative licensed ANSYS exports, with numerical errors and unresolved fractions reported alongside speed-up. No blanket claim of one unique LTB mode per physical member.

## 4. Member resistance and reporting

Only after numerical and physical coverage are established, add an explicitly selected design-code layer, code applicability checks, member parameters and load-combination processing. A full-frame combined-loading factor must not be relabelled as pure-bending Mcr. Resistance approval is separate from ideal elastic bifurcation.

No automated or background work is scheduled by this roadmap.
