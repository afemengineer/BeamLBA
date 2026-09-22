# Validation record and acceptance gates

## Executed locally

Date: 2026-09-22. Environment: CPython 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pytest 9.0.2. Single-threaded OpenBLAS.

**68 automated tests passed.** Tests exercise the actual package source. Mechanical tests use mock objects, not a licensed runtime. A separate Windows/.NET worker test has not been executed. CI is configured for Linux and Windows with Python 3.10/3.12; consult actual workflow runs rather than assuming this configuration has executed.

### Physical analytical benchmark

Prismatic doubly symmetric beam under uniform moment, fork-ended lateral/torsional restraint and free end warping. Independent two-field cubic-Hermite Vlasov FE; four-point integration. E=210 GPa, G=E/2.6, Iz=8e-6 m4, J=2e-7 m4, Iw=2e-7 m6, L=6 m, M0=20000 N m.

Analytical Mcr: 112887.32064731122 N m.
16-element FE Mcr: 112887.48532687814 N m.
Relative difference: 1.458795956699177e-6.
Convergence is checked at 2, 4, 8 and 16 elements. This validates the fixture and numerical kernel, not an ANSYS beam implementation.

### Crowded-spectrum targeting

600 equation coordinates, 300 weakly elastically coupled algebraic components. Target component 281; its dominant positive mode is rank 281 in the independently computed dense spectrum. Reference multiplier: 18.867662966267417. Targeted recovery: approximately 18.867662966267442 with one shifted search requesting eight eigenpairs. The first 20 modes have negligible target participation in this fixture.

The full dense spectrum is used only to verify the demonstration; `assess` does not compute it. This fixture is not a 3D frame, does not validate load-height effects, and provides no benchmark of industrial ANSYS runtime or RAM savings.

### Negative and integration checks

The suite covers invalid symmetry/definiteness, coordinate-scaling invariance, singular/indefinite G, invalid DOFs, shared-node versus overlapping-element ownership, transformed element assembly, native sign/mapping mismatch, repeated independent eigenvectors, exhausted budgets, sparse native-vector import, tampered files, explicit triangular storage, output preservation, HTML escaping and diagnostic-only import gates.

## Required ANSYS release validation

Start with customer-independent synthetic APDL/Mechanical fixtures and publish only data authorized for this public repository. For each supported configuration:

1. Compare native and exported eigenpairs in the same independent equation space.
2. Reconstruct the complete assembly including non-beam and load-stiffness contributions.
3. Compare reconstructed physical mode shapes to ANSYS, not just eigenvalue scalars.
4. Check BEAM188/189 warping, unequal DOF layouts, offsets, orientations, releases, nonconsecutive IDs and CP/CE/MPC constraints separately.
5. Compare targeted members with a sufficiently deep native spectrum on a manageable frame; report missed/mixed/repeated modes.
6. Demonstrate mesh convergence and classification on symmetric beams, channels, unequal-flange sections and coupled frames before broadening advertised support.
7. Benchmark peak RAM, factorization count and wall time on representative industrial models. Never infer a speed-up from sparse mode count alone.

Until these gates are met, 'ANSYS integration complete', 'each member validated', and 'production-ready' are incorrect descriptions of this release.
