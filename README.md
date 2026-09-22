# BeamLBA

**Member-informed elastic buckling searches for ANSYS beam models.**

BeamLBA addresses a specific problem: the first few global LBA modes of a complex frame often provide little information about some individual members. Instead of requesting the entire preceding spectrum, this package uses inexpensive member predictors to target relevant regions of the original full-system buckling problem.

**Version 0.1.0a1 is a research alpha.** The numerical engine is implemented and tested. The Mechanical inventory/launcher and an experimental MAPDL export helper are included, but a universally automatic, validated ANSYS element/constraint mapping adapter is **not yet implemented**. Dropping the script into an arbitrary project does not yet produce a complete member validation. No licensed ANSYS end-to-end execution has been performed.

## What works now

| Capability | Status |
|---|---|
| Symmetric `K + lambda G = 0` with positive-definite `K` | Implemented; congruence scaling and input checks |
| Member predictors with surrounding elastic restraint | Implemented; cached LU and bounded block solves |
| Clustered shifted full-frame searches | Implemented; residual filtering and unresolved-member refinement |
| Reuse of supplied native eigenpairs | Implemented; algebraic parity is required |
| Member energy attribution | Implemented as a **candidate** indicator, not an LTB classifier |
| Transformed element assembly and assembly residual checks | Implemented; requires explicitly supplied transformations |
| Matrix Market import and portable bundles | Implemented; explicit signs/storage, SHA-256 checks, no pickle |
| HTML/JSON reports and mode-vector export | Implemented; no implicit design pass or arbitrary member Mcr |
| Mechanical Named Selection inventory | Implemented and mock-tested; read-only by default |
| Mechanical-launched numerical worker | Implemented; reviewed bundle required; runtime unverified |
| Automatic ANSYS release-specific element/DOF adapter | **Open integration gate** |
| Certified lowest-LTB coverage, nonlinear analysis, code resistance checks | **Not implemented** |

## Install and run the verification fixtures

Use a separate CPython environment, not Mechanical's embedded IronPython environment. Python 3.10 or newer is required for the numerical worker.

```powershell
git clone https://github.com/afemengineer/BeamLBA.git
cd BeamLBA
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\beamlba.exe demo --kind beam --out results/beam
.\.venv\Scripts\beamlba.exe demo --kind crowded --out results/crowded
```

On Linux, use `python3 -m venv .venv` and the executables in `.venv/bin/`. Output directories must be new; existing results are never silently overwritten.

Open `results/beam/assessment/report.html` or `results/crowded/assessment/report.html`. Machine-readable reports and signed equation-space vectors are stored alongside the HTML.

## Numerical verification performed

Local environment: CPython 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pytest 9.0.2. **68 tests passed.** The suite covers numerical invariance, convergence, invalid models, native parity, member ownership, import gates, reporting, and mocked Mechanical calls.

- A 16-element Vlasov beam gives `Mcr = 112887.4853 N m`, against the analytical `112887.3206 N m`: relative error `1.4588e-6` (about `0.000146%`).
- The algebraic crowded-spectrum fixture recovers the target at positive spectral rank **281** with **one shifted search requesting eight eigenpairs**. Its factor agrees with the independent dense reference, `18.8676629663`.

The crowded fixture is **not an ANSYS frame model**. The dense calculation is used only as the demonstration's independent reference, not by the targeting workflow. These results establish numerical checks, not industrial speed-up or structural-design validation. CI is configured for Linux/Windows and Python 3.10/3.12; local results do not imply those CI jobs have run.

## First use inside Mechanical

1. Work on a saved project copy. Keep the existing static analysis and LBA unchanged.
2. Define one whole-member Named Selection per physical member, such as `MEMBER_RAFTER_001`. A member may contain several line bodies. Overlapping element ownership is rejected; shared nodes are allowed.
3. Open `mechanical/beamlba_mechanical.py` in the Scripting panel. Set the exact static and LBA analysis names and run with its default `RUN_ASSESSMENT = False`.
4. Inspect the generated `inventory.json` under your home directory's `BeamLBA` folder. This step does **not** run a buckling check.

To run targeted numerical assessment, first establish a reviewed equation-space bundle using the procedure in [ANSYS integration](docs/ANSYS_INTEGRATION.md). Then configure the CPython executable, bundle directory and `RUN_ASSESSMENT = True`. The launcher never changes supports, section definitions, beam KEYOPTs or prestress links.

## Assessment of a prepared bundle

```powershell
beamlba verify path/to/bundle
beamlba run path/to/bundle --out results/case_01 --config examples/search_config.json
```

ANSYS-labelled bundles require native eigenpair parity plus explicit mapping and element-formulation review. `--diagnostic` permits an unreviewed bundle only with **DIAGNOSTIC_ONLY** results; it does not bypass a failed parity check.

`CANDIDATE_REQUIRES_REVIEW` means an original-system eigenpair with significant member energy has been recovered. It does not establish that the mode is LTB, that it is the lowest relevant mode, or that the member passes a resistance check. `UNRESOLVED` never means safe.

## Documentation and development

- [ANSYS integration and release gates](docs/ANSYS_INTEGRATION.md)
- [Numerical method and limits](docs/NUMERICAL_METHOD.md)
- [Bundle and mapping contract](docs/BUNDLE_FORMAT.md)
- [Validation and acceptance plan](docs/VALIDATION.md)
- [Implementation roadmap](docs/ROADMAP.md)
- [Contributor/agent instructions](AGENTS.md)

No customer models, solver results, or credentials are included. Keep actual project data out of this public repository. A distribution license has not been selected by the repository owner.
