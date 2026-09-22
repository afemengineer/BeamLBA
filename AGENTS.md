# BeamLBA engineering rules

- Preserve the original full-system pencil for accepted eigenpairs. Elastic-restraint predictors are search aids, not original-frame results.
- Never report unresolved as safe, a small residual as complete coverage, or an elastic multiplier as design resistance.
- Never convert a full-frame combined-loading factor into per-member pure-bending Mcr without a separately justified analysis.
- Treat matrices, native vectors and mappings as one versioned analysis-stage dataset. Do not guess equation indices from node IDs or assume every node has six/seven DOFs.
- Use actual element contributions with the full coordinate/offset/release/constraint transforms. Global principal submatrices are not member contributions.
- Keep customer models, solver output, credentials and proprietary section/code tables out of this public repository.
- Never change Mechanical KEYOPTs, restraints, analysis links or solve state without an explicit user action and a project-copy workflow.
- Run `python -m pytest`, the beam demo and crowded demo. Add a regression for each numerical/integration fix. Mock tests do not establish licensed ANSYS runtime support.
- Keep declared scope, status tables and validation records aligned with implemented behaviour. Do not replace missing integrations with success-shaped placeholders.
- Read docs/ANSYS_INTEGRATION.md before editing the adapter. Read docs/NUMERICAL_METHOD.md before changing the eigensolver.
- Do not select a redistribution license on behalf of the owner without an explicit decision.
