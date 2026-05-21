# Refactor Rules

Purpose: define repo-specific ownership and edit-scope rules.

Use this file before refactoring, extracting logic, or touching multiple owner layers.

## Rules

- Keep zebrafish orchestration in `ants_toRef.sh` unless a deliberate migration is planned.
- Keep ANTs algorithm changes in C++ owners, not in zebrafish or shipped shell wrappers.
- Keep build behavior in CMake/SuperBuild owners.
- Fix bugs at the narrowest owning layer.
- Preserve manifest columns, role semantics, output paths, and registration stage order unless intentionally migrating them.
- Preserve public ANTs CLI compatibility unless the task explicitly asks for interface change.
- Do not patch downstream consumers to compensate for upstream semantic bugs.
- Do not treat generated job scripts under `$SCRATCH` as source-of-truth files.
- Avoid broad formatting churn in upstream ANTs files.
- If public workflow behavior changes, update the relevant router/reference and recent-changes log in the same pass.

## Validation expectation

Every meaningful edit needs validation. For zebrafish wrapper edits, default to shell/static checks and dry-run behavior. For ANTs source or build edits, default to the narrowest build or CTest surface available.
