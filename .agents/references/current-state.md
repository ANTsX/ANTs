# Current State

Purpose: record active local state and caveats that future agents should not misread.

Use this file before migrations, broad refactors, or when git status looks surprising.

## Notes

- This checkout is a full ANTs source tree.
- `ants_toRef.sh` is a local zebrafish registration wrapper at the repo root and is currently the canonical operational entrypoint.
- `setupAgents.md` is the local source document for creating this routed workflow.
- At the time this workflow was created, `ants_toRef.sh` and `setupAgents.md` were untracked local files.
- The workflow docs intentionally cover the full ANTs tree, but zebrafish registration behavior is routed separately from upstream ANTs internals.
- `coding.md` should be regenerated from the appendix block in `setupAgents.md` if that baseline guidance changes.

## Practical caveats

- Do not assume local data exists outside the repo; zebrafish inputs live under environment-dependent `$SCRATCH` and `$NAS` paths.
- Do not submit cluster jobs unless explicitly asked.
- Real registration quality validation requires sample data and cluster/runtime context from the user.
