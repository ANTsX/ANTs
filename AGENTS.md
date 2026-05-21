# Agent Entrypoint

This repository uses a routed instruction system under `.agents/`. Start with the generic coding baseline, then use the workflow routers to decide which small reference files and owning code areas to read before making changes.

## Required startup order

1. Open `coding.md` for baseline coding behavior.
2. Open `.agents/workflows/ants-router.md`.
3. Follow its dispatch table to the right workflow router.
4. Read the smallest relevant reference doc.
5. Open symbol or stage maps only if needed.
6. Open owning modules or scripts before broad source-tree searches.
7. Open large scripts or generated job bodies only when owner context is insufficient.

## Non-negotiable repo rules

- Treat `ants_toRef.sh` as the canonical zebrafish registration entrypoint.
- Treat the ANTs C++/CMake tree as upstream project source unless the task is explicitly about zebrafish operations.
- Preserve canonical manifest columns, role names, output paths, log names, and registration stage semantics unless the task is a deliberate migration.
- Do not use wrapper scripts as authority for ANTs algorithm internals; wrappers own orchestration, paths, resources, and operational defaults.
- Validate after edits with the smallest relevant check.

## Reference files

- `.agents/workflows/ants-router.md`: top-level workflow dispatcher.
- `.agents/workflows/zebrafish-registration-router.md`: zebrafish manifest, Slurm, path, output, and profile routing.
- `.agents/workflows/ants-core-router.md`: ANTs C++/CLI source routing.
- `.agents/workflows/build-test-release-router.md`: CMake, SuperBuild, CI, and release routing.
- `.agents/references/repo-map.md`: repository ownership map.
- `.agents/references/zebrafish-stage-map.md`: end-to-end zebrafish registration flow.
- `.agents/references/zebrafish-data-contracts.md`: manifest, role, path, and output contracts.
- `.agents/references/zebrafish-registration-profile.md`: default ANTs registration profile in `ants_toRef.sh`.
- `.agents/references/symbol-index.md`: compact public and owning surfaces.
- `.agents/references/refactor-rules.md`: repo-specific editing rules.
- `.agents/references/refactor-loop-policy.md`: working-slice and handoff policy.
- `.agents/references/current-state.md`: active local state and caveats.
- `.agents/references/recent-changes.md`: change-log router.
- `.agents/references/remaining-work.md`: remaining-work router.

## Scope note

Keep this file short. Details live under `.agents/`; `coding.md` is always loaded first as generic behavior, then workflow routers add repo-specific ownership and validation rules.
