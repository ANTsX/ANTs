# ANTs Router

Purpose: dispatch future agent work to the smallest useful workflow docs for this full ANTs source tree plus the local zebrafish registration wrapper.

Use this file when starting any task in this repository.

## Read this first

1. `coding.md`
2. This router
3. The matching workflow router
4. The smallest relevant file under `.agents/references/`
5. Owning code, script, or CMake module

## Workflow profile dispatch

| Primary target or task language | Open first |
| --- | --- |
| `ants_toRef.sh`, zebrafish, fish IDs, 2P, confocal, manifests, roles, Slurm jobs, `$SCRATCH`, `$NAS`, average 2P reference, registration outputs | `.agents/workflows/zebrafish-registration-router.md` |
| ANTs command behavior, C++ registration internals, segmentation, image utilities, tensor handling, CLI examples, `antsRegistrationCLP` | `.agents/workflows/ants-core-router.md` |
| CMake, SuperBuild, Docker, CI, CTest, ExternalData, build failures, install paths, release or packaging mechanics | `.agents/workflows/build-test-release-router.md` |
| Cross-cutting refactor or unclear ownership | Read `.agents/references/repo-map.md`, then choose one profile before editing |

## Cross-workflow invariants

- Apply `coding.md` first, then layer workflow-specific rules.
- Preserve canonical outputs and stage semantics unless the user requests a migration.
- Prefer owner-level fixes over downstream compensation.
- Keep zebrafish operational behavior in `ants_toRef.sh` and its references; keep ANTs algorithm changes in the upstream source owners.
- Validate every meaningful edit with the smallest relevant check.
- Update recent-changes logs for behavior, public workflow, validation, output, or ownership changes.

## Compact scaling rule

Add tasks to an existing profile when they share entrypoints, validation, outputs, and owner boundaries. Create a new workflow profile only when a task family needs its own stage map, reference docs, validation cadence, and handoff tracker.
