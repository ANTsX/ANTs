# Zebrafish Registration Router

Purpose: route work for registering zebrafish brains between individuals or within an individual using the root `ants_toRef.sh` cluster wrapper.

Use this file when the task mentions `ants_toRef.sh`, zebrafish brains, 2P anatomy, confocal rounds, manifests, Slurm rows, registration output layout, or the default ANTs profile used for these jobs.

## Read order

1. `coding.md`
2. `.agents/workflows/ants-router.md`
3. This router
4. `.agents/references/zebrafish-data-contracts.md` for manifest/path/output questions
5. `.agents/references/zebrafish-registration-profile.md` for ANTs parameter questions
6. `.agents/references/zebrafish-stage-map.md` for execution flow or cluster behavior
7. `ants_toRef.sh`

## Task routing table

| Task language | Read first | Owner layer |
| --- | --- | --- |
| Add or change manifest columns, row parsing, comments, CRLF/BOM handling, row selection | `zebrafish-data-contracts.md` | `ants_toRef.sh` manifest parsing |
| Add a role such as a new imaging channel, confocal naming pattern, or reference type | `zebrafish-data-contracts.md` | `resolve_role_path`, `role_label`, `role_group` |
| Change output folders, filenames, provenance, resolved manifests, logs, or collision behavior | `zebrafish-data-contracts.md` and `zebrafish-stage-map.md` | output-prefix and logging sections in `ants_toRef.sh` |
| Tune registration quality, metrics, transforms, iterations, shrink factors, smoothing, interpolation, histogram matching, float precision | `zebrafish-registration-profile.md` | `register_pair` |
| Change cluster resources, partition behavior, job names, Slurm flags, mail, threading, dry-run behavior | `zebrafish-stage-map.md` | outer submit wrapper and generated job header |
| Debug missing file errors or wrong paths | `zebrafish-data-contracts.md` | role resolution and environment defaults |
| Convert wrapper into a cleaner CLI or package | `current-state.md` and `refactor-rules.md` | plan migration first; preserve `ants_toRef.sh` compatibility |

## Ownership guidance

- `ants_toRef.sh` owns operational orchestration, not ANTs algorithm internals.
- The generated job body inside `ants_toRef.sh` owns row-level execution and per-row outputs.
- The manifest owns pair intent: moving fish, moving role, fixed fish, fixed role, and optional path overrides.
- Path resolution owns only default conventions; explicit manifest overrides take precedence.
- Registration profile edits must keep fixed/moving semantics clear: fixed is the target space, moving is warped into fixed.
- Do not submit real Slurm jobs by default. Validate with `bash -n`, `shellcheck` when available, and dry-run-style checks unless the user asks for cluster execution.
- Log meaningful completed changes in `.agents/references/recent-changes-zebrafish-registration.md`; track unresolved follow-up in `.agents/references/remaining-work-zebrafish-registration.md`.
