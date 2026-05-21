# Refactor Loop Policy

Purpose: keep future agents working through a coherent ownership slice instead of stopping after a superficial change.

Use this file for refactors or multi-step cleanup.

## Default working unit

One ownership slice is the smallest owner that can be validated independently: one wrapper behavior in `ants_toRef.sh`, one ANTs command implementation, one CMake/build area, or one semantic reference update.

## Required loop

1. Identify owner and validation surface.
2. Make the narrowest coherent change.
3. Remove only orphaned code introduced by the change.
4. Run the smallest relevant validation.
5. Update references/logs when public behavior or future routing changed.

## Valid stop conditions

- The slice is complete and validated.
- A missing data, cluster, dependency, or build context blocks further validation.
- The next step belongs to a different owner or would change public behavior beyond the request.

## Invalid stop conditions

- One helper was cleaned while adjacent caller behavior in the same owner remains inconsistent.
- Docs were updated without matching code when the task required behavior change.
- Static reasoning was used as the only proof for an executable change.

## Handoff

If stopping with unresolved work, update the matching `remaining-work-*` tracker with the status, next likely breakpoint, and validation still needed.
