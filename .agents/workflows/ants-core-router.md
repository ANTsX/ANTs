# ANTs Core Router

Purpose: route work on upstream ANTs command-line tools, C++ image registration and segmentation internals, utilities, tensor code, and examples.

Use this file when the task is about ANTs behavior itself rather than the local zebrafish cluster wrapper.

## Read order

1. `coding.md`
2. `.agents/workflows/ants-router.md`
3. This router
4. `.agents/references/repo-map.md`
5. `.agents/references/symbol-index.md` if entrypoints or owners are unclear
6. Owning C++/script files
7. Relevant tests in `Examples/TestSuite/`

## Task routing table

| Task language | Read first | Owner layer |
| --- | --- | --- |
| `antsRegistration`, transforms, metrics, convergence, optimizer behavior, fixed/moving semantics | `symbol-index.md` | `Examples/antsRegistration*.cxx`, `Examples/itkantsRegistrationHelper.*`, `ImageRegistration/` |
| `antsApplyTransforms`, transform order, image warping, interpolation | `symbol-index.md` | `Examples/antsApplyTransforms.cxx`, `Utilities/itkWarp*` |
| segmentation, Atropos, label fusion, priors | `repo-map.md` | `ImageSegmentation/`, relevant `Scripts/` wrappers |
| tensor or vector image behavior | `repo-map.md` | `Tensor/`, `Utilities/itkWarpTensor*` |
| command-line wrapper scripts shipped with ANTs | `refactor-rules.md` | `Scripts/` |
| public CLI interface or generated CLP behavior | `symbol-index.md` | `antsRegistrationCLP/`, `Examples/include/` |

## Ownership guidance

- Make algorithmic fixes in the C++ owner, not in a shell wrapper.
- Keep script changes scoped to orchestration, argument handling, and shipped workflow defaults.
- Preserve public CLI semantics unless the task explicitly asks for an interface change.
- Prefer existing CTest coverage in `Examples/TestSuite/` for validation.
- Log meaningful completed changes in `.agents/references/recent-changes-ants-core.md`; track unresolved work in `.agents/references/remaining-work-ants-core.md`.
