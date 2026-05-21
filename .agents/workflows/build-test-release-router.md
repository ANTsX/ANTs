# Build Test Release Router

Purpose: route work on build configuration, dependency superbuilds, CI, Docker, CTest, ExternalData, and release mechanics.

Use this file when the task mentions CMake, SuperBuild, compilation, installation, CI, test data fetches, Docker, packaging, or release metadata.

## Read order

1. `coding.md`
2. `.agents/workflows/ants-router.md`
3. This router
4. `.agents/references/repo-map.md`
5. Owning build, CI, or test configuration

## Task routing table

| Task language | Read first | Owner layer |
| --- | --- | --- |
| top-level configuration, compiler flags, project options | `repo-map.md` | `CMakeLists.txt`, `Common.cmake`, `Version.cmake`, `CMake/` |
| dependency build, ITK/VTK/SlicerExecutionModel, ExternalProject issues | `repo-map.md` | `SuperBuild/`, `SuperBuild.cmake` |
| CTest, ExternalData, regression baselines | `symbol-index.md` | `Examples/TestSuite/`, `CTest*.cmake`, `TestData/` |
| GitHub Actions, AppVeyor, Docker | `repo-map.md` | `.github/workflows/`, `appveyor.yml`, `Dockerfile` |
| local build or install documentation | `repo-map.md` | `README.md`, build scripts, `Utilities/SetupForDevelopment.sh` |

## Ownership guidance

- Keep build-system changes separate from zebrafish operational docs unless the task explicitly connects them.
- Do not ignore or delete checked-in `TestData/*.sha512` or `.md5` fixtures.
- Prefer the narrowest CMake configure, build, or CTest command that proves the touched area.
- Log meaningful completed changes in `.agents/references/recent-changes-build-test-release.md`; track unresolved work in `.agents/references/remaining-work-build-test-release.md`.
