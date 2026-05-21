# Symbol Index

Purpose: compact index of public or owner-level surfaces future agents should know before editing.

Use this file when a task mentions a function, script helper, CLI, or build surface and the owner is not obvious.

## Zebrafish wrapper

- `ants_toRef.sh`: canonical zebrafish registration entrypoint.
- `register_pair`: row-level ANTs registration profile owner.
- `resolve_role_path`: role-to-file resolver for 2P, average reference, and confocal rounds.
- `role_round`: extracts confocal round from role strings.
- `role_label`: maps roles to run-tag labels.
- `role_group`: maps roles to output grouping folders.
- `sanitize_tag`: sanitizes output and job-name fragments.

## ANTs command surfaces

- `Examples/antsRegistration*.cxx`: installed `antsRegistration` variants.
- `Examples/itkantsRegistrationHelper.h` and `.hxx`: registration helper implementation.
- `Examples/antsApplyTransforms.cxx`: transform application CLI.
- `Examples/include/antsRegistration.h`: public command entrypoint declaration.
- `antsRegistrationCLP/`: command-line parser wrapper and XML surface.

## Build and test surfaces

- `CMakeLists.txt`: top-level build configuration.
- `SuperBuild.cmake` and `SuperBuild/`: external dependency build orchestration.
- `Examples/TestSuite/CMakeLists.txt`: regression tests for command-line tools.
- `TestData/`: hashed data references and fixtures.
- `.github/workflows/`, `appveyor.yml`, `Dockerfile`: CI and container surfaces.
