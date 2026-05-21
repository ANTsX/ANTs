# Repository Map

Purpose: identify the main owner layers in this repository before opening broad source files.

Use this file when task ownership is unclear or a task crosses zebrafish workflow and upstream ANTs source boundaries.

## Owner layers

- Zebrafish operational wrapper: `ants_toRef.sh`
  Owns manifest discovery, Slurm submission, role-to-path defaults, per-row registration execution, job scripts, logs, and output naming for zebrafish brain registration.
- ANTs command-line examples: `Examples/`
  Owns many installed command implementations, including `antsRegistration`, `antsApplyTransforms`, image utilities, and test wrappers.
- Registration internals: `ImageRegistration/` and `Examples/itkantsRegistrationHelper.*`
  Own core metrics, transforms, optimizers, and registration plumbing.
- Segmentation internals: `ImageSegmentation/`
  Own segmentation filters and list-sample behavior.
- Shared utilities: `Utilities/`
  Own image IO, transform helpers, warping filters, command-line parsing helpers, and common ITK utilities.
- Tensor handling: `Tensor/`
  Own tensor-specific filters and reorientation behavior.
- Shell workflow scripts: `Scripts/`
  Own shipped ANTs orchestration wrappers; do not treat them as algorithmic truth.
- Build and dependency machinery: `CMake/`, `SuperBuild/`, `CMakeLists.txt`, `SuperBuild.cmake`
  Own configuration, external dependencies, compiler behavior, and build layout.
- Tests and fixtures: `Examples/TestSuite/`, `TestData/`
  Own regression command definitions and hashed fixture references.

## Local project boundary

This checkout is an ANTs source tree with local zebrafish registration workflow files at the root. Future agents should route full ANTs source tasks normally, but treat `ants_toRef.sh` and its `.agents` references as the operational zebrafish layer.
