# Zebrafish Registration Profile

Purpose: document the default `antsRegistration` profile used by `ants_toRef.sh`.

Use this file before changing registration quality, transform stages, metrics, convergence, shrink factors, smoothing, interpolation, or precision.

## Current profile owner

The profile is owned by `register_pair` in `ants_toRef.sh`.

Call shape:

```bash
register_pair <fixed> <moving> <outprefix> <logfile>
```

Fixed is the target space. Moving is warped into fixed.

## Global options

- Dimensionality: `-d 3`
- Precision: `--float 1`
- Verbosity: `--verbose 1`
- Output: `-o [<outprefix>,<outprefix>_aligned.nrrd]`
- Interpolation: `WelchWindowedSinc`
- Winsorization: `[0.05,0.95]`
- Histogram matching: enabled
- Initial moving transform: `[fixed,moving,1]`

## Stage sequence

1. Rigid
   - Transform: `Rigid[0.1]`
   - Metric: `MI[fixed,moving,1,32,Regular,0.25]`
   - Convergence: `[200x200x200x0,1e-8,10]`
   - Shrink factors: `12x8x4x2`
   - Smoothing sigmas: `4x3x2x1vox`

2. Affine
   - Transform: `Affine[0.1]`
   - Metric: `MI[fixed,moving,1,32,Regular,0.25]`
   - Convergence: `[200x200x200x0,1e-8,10]`
   - Shrink factors: `12x8x4x2`
   - Smoothing sigmas: `4x3x2x1vox`

3. SyN
   - Transform: `SyN[0.25,6,0.1]`
   - Metric: `CC[fixed,moving,1,4]`
   - Convergence: `[200x200x200x200x10,1e-7,10]`
   - Shrink factors: `12x8x4x2x1`
   - Smoothing sigmas: `4x3x2x1x0vox`

## Operational rules

- Preserve Rigid -> Affine -> SyN ordering unless the task is explicitly profile redesign.
- Keep fixed/moving order visible in logs and dry-run output.
- If profile behavior changes, update this file and log the change in `recent-changes-zebrafish-registration.md`.
- Profile changes should be validated at least with shell checks and dry-run output; real registration validation requires user-approved data/cluster context.
