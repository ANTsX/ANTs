# Zebrafish Stage Map

Purpose: describe the current end-to-end flow of `ants_toRef.sh`.

Use this file before changing execution order, Slurm behavior, job resources, provenance, logs, or outputs.

## Stages

1. Submit wrapper startup
   - Requires `SCRATCH`.
   - Defaults `ANTSPATH` to `$HOME/ANTs/antsInstallExample/install/bin`.
   - Defaults `REF_AVG_2P` to `/scratch/ddharmap/refBrains/ref_05_LB_Perrino_2p/average_2p.nrrd`.
   - Accepts only `--dry-run`/`-n`; unknown flags are ignored.

2. Partition and manifest discovery
   - Prompts for partition.
   - Test partition runs locally with 1 CPU, 8G, 30 minutes.
   - Manifest search order is explicit `MANIFEST_CSV`, `MANIFEST_DIR/regManifest.csv`, `$SCRATCH/regManifest.csv`, then `$SCRATCH/experiments/_manifests/register_pairs.csv`.
   - If no manifest exists, legacy mode asks for fish IDs.

3. Job script materialization
   - Writes a timestamped generated job under `$SCRATCH/experiments/_jobs`.
   - Copies the manifest beside the job for provenance and records a SHA-256 checksum when possible.
   - Injects reference path, manifest path, checksum, dry-run mode, and legacy fish list.

4. Row dispatch
   - Manifest mode submits or runs one job per parsed row using `MANIFEST_ROW`.
   - Test partition runs each row interactively from `$SCRATCH/experiments`.
   - Other partitions submit each row with `sbatch`, configured CPU, memory, time, mail, output log, and error log.

5. Row execution
   - Normalizes manifest lines for BOM, CRLF, and tabs.
   - Resolves moving and fixed paths from role defaults unless overrides are present.
   - Writes a resolved manifest under `$SCRATCH/experiments/_jobs`.
   - Skips rows with missing moving or fixed files.

6. Registration
   - Calls `register_pair fixed moving outprefix logfile`.
   - Dry-run writes the intended registration command to the row log.
   - Real runs execute the default Rigid -> Affine -> SyN profile.

7. Outputs
   - `anatomy_2p -> avg_2p` rows write under `subjects/<fish>/reg_to_avg2p/rowN/`.
   - Other role pairs write under `subjects/<moving_fish>/reg/<moving_group>_to_<fixed_group>/<run_tag>/`.
   - Logs live under each run directory's `logs/`.

## Validation notes

Default validation is `bash -n ants_toRef.sh`, `shellcheck ants_toRef.sh` when available, and dry-run/path-resolution checks. Do not submit Slurm jobs unless requested.
