# Zebrafish Data Contracts

Purpose: record the manifest, role, path, and output contracts owned by `ants_toRef.sh`.

Use this file before changing path resolution, role names, output naming, or manifest semantics.

## Environment contracts

- `SCRATCH` is required.
- `NAS` is optional and helps form the default `MANIFEST_DIR`.
- `ANTSPATH` may override the ANTs binary directory.
- `REF_AVG_2P` may override the default average 2P reference.
- `CPUS`, `MEM`, `MAIL_TYPE`, `MAIL_USER`, `JOB_CHDIR`, `MANIFEST_DIR`, and `MANIFEST_CSV` may override operational defaults.

## Manifest contract

Rows are parsed as:

```text
moving_fish,moving_role,fixed_fish,fixed_role,mov_override,fix_override
```

- Header row is skipped when `moving_fish` is `moving_fish_id`.
- Blank rows and `#` comment rows are skipped.
- UTF-8 BOM, Windows CRLF, and tabs are normalized.
- `mov_override` and `fix_override` take precedence over role-based path resolution.
- Empty `fixed_fish` defaults to `moving_fish` unless `fixed_role` is `avg_2p`.

## Supported roles

- `anatomy_2p`: preferred path is `subjects/<fish>/raw/2pA/<fish>_anatomy_2P_GCaMP.nrrd`, then lowercase `2p`, then `subjects/<fish>/fixed/anatomy_2P_ref_GCaMP.nrrd`.
- `avg_2p`: global average reference from `REF_AVG_2P`.
- `confocal_r1` or `confocal_round1`: files under `subjects/<fish>/raw/Rbest`, preferring `*GCaMP*.nrrd`.
- `confocal_rN` or `confocal_roundN`: files under `subjects/<fish>/raw/Rn`, preferring round/channel-specific `*GCaMP*.nrrd` patterns.

## Output contract

- Special case `anatomy_2p -> avg_2p` writes `reg_to_avg2p/rowN/2P_to_avg2p_` outputs and copies `_aligned.nrrd` to `anatomy_2P_in_avg2p.nrrd` after real runs.
- General role pairs write under `reg/<moving_group>_to_<fixed_group>/<moving_label>_to_<fixed_label>_rowN/`.
- If fixed fish differs from moving fish and is not global, fixed fish is included in the fixed label.
- Each run writes logs under `logs/`.
- Resolved manifests include row index, roles, resolved paths, output prefix, and status.

## Change rules

- Do not rename roles or canonical outputs without updating the router, stage map, and recent-changes log.
- Add new roles in all three places: path resolution, label mapping, and group mapping.
- Keep explicit path overrides higher precedence than inferred defaults.
