#!/usr/bin/env bash
# sync_min.sh  — minimal pull → canonicalize → publish
#
# WHAT IT DOES
#   PULL:
#     NAS/<owner>/<fish>/02_reg/00_preprocessing → WORK mirror (2p_anatomy + rbest + rn)
#     WORK mirror → SCRATCH/raw + SCRATCH/fixed:
#       2p_anatomy → raw/2pA, rbest → raw/Rbest, rn → raw/Rn, plus nas symlink + .origin.json
#   CANONICALIZE:
#     Read SCRATCH/reg (+ optional SCRATCH/reg_to_avg2p)
#     Write canonical files into WORK/subjects/<fish>/02_reg/_canonical:
#       <fish>_<source>_in_<space>.nrrd, where <space> ∈ {2p, ref, round1 (Rbest), roundN (Rn)}
#   PUBLISH:
#     Use SCRATCH/<fish>/nas symlink (truth) to find NAS/<owner>/<fish>
#     Copy canonical files → NAS/02_reg stages:
#       01_rbest-2p, 02_rn-rbest, 03_rn-2p, 04_rbest-ref, 05_rn-ref, 06_total-ref, 08_2pa-ref
#
# USAGE
#   sync_min.sh [--pull|--push|--both] [--owner NAME | --nas-project-root PATH]
#               [--force] [--dry-run]
#               <fishID1> [fishID2 ...]
#
# Examples
#   Pull from Matilde and stage compute:
#     NAS="/nas/FAC/FBM/CIG/jlarsch/default/D2c/07 Data" \
#     WORK="/work/FAC/FBM/CIG/jlarsch/default/Danin" \
#     SCRATCH="/scratch/$USER" \
#     ./sync_min.sh --pull --owner Matilde L395_f10
#
#   Push results back (nas symlink already present from pull):
#     ./sync_min.sh --push L395_f10
#
# ENV
#   NAS (required for first pull if no symlink exists), WORK, SCRATCH
#
# Minimal dependencies: bash, rsync, readlink, sed

set -euo pipefail
IFS=$'\n\t'

# ---------- Base roots ----------
NAS_BASE="${NAS:-}"        # e.g. /nas/.../07 Data   (used with --owner on first pull)
WORK_BASE="${WORK:-$HOME/WORK}/experiments"
SCRATCH_BASE="${SCRATCH:-/scratch/$USER}/experiments"
# Stage names (NAS/WORK)
STAGE_RBEST_2P="01_rbest-2p"
STAGE_RN_RBEST="02_rn-rbest"
STAGE_RN_2P="03_rn-2p"

# ---------- CLI ----------
MODE="both"                # pull | push | both
FORCE=0
DRY=0
OWNER=""                   # e.g. Matilde
NAS_PROJECT_ROOT=""        # e.g. "/nas/.../07 Data/Matilde" (alternative to --owner)

usage() {
  cat <<USAGE
Usage: $0 [--pull|--push|--both] [--owner NAME | --nas-project-root PATH]
          [--force] [--dry-run]
          <fishID1> [fishID2 ...]

Options:
  --pull               Only NAS → WORK → SCRATCH
  --push               Only SCRATCH → WORK/_canonical → NAS
  --both               (default) do both
  --owner NAME         When pulling first time, use \$NAS/NAME/<fish> as NAS subject
  --nas-project-root P Alternative base path like ".../07 Data/Matilde" (overrides --owner)
  --force              Allow overwrites at destinations
  --dry-run            Print actions without writing
Notes:
  - 00_preprocessing expects rbest (best round) and rn (remaining rounds).
  - rbest populates raw/Rbest; rn populates raw/Rn; 2p_anatomy populates raw/2pA.
USAGE
  exit 1
}

args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pull) MODE="pull"; shift ;;
    --push) MODE="push"; shift ;;
    --both) MODE="both"; shift ;;
    --owner) shift; OWNER="${1:-}"; [[ -n "$OWNER" ]] || { echo "ERR: --owner needs a name" >&2; exit 2; }; shift ;;
    --nas-project-root) shift; NAS_PROJECT_ROOT="${1:-}"; [[ -n "$NAS_PROJECT_ROOT" ]] || { echo "ERR: --nas-project-root needs a path" >&2; exit 2; }; shift ;;
    --force) FORCE=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) usage ;;
    --*) echo "ERR: unknown option $1" >&2; usage ;;
    *) args+=( "$1" ); shift ;;
  esac
done
[[ ${#args[@]} -ge 1 ]] || usage
FISH_IDS=( "${args[@]}" )

# ---------- utils ----------
log(){ echo "[$(date -Iseconds)] $*"; }
ensure_dir() {
  local p
  for p in "$@"; do
    if (( DRY )); then echo "DRY: mkdir -p $p"; else mkdir -p "$p"; fi
  done
}
rsync_cp() {
  local add=( -a --no-owner --no-group --chmod=ugo=rwX )
  (( FORCE )) || add+=( --ignore-existing )
  if (( DRY )); then
    echo "DRY: rsync ${add[*]} $*"
    return 0
  fi
  rsync "${add[@]}" "$@"
}
count_files() {
  local dir="$1"
  find "$dir" -type f 2>/dev/null | wc -l | tr -d ' '
}
extract_row_tag() {
  local p="$1"
  if [[ "$p" =~ /row([0-9]+)(/|$) ]]; then
    echo "row${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}
extract_run_tag() {
  local p="$1"
  if [[ "$p" =~ /reg/[^/]+/([^/]+)/ ]]; then
    local tag="${BASH_REMATCH[1]}"
    if [[ "$tag" != "logs" ]]; then
      echo "$tag"
      return 0
    fi
  fi
  if [[ "$p" =~ /_canonical/([^/]+)/ ]]; then
    local tag="${BASH_REMATCH[1]}"
    if [[ "$tag" != "logs" ]]; then
      echo "$tag"
      return 0
    fi
  fi
  echo ""
  return 1
}
extract_subtag() {
  local p="$1" t=""
  t="$(extract_run_tag "$p")"
  [[ -n "$t" ]] && { echo "$t"; return 0; }
  t="$(extract_row_tag "$p")"
  echo "$t"
}
reg_stage_for_path() {
  local p="$1"
  if [[ "$p" == *"/Rbest_to_2pA/"* ]]; then
    echo "$STAGE_RBEST_2P"
    return 0
  fi
  if [[ "$p" == *"/Rn_to_Rbest/"* ]]; then
    echo "$STAGE_RN_RBEST"
    return 0
  fi
  if [[ "$p" == *"/Rn_to_2pA/"* ]]; then
    echo "$STAGE_RN_2P"
    return 0
  fi
  echo ""
  return 1
}
find_scratch_subj() {
  local fish="$1" top="$SCRATCH_BASE/subjects/$fish"
  [[ -d "$top" ]] && { echo "$top"; return 0; }
  mapfile -t cand < <(ls -d "$SCRATCH_BASE"/*/subjects/"$fish" 2>/dev/null || true)
  if (( ${#cand[@]} )); then
    local best="" best_m=0 m
    for c in "${cand[@]}"; do
      if   [[ -d "$c/reg" ]];          then m=$(stat -c %Y "$c/reg" 2>/dev/null || echo 0)
      elif [[ -d "$c/reg_to_avg2p" ]]; then m=$(stat -c %Y "$c/reg_to_avg2p" 2>/dev/null || echo 0)
      else m=0; fi
      (( m > best_m )) && { best_m=$m; best="$c"; }
    done
    [[ -n "$best" ]] && { echo "$best"; return 0; }
  fi
  echo "$top"
}
resolve_nas_subject() {
  local subj="$1" fish="$2"
  if [[ -L "$subj/nas" ]]; then
    local t; t="$(readlink -f "$subj/nas" || true)"
    [[ -n "$t" && -d "$t" ]] && { echo "$t"; return 0; }
  fi
  if [[ -f "$subj/.origin.json" ]]; then
    local t; t="$(sed -n -E 's/.*"nas_subject_root"[[:space:]]*:[[:space:]]*\"([^\"]+)\".*/\1/p' "$subj/.origin.json" | head -n1)"
    [[ -n "$t" && -d "$t" ]] && { echo "$t"; return 0; }
  fi
  # Explicit override wins
  if [[ -n "$NAS_PROJECT_ROOT" ]]; then
    echo "$NAS_PROJECT_ROOT/$fish"
    return 0
  fi
  # Owner-based defaults
  if [[ -n "$NAS_BASE" && -n "$OWNER" ]]; then
    # Special-case: Matilde subjects live under an extra 'Microscopy/' level
    if [[ "$OWNER" == "Matilde" ]]; then
      echo "$NAS_BASE/$OWNER/Microscopy/$fish"
      return 0
    fi
    # Generic case for all other owners
    echo "$NAS_BASE/$OWNER/$fish"
    return 0
  fi
  return 1
}

# ---------- PULL ----------
pull_one() {
  local fish="$1"
  log "=== PULL $fish ==="
  local SCR_SUBJ; SCR_SUBJ="$(find_scratch_subj "$fish")"
  local NAS_SUBJ
  if ! NAS_SUBJ="$(resolve_nas_subject "$SCR_SUBJ" "$fish")"; then
    echo "ERR: cannot resolve NAS subject for $fish (provide --owner or --nas-project-root, and set \$NAS)." >&2
    return 0
  fi
  local PRE="$NAS_SUBJ/02_reg/00_preprocessing"
  if [[ ! -d "$PRE" ]]; then
    log "  WARN: missing $PRE — nothing to pull."
    return 0
  fi
  local WORK_SUBJ="$WORK_BASE/subjects/$fish"
  local WORK_PRE="$WORK_SUBJ/02_reg/00_preprocessing"
  local RAW="$SCR_SUBJ/raw" FIXED="$SCR_SUBJ/fixed"

  # 1) Mirror preprocessing to WORK (rbest + rn)
  ensure_dir "$WORK_PRE"
  for sub in 2p_anatomy rbest rn; do
    if [[ -d "$PRE/$sub" ]]; then
      ensure_dir "$WORK_PRE/$sub"
      rsync_cp "$PRE/$sub/" "$WORK_PRE/$sub/"
    else
      log "  INFO: no $sub in preprocessing"
    fi
  done

  # Choose best-round source (rbest only)
  local BEST_SRC=""
  if [[ -d "$WORK_PRE/rbest" ]]; then
    BEST_SRC="$WORK_PRE/rbest"
  else
    log "  WARN: no rbest found under 00_preprocessing."
  fi

  # 2) Stage SCRATCH raw/fixed (respect names for downstream)
  ensure_dir "$RAW/2pA" "$RAW/Rbest" "$RAW/Rn" "$FIXED"
  [[ -d "$WORK_PRE/2p_anatomy" ]] && rsync_cp "$WORK_PRE/2p_anatomy/" "$RAW/2pA/" || true
  [[ -n "$BEST_SRC"           ]] && rsync_cp "$BEST_SRC/"             "$RAW/Rbest/" || true
  [[ -d "$WORK_PRE/rn"         ]] && rsync_cp "$WORK_PRE/rn/"         "$RAW/Rn/" || true

  # Fixed refs (first GCaMP found; skip if already present unless --force)
  shopt -s nullglob
  local gA=( "$RAW/2pA/"*GCaMP*.nrrd )
  local gB=( "$RAW/Rbest/"*GCaMP*.nrrd )
  shopt -u nullglob
  if (( FORCE )) || [[ ! -f "$FIXED/anatomy_2P_ref_GCaMP.nrrd" ]]; then
    [[ ${#gA[@]} -gt 0 ]] && rsync_cp "${gA[0]}" "$FIXED/anatomy_2P_ref_GCaMP.nrrd"
  fi
  if (( FORCE )) || [[ ! -f "$FIXED/round1_ref_GCaMP.nrrd" ]]; then
    [[ ${#gB[@]} -gt 0 ]] && rsync_cp "${gB[0]}" "$FIXED/round1_ref_GCaMP.nrrd"
  fi

  # 3) Traceability on SCRATCH
  if (( DRY )); then
    echo "DRY: ln -sfn \"$NAS_SUBJ\" \"$SCR_SUBJ/nas\""
    echo "DRY: write $SCR_SUBJ/.origin.json"
  else
    ln -sfn "$NAS_SUBJ" "$SCR_SUBJ/nas"
    cat > "$SCR_SUBJ/.origin.json" <<JSON
{
  "fish_id": "$fish",
  "nas_subject_root": "$NAS_SUBJ",
  "work_subject_root": "$WORK_SUBJ",
  "scratch_subject_root": "$SCR_SUBJ",
  "created_utc": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "created_by": "${USER:-unknown}",
  "script": "sync_min.sh"
}
JSON
  fi
  log "  ✔ pulled $fish"
}

# ---------- NEW: mirror SCRATCH → WORK/reg (safety/visibility before push) ----------
mirror_to_work_one() {
  local fish="$1"
  local SCR_SUBJ; SCR_SUBJ="$(find_scratch_subj "$fish")"
  local WORK_SUBJ="$WORK_BASE/subjects/$fish"
  ensure_dir "$WORK_SUBJ"
  if [[ -d "$SCR_SUBJ/reg" ]]; then
    ensure_dir "$WORK_SUBJ/reg"
    rsync_cp "$SCR_SUBJ/reg/" "$WORK_SUBJ/reg/"
  fi
  if [[ -d "$SCR_SUBJ/reg_to_avg2p" ]]; then
    ensure_dir "$WORK_SUBJ/reg_to_avg2p"
    rsync_cp "$SCR_SUBJ/reg_to_avg2p/" "$WORK_SUBJ/reg_to_avg2p/"
  fi
}

# ---------- CANONICALIZE (NRRDs + MATRICES; logs excluded) ----------
canonicalize_one() {
  local fish="$1"
  log "=== CANONICALIZE $fish ==="
  local SCR_SUBJ; SCR_SUBJ="$(find_scratch_subj "$fish")"
  local REG="$SCR_SUBJ/reg"
  local REG_AVG="$SCR_SUBJ/reg_to_avg2p"
  if [[ ! -d "$REG" && ! -d "$REG_AVG" ]]; then
    log "  INFO: no reg/ or reg_to_avg2p/ — skip."
    return 0
  fi
  local CANON="$WORK_BASE/subjects/$fish/02_reg/_canonical"
  ensure_dir "$CANON"

  # Collect only NRRDs + transform files (mat/nii.gz). Do NOT include logs here.
  shopt -s nullglob globstar
  local files=()
  [[ -d "$REG"     ]] && files+=( "$REG"/**/*.nrrd "$REG"/**/*GenericAffine.mat "$REG"/**/*Warp.nii.gz "$REG"/**/*InverseWarp.nii.gz )
  [[ -d "$REG_AVG" ]] && files+=( "$REG_AVG"/**/*.nrrd "$REG_AVG"/**/*GenericAffine.mat "$REG_AVG"/**/*Warp.nii.gz "$REG_AVG"/**/*InverseWarp.nii.gz )
  shopt -u globstar
  (( ${#files[@]} )) || { log "  INFO: nothing to canonicalize"; return 0; }

  for src in "${files[@]}"; do
    local base="$(basename "$src")"
    local out="$base"
    local subtag=""
    subtag="$(extract_subtag "$src")"

    # Keep these two distinct
    if [[ "$base" == "2P_to_avg2p__aligned.nrrd" ]]; then
      out="${fish}_2P_in_ref.nrrd"
    elif [[ "$base" == "anatomy_2P_in_avg2p.nrrd" ]]; then
      out="${fish}_anatomy_2P_in_ref.nrrd"
    fi

    # New reg outputs (Rbest/Rn layout)
    if [[ "$src" == *"/reg/"* ]]; then
      tag="$base"
      case "$base" in
        *__aligned.nrrd)        tag="${base%__aligned.nrrd}" ;;
        *_0GenericAffine.mat)   tag="${base%_0GenericAffine.mat}" ;;
        *_1Warp.nii.gz)         tag="${base%_1Warp.nii.gz}" ;;
        *_1InverseWarp.nii.gz)  tag="${base%_1InverseWarp.nii.gz}" ;;
      esac
      tag_lc="${tag,,}"
      if [[ "$tag_lc" =~ ^(.*)_row[0-9]+$ ]]; then
        tag_lc="${BASH_REMATCH[1]}"
      fi

      round=""
      target=""
      if [[ "$tag_lc" =~ ^rbest_to_2pa(_.*)?$ ]]; then
        round="1"; target="2p"
      elif [[ "$tag_lc" =~ ^r([0-9]+)_to_rbest(_.*)?$ ]]; then
        round="${BASH_REMATCH[1]}"; target="r1"
      elif [[ "$tag_lc" =~ ^r([0-9]+)_to_2pa(_.*)?$ ]]; then
        round="${BASH_REMATCH[1]}"; target="2p"
      fi

      if [[ -n "$round" && -n "$target" ]]; then
        case "$base" in
          *__aligned.nrrd)        out="${fish}_round${round}_GCaMP_in_${target}.nrrd" ;;
          *_0GenericAffine.mat)   out="${fish}_round${round}_GCaMP_to_${target}_0GenericAffine.mat" ;;
          *_1Warp.nii.gz)         out="${fish}_round${round}_GCaMP_to_${target}_1Warp.nii.gz" ;;
          *_1InverseWarp.nii.gz)  out="${fish}_round${round}_GCaMP_to_${target}_1InverseWarp.nii.gz" ;;
        esac
      fi
    fi

    # ---------- NRRD renaming (generic) ----------
    out="${out/_in_avg2p.nrrd/_in_ref.nrrd}"
    out="${out/_to_avg2p__aligned.nrrd/_in_ref.nrrd}"
    out="${out/_to_ref__aligned/_in_ref}"
    out="${out/_to_ref_aligned/_in_ref}"
    out="$(echo "$out" | sed -E 's/_to_ref(_[^.]*)?_aligned\.nrrd$/_in_ref.nrrd/')"
    out="$(echo "$out" | sed -E 's/_aligned_2P\.nrrd$/_in_2p.nrrd/')"
    out="$(echo "$out" | sed -E 's/_to_2p([^.]*)?_aligned\.nrrd$/_in_2p.nrrd/')"
    out="$(echo "$out" | sed -E 's/_to_r1([^.]*)?_aligned\.nrrd$/_in_r1.nrrd/')"

    # infer for generic per-channel aligned in reg/
    if [[ "$out" =~ ^${fish}_round1_.*_aligned\.nrrd$ ]]; then
        out="${out/_aligned.nrrd/_in_2p.nrrd}"
    elif [[ "$out" =~ ^${fish}_round2_.*_aligned\.nrrd$ ]]; then
        out="${out/_aligned.nrrd/_in_r1.nrrd}"
    fi

    # ---------- MATRICES (mat/nii.gz) ----------
    if [[ "$base" == *GenericAffine.mat || "$base" == *Warp.nii.gz || "$base" == *InverseWarp.nii.gz ]]; then
      if [[ "$src" == *"/reg_to_avg2p/"* ]]; then
        # global: to_avg2p → to_ref
        out="${out/to_avg2p_/to_ref_}"
        [[ "$out" == ${fish}_* ]] || out="${fish}_$out"
      else
        # relative reg/: round1 (Rbest) → in_2p, roundN → in_r1
        if   [[ "$out" =~ ^round1_GCaMP_to_ref ]]; then out="${out/round1_GCaMP_to_ref/round1_GCaMP_in_2p}"
        elif [[ "$out" =~ ^round2_GCaMP_to_ref ]]; then out="${out/round2_GCaMP_to_ref/round2_GCaMP_in_r1}"; fi
        [[ "$out" == ${fish}_* ]] || out="${fish}_$out"
      fi
    fi

    local dst="$CANON/$out"
    if [[ -n "$subtag" ]]; then
      ensure_dir "$CANON/$subtag"
      dst="$CANON/$subtag/$out"
    fi
    if (( FORCE )) || [[ ! -f "$dst" ]]; then
      rsync_cp "$src" "$dst"
      log "  + $(basename "$src") -> $(basename "$dst")"
    else
      log "  = exists: $(basename "$dst")"
    fi
  done
}

# ---------- STAGE MATRICES & LOGS into WORK stage folders ----------
stage_mats_logs_one() {
  local fish="$1"
  local SCR_SUBJ; SCR_SUBJ="$(find_scratch_subj "$fish")"
  [[ -n "$SCR_SUBJ" ]] || { log "  WARN: no SCRATCH subject; skip mats/logs staging"; return 0; }
  local REG="$SCR_SUBJ/reg"
  local REG_AVG="$SCR_SUBJ/reg_to_avg2p"
  local REG_WORK="$WORK_BASE/subjects/$fish/02_reg"

  # ---- from reg/ (Rbest/Rn layout) ----
  if [[ -d "$REG" ]]; then
    shopt -s nullglob globstar
    for m in "$REG"/**/*0GenericAffine.mat "$REG"/**/*1Warp.nii.gz "$REG"/**/*1InverseWarp.nii.gz; do
      [[ -f "$m" ]] || continue
      stg="$(reg_stage_for_path "$m")"
      [[ -n "$stg" ]] || continue
      run_tag="$(extract_run_tag "$m")"
      dest_dir="$REG_WORK/$stg/transMatrices"
      [[ -n "$run_tag" ]] && dest_dir="$dest_dir/$run_tag"
      ensure_dir "$dest_dir"
      rsync_cp "$m" "$dest_dir/${fish}_$(basename "$m")"
    done
    for lf in "$REG"/**/logs/*; do
      [[ -f "$lf" ]] || continue
      stg="$(reg_stage_for_path "$lf")"
      [[ -n "$stg" ]] || continue
      run_tag="$(extract_run_tag "$lf")"
      dest_dir="$REG_WORK/$stg/logs"
      [[ -n "$run_tag" ]] && dest_dir="$dest_dir/$run_tag"
      ensure_dir "$dest_dir"
      rsync_cp "$lf" "$dest_dir/$(basename "$lf")"
    done
    shopt -u globstar nullglob
  fi

  # ---- from reg_to_avg2p/ (global ref: 08 / 04 / 05) ----
  if [[ -d "$REG_AVG" ]]; then
    shopt -s nullglob globstar
    # 2P → ref
    for m in "$REG_AVG"/**/2P_to_avg2p_*; do
      [[ -f "$m" ]] || continue
      bn="$(basename "$m")"; bn="${bn/to_avg2p_/to_ref_}"
      row_tag="$(extract_row_tag "$m")"
      dest_dir="$REG_WORK/08_2pa-ref/transMatrices"
      [[ -n "$row_tag" ]] && dest_dir="$dest_dir/$row_tag"
      ensure_dir "$dest_dir"
      rsync_cp "$m" "$dest_dir/${fish}_$bn"
    done
    # roundN GCaMP → ref
    for m in "$REG_AVG"/**/round*_GCaMP_to_avg2p_*; do
      [[ -f "$m" ]] || continue
      b="$(basename "$m")"
      r="1"; [[ "$b" =~ ^round([0-9]+)_ ]] && r="${BASH_REMATCH[1]}"
      bn="${b/to_avg2p_/to_ref_}"
      stg="04_rbest-ref"; [[ "$r" != "1" ]] && stg="05_rn-ref"
      row_tag="$(extract_row_tag "$m")"
      dest_dir="$REG_WORK/$stg/transMatrices"
      [[ -n "$row_tag" ]] && dest_dir="$dest_dir/$row_tag"
      ensure_dir "$dest_dir"
      rsync_cp "$m" "$dest_dir/${fish}_$bn"
    done
    # logs (global)
    for lf in "$REG_AVG"/**/logs/*; do
      [[ -f "$lf" ]] || continue
      row_tag="$(extract_row_tag "$lf")"
      dest_dir="$REG_WORK/08_2pa-ref/logs"
      [[ -n "$row_tag" ]] && dest_dir="$dest_dir/$row_tag"
      ensure_dir "$dest_dir"
      rsync_cp "$lf" "$dest_dir/$(basename "$lf")"
    done
    shopt -u globstar nullglob
  fi
}

# ---------- PUBLISH ----------
publish_one() {
  local fish="$1"
  log "=== PUBLISH $fish ==="
  local SCR_SUBJ; SCR_SUBJ="$(find_scratch_subj "$fish")"
  local NAS_SUBJ
  if ! NAS_SUBJ="$(resolve_nas_subject "$SCR_SUBJ" "$fish")"; then
    echo "ERR: cannot resolve NAS target for $fish (missing nas symlink?)." >&2
    return 0
  fi
  local CANON="$WORK_BASE/subjects/$fish/02_reg/_canonical"
  [[ -d "$CANON" ]] || { log "  INFO: no _canonical — skip."; return 0; }

  local ROOT="$NAS_SUBJ/02_reg"
  ensure_dir "$ROOT"

  shopt -s nullglob globstar
  local nrrds=( "$CANON"/**/*.nrrd )
  if (( ${#nrrds[@]} == 0 )); then
    if (( DRY )); then
      log "  INFO: no canonical NRRDs in $CANON (dry-run doesn't create them)."
    else
      log "  INFO: no canonical NRRDs in $CANON."
    fi
  else
    local nrrd_copied=0 nrrd_skipped=0
    for f in "${nrrds[@]}"; do
    [[ -f "$f" ]] || continue
    local bn="$(basename "$f")" dest="" stage="" sub="aligned"

    # per-channel in ref → 06_total-ref/aligned/roundN
    if [[ "$bn" =~ ^${fish}_round([0-9]+)_channel[0-9]+_.*_in_ref\.nrrd$ ]]; then
      local R="${BASH_REMATCH[1]}"
      stage="06_total-ref"; ensure_dir "$ROOT/$stage/aligned/round${R}"
      dest="$ROOT/$stage/aligned/round${R}/$bn"
    # 2P anatomy (either name) → 08_2pa-ref/aligned
    elif [[ "$bn" == "${fish}_anatomy_2P_in_ref.nrrd" || "$bn" == "${fish}_2P_in_ref.nrrd" ]]; then
      stage="08_2pa-ref"; ensure_dir "$ROOT/$stage/aligned"
      dest="$ROOT/$stage/aligned/$bn"
    # rbest in other rounds (best → rn)
    elif [[ "$bn" =~ ^${fish}_round1_.*_in_r[0-9]+\.nrrd$ ]]; then
      stage="$STAGE_RN_RBEST"; ensure_dir "$ROOT/$stage/aligned"
      dest="$ROOT/$stage/aligned/$bn"
    # rbest in 2p
    elif [[ "$bn" =~ ^${fish}_round1_.*_in_2p\.nrrd$ ]]; then
      stage="$STAGE_RBEST_2P"; ensure_dir "$ROOT/$stage/aligned"
      dest="$ROOT/$stage/aligned/$bn"
    # remaining rounds in rbest
    elif [[ "$bn" =~ ^${fish}_round([2-9][0-9]*)_.*_in_r1\.nrrd$ ]]; then
      stage="$STAGE_RN_RBEST"; ensure_dir "$ROOT/$stage/aligned"
      dest="$ROOT/$stage/aligned/$bn"
    # remaining rounds in 2p (if present)
    elif [[ "$bn" =~ ^${fish}_round([2-9][0-9]*)_.*_in_2p\.nrrd$ ]]; then
      stage="$STAGE_RN_2P"; ensure_dir "$ROOT/$stage/aligned"
      dest="$ROOT/$stage/aligned/$bn"
    # rbest ref (stage root)
    elif [[ "$bn" =~ ^${fish}_round1_.*_in_ref\.nrrd$ ]]; then
      stage="04_rbest-ref"; ensure_dir "$ROOT/$stage"
      dest="$ROOT/$stage/$bn"
    # roundN (N>=2) ref (stage root)
    elif [[ "$bn" =~ ^${fish}_round([2-9][0-9]*)_.*_in_ref\.nrrd$ ]]; then
      stage="05_rn-ref"; ensure_dir "$ROOT/$stage"
      dest="$ROOT/$stage/$bn"
    else
      log "  WARN: no NAS mapping for $bn (skipped)."
      continue
    fi

    if (( FORCE )) || [[ ! -f "$dest" ]]; then
      rsync_cp "$f" "$dest"
      nrrd_copied=$((nrrd_copied + 1))
      log "  → $stage/${dest##*/}"
    else
      nrrd_skipped=$((nrrd_skipped + 1))
      log "  = exists: $stage/${dest##*/}"
    fi
    done
    log "  INFO: NRRDs: $nrrd_copied copied, $nrrd_skipped skipped."
  fi
  shopt -u globstar nullglob
  # Publish staged transMatrices & logs mirrored from WORK → NAS
  local WORK_02="$WORK_BASE/subjects/$fish/02_reg"
  shopt -s nullglob globstar
  for stgdir in "$WORK_02"/*; do
    [[ -d "$stgdir" ]] || continue
    stg="$(basename "$stgdir")"

    if [[ -d "$stgdir/transMatrices" ]]; then
      local mcount
      mcount="$(count_files "$stgdir/transMatrices")"
      ensure_dir "$ROOT/$stg/transMatrices"
      rsync_cp "$stgdir/transMatrices/" "$ROOT/$stg/transMatrices/"
      log "  → $stg/transMatrices (${mcount} files)"
    fi
    if [[ -d "$stgdir/logs" ]]; then
      local lcount
      lcount="$(count_files "$stgdir/logs")"
      ensure_dir "$ROOT/$stg/logs"
      rsync_cp "$stgdir/logs/" "$ROOT/$stg/logs/"
      log "  → $stg/logs (${lcount} files)"
    fi
  done
  shopt -u globstar nullglob
  log "  ✔ published $fish"
}

# ---------- MAIN ----------
echo "Roots:
  NAS     : ${NAS_BASE:-"(not used unless first pull + --owner)"} 
  WORK    : $WORK_BASE
  SCRATCH : $SCRATCH_BASE
Mode=$MODE  Force=$FORCE  Dry=$DRY  Owner=${OWNER:-"-"}  NAS_ROOT=${NAS_PROJECT_ROOT:-"-"}"

for fish in "${FISH_IDS[@]}"; do
  [[ "$MODE" == "pull" || "$MODE" == "both" ]] && pull_one "$fish"
  # Mirror reg trees into WORK before canonicalizing/publishing
  [[ "$MODE" == "push" || "$MODE" == "both" ]] && mirror_to_work_one "$fish"
  canonicalize_one "$fish"
  stage_mats_logs_one "$fish"
  [[ "$MODE" == "push" || "$MODE" == "both" ]] && publish_one "$fish"
done

log "All done."
