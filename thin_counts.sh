#!/usr/bin/env bash
#SBATCH --job-name=perpop_thin_counts
#SBATCH --output=perpop_thin_counts_%A_%a.out
#SBATCH --error=perpop_thin_counts_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-21
#SBATCH --requeue

set -euo pipefail

# ---------- CONFIG ----------
IN_DIR="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices"
POP_MANIFEST="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices/filter_byPopulation/parents_byPopulation.txt"
BASE="${IN_DIR}/filter_byPopulation/per_population_filters"
OUT_ROOT="${IN_DIR}/filter_byPopulation/thinning_scenarios"

# Set any number of thresholds (in base pairs) in ONE place:
THRESHOLDS_BP=(50000 100000 150000 200000 250000)
mkdir -p "${OUT_ROOT}"

# Chromosomes (array index driven by SLURM_ARRAY_TASK_ID)
CHRS=(1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D)
CHR="${CHRS[$((SLURM_ARRAY_TASK_ID-1))]}"

# ---------- helpers ----------
# Label like 100kb / 1Mb for headers
label_bp() {
  local bp=$1
  if (( bp % 1000000 == 0 )); then
    printf "%dMb" $((bp/1000000))
  elif (( bp % 1000 == 0 )); then
    printf "%dkb" $((bp/1000))
  else
    printf "%dbp" "$bp"
  fi
}

# Build dynamic header suffixes once
build_header_suffixes() {
  local suf=()
  for t in "${THRESHOLDS_BP[@]}"; do
    suf+=("kept_at_$(label_bp "$t")")
  done
  printf "%s" "$(IFS=$'\t'; echo "${suf[*]}")"
}

HDR_SUFFIXES="$(build_header_suffixes)"

# Greedy count for an arbitrary set of thresholds (passed as -v THRSTR="bp1,bp2,...")
# Input: TSV with POS in column 2, sorted numerically.
greedy_count() {
  local tsv="$1"
  local thr_csv
  thr_csv="$(IFS=, ; echo "${THRESHOLDS_BP[*]}")"
  awk -F'\t' -v THRSTR="$thr_csv" '
    function count_for_thresh(thresh,   last, cnt, i){
      last = -1
      cnt  = 0
      for(i=1;i<=N;i++){
        if(pos[i] == "") continue
        if(last < 0 || pos[i] - last >= thresh){
          cnt++
          last = pos[i]
        }
      }
      return cnt
    }
    BEGIN{
      split(THRSTR, TH, /,/)
      M = length(TH)
    }
    {
      if(NR>0 && NF>=2 && $2 ~ /^[0-9]+$/){
        N++
        pos[N]=$2+0
      }
    }
    END{
      if(N==0){
        for(i=1;i<=M;i++){ printf (i<M? "0\t" : "0\n") }
        exit
      }
      for(i=1;i<=M;i++){
        c = count_for_thresh(TH[i]+0)
        printf (i<M ? c "\t" : c "\n")
      }
    }
  ' "$tsv"
}

# ---------- main: iterate manifest populations ----------
# Grab just Population column (header robust)
awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){ h[$i]=i }
    if(!("Population" in h)){
      msg="ERROR: Required column \"Population\" not found in manifest header"
      print msg | "cat 1>&2"
      exit 1
    }
    next
  }
  $0 !~ /^[[:space:]]*$/ && $0 !~ /^[[:space:]]*#/ {
    print $(h["Population"])
  }
' "$POP_MANIFEST" | sort -u | while IFS=$'\t' read -r POP; do
  POP_DIR="${BASE}/${POP}"
  KEEP_DIR="${POP_DIR}/kept"
  OUT_DIR="${OUT_ROOT}/${POP}"
  mkdir -p "$OUT_DIR"

  KEEP_FILE="${KEEP_DIR}/chr${CHR}_${POP}_kept.tsv"
  OUT_PERCHR="${OUT_DIR}/chr${CHR}_${POP}_thin_counts.tsv"
  OUT_ALL="${OUT_DIR}/${POP}_thin_counts_all_chr.tsv"

  # Write dynamic headers if files are new
  if [[ ! -s "$OUT_PERCHR" ]]; then
    echo -e "chr\tpopulation\ttotal_markers\t${HDR_SUFFIXES}" > "$OUT_PERCHR"
  fi
  if [[ ! -s "$OUT_ALL" ]]; then
    echo -e "chr\tpopulation\ttotal_markers\t${HDR_SUFFIXES}" > "$OUT_ALL"
  fi
  if [[ ! -s "$KEEP_FILE" ]]; then
    echo "WARN: No kept file for ${POP} chr${CHR}: ${KEEP_FILE}" >&2
    echo -e "chr${CHR}\t${POP}\t0\t$(printf '0\t%.0s' $(seq 1 ${#THRESHOLDS_BP[@]} ) | sed 's/\t$//')" >> "$OUT_PERCHR"
    echo -e "chr${CHR}\t${POP}\t0\t$(printf '0\t%.0s' $(seq 1 ${#THRESHOLDS_BP[@]} ) | sed 's/\t$//')" >> "$OUT_ALL"
    continue
  fi

  # Extract POS from Marker (col1): take last numeric block; sort by POS
  TMP_POS="$(mktemp)"
  awk -F'\t' '
    NR==1{ next }
    {
      m=$1
      pos=""
      n=match(m, /([0-9]+)[^0-9]*$/)
      if(n){ pos=substr(m, RSTART, RLENGTH) }
      if(pos ~ /^[0-9]+$/){ print m "\t" pos }
    }
  ' "$KEEP_FILE" | sort -k2,2n > "$TMP_POS"

  total=0
  if [[ -s "$TMP_POS" ]]; then
    total=$(wc -l < "$TMP_POS" | tr -d ' ' || echo 0)
  fi

  # Default zeroes for all thresholds
  mapfile -t zero_arr < <(printf '0\n%.0s' $(seq 1 ${#THRESHOLDS_BP[@]} ))
  counts=("${zero_arr[@]}")
  if (( total > 0 )); then
    out="$(greedy_count "$TMP_POS" 2>/dev/null || true)"
    # Split on tabs into counts[], keep array length flexible
    IFS=$'\t' read -r -a parsed <<< "${out:-}"
    # Fill counts from parsed (only if present)
    for i in "${!parsed[@]}"; do
      # guard numeric
      if [[ "${parsed[$i]}" =~ ^[0-9]+$ ]]; then
        counts[$i]="${parsed[$i]}"
      fi
    done
  fi
  rm -f "$TMP_POS"

  # Join counts with tabs
  counts_tsv="$(IFS=$'\t'; echo "${counts[*]}")"
  echo -e "chr${CHR}\t${POP}\t${total}\t${counts_tsv}" >> "$OUT_PERCHR"
  echo -e "chr${CHR}\t${POP}\t${total}\t${counts_tsv}" >> "$OUT_ALL"
  echo "chr${CHR} ${POP}: total=${total}; $(for i in "${!THRESHOLDS_BP[@]}"; do printf '%s=%s ' "$(label_bp "${THRESHOLDS_BP[$i]}")" "${counts[$i]:-0}"; done)"
done
