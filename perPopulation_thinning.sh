#!/usr/bin/env bash
#SBATCH --job-name=perpop_thin_apply
#SBATCH --output=perpop_thin_apply_%A_%a.out
#SBATCH --error=perpop_thin_apply_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G
#SBATCH --array=1-21
#SBATCH --requeue

set -euo pipefail

# ========= CONFIG =========
IN_DIR="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices"
POP_MANIFEST="${IN_DIR}/filter_byPopulation/parents_byPopulation.txt"
BASE="${IN_DIR}/filter_byPopulation/per_population_filters"   # where kept/ live
OUT_ROOT="${IN_DIR}/filter_byPopulation/thinned_markers"     # where thinned/ will go

# Default bin size (bp) applied to ALL populations:
DEFAULT_BIN_BP=50000

# Populations that should use 250 kb bins (edit names to match your manifest exactly).
# Include whatever canonical label your manifest uses for “population 41” (e.g., TdGD41).
POPS_250KB=("")

# ==========================
mkdir -p "${OUT_ROOT}"

# Chromosome list (array index from SLURM_ARRAY_TASK_ID)
CHRS=(1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D)
CHR="${CHRS[$((SLURM_ARRAY_TASK_ID-1))]}"

# Helper: return 1 if $1 is in POPS_250KB
is_250kb_pop() {
  local p="$1"
  for x in "${POPS_250KB[@]}"; do
    [[ "$p" == "$x" ]] && return 0
  done
  return 1
}

# Helper: label for filenames/logs
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

# ---- Iterate populations from manifest (header-robust, only need Population) ----
awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{for(i=1;i<=NF;i++) h[$i]=i; if(!("Population" in h)){print "ERROR: Population col missing" | "cat 1>&2"; exit 1}; next}
  $0 !~ /^[[:space:]]*$/ && $0 !~ /^[[:space:]]*#/ {print $(h["Population"])}
' "$POP_MANIFEST" \
| sort -u \
| while IFS=$'\t' read -r POP; do
  KEEP_DIR="${BASE}/${POP}/kept"
  KEEP_FILE="${KEEP_DIR}/chr${CHR}_${POP}_kept.tsv"

  # Decide bin size for this population
  BIN_BP="$DEFAULT_BIN_BP"
  if is_250kb_pop "$POP"; then
    BIN_BP=250000
  fi
  BIN_LABEL="$(label_bp "$BIN_BP")"
  OUT_DIR="${OUT_ROOT}/${POP}"
  mkdir -p "$OUT_DIR"

  OUT_FILE="${OUT_DIR}/chr${CHR}_${POP}_thinned_${BIN_LABEL}.tsv"
  SUM_FILE="${OUT_DIR}/${POP}_thinned_${BIN_LABEL}_summary.tsv"

  # Write headers if new
  if [[ ! -s "$SUM_FILE" ]]; then
    echo -e "chr\tpopulation\tbin_bp\tkept_input\tkept_output" > "$SUM_FILE"
  fi

  if [[ ! -s "$KEEP_FILE" ]]; then
    echo "WARN: No kept file for ${POP} chr${CHR}: ${KEEP_FILE}" >&2
    echo -e "chr${CHR}\t${POP}\t${BIN_BP}\t0\t0" >> "$SUM_FILE"
    # still create an empty thinned file with header for consistency
    echo -e "Marker\tREF\tALT\tRP\tDONOR" > "$OUT_FILE"
    continue
  fi

# Thinning:
# 1) skip header, extract POS as last numeric block from Marker, prepend POS
# 2) sort by POS
# 3) greedy select: keep first, then keep if pos - last_kept >= BIN_BP
# 4) print original row (without POS) and restore header
TMP="$(mktemp)"
awk -F'\t' '
  NR==1{next}
  {
    m=$1
    n=match(m,/([0-9]+)[^0-9]*$/)
    if(n){
      pos=substr(m,RSTART,RLENGTH)+0
      if(pos>0){
        print pos "\t" $0
      }
    }
  }
' "$KEEP_FILE" | sort -k1,1n > "$TMP"

input_count=$(wc -l < "$TMP" | tr -d ' ' || echo 0)

echo -e "Marker\tREF\tALT\tRP\tDONOR" > "$OUT_FILE"

COUNT_FILE="${OUT_FILE}.count"
: > "$COUNT_FILE"

if (( input_count > 0 )); then
  awk -F'\t' -v OFS='\t' -v BIN="$BIN_BP" -v CF="$COUNT_FILE" '
    BEGIN{last=-1; kept=0}
    {
      pos=$1+0
      # row begins at field 2 (original kept line)
      if(last<0 || pos-last>=BIN){
        print $2,$3,$4,$5,$6
        last=pos
        kept++
      }
    }
    END{
      print kept > CF
      close(CF)
    }
  ' "$TMP" >> "$OUT_FILE"
else
  echo 0 > "$COUNT_FILE"
fi

rm -f "$TMP"

kept_output=$(cat "$COUNT_FILE" 2>/dev/null || echo 0)
rm -f "$COUNT_FILE"

echo -e "chr${CHR}\t${POP}\t${BIN_BP}\t${input_count}\t${kept_output}" >> "$SUM_FILE"

echo "chr${CHR} ${POP}: bin=${BIN_LABEL} input=${input_count} -> thinned=${kept_output} -> ${OUT_FILE}"
done
