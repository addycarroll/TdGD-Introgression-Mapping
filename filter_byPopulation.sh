#!/usr/bin/env bash
#SBATCH --job-name=perpop_parent_filters
#SBATCH --output=perpop_parent_filters_%A_%a.out
#SBATCH --error=perpop_parent_filters_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G
#SBATCH --array=1-21
#SBATCH --requeue

set -euo pipefail

# ---------- CONFIG ----------
IN_DIR="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices"
POP_MANIFEST="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices/filter_byPopulation/parents_byPopulation.txt"
OUT_BASE="${IN_DIR}/filter_byPopulation/per_population_filters"

# 21 wheat chromosomes in array order
CHRS=(1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D)
CHR="${CHRS[$((SLURM_ARRAY_TASK_ID-1))]}"
CHR_LETTER="${CHR: -1}"                      # "A", "B", or "D"
IN_FILE="${IN_DIR}/Parents_chr${CHR}_D250_DP10_60_qual_miss_filt_alleleGenos.tsv"
mkdir -p "${OUT_BASE}"

# If input for this chromosome is missing, exit gracefully
if [[ ! -s "$IN_FILE" ]]; then
  echo "WARN: Missing input for chr${CHR}: $IN_FILE" >&2
  exit 0
fi

# ---------- helper: subset two columns (RP & DONOR) + first 3 meta cols ----------
subset_two_cols() {
  # $1 = input file, $2 = rp name, $3 = donor name, prints to stdout
  awk -F'\t' -v RP="$2" -v DN="$3" '
    BEGIN{OFS="\t"}
    NR==1{
      rpC=0; dnC=0
      for(i=1;i<=NF;i++){
        if($i==RP) rpC=i
        if($i==DN) dnC=i
      }
      if(!rpC || !dnC) exit 2
      print $1,$2,$3,$(rpC),$(dnC)  # header
      next
    }
    { print $1,$2,$3,$(rpC),$(dnC) }
  ' "$1"
}

# ---------- read manifest, robust to column order ----------
# Required header names (tab-delimited): Population   RP   DONOR_AB   DONOR_D
awk -F'\t' '
  BEGIN{OFS="\t"}
  NR==1{
    for(i=1;i<=NF;i++){
      h[$i]=i
    }
    need="Population\tRP\tDONOR_AB\tDONOR_D"
    split(need, req, "\t")
    for(j in req){
      if(!(req[j] in h)){
        msg = "ERROR: Required column \"" req[j] "\" not found in manifest header"
        print msg | "cat 1>&2"
        exit 1
      }
    }
    next
  }
  NF>=4 {
    # skip comment/blank lines
    if($0 ~ /^[[:space:]]*#/ || $0 ~ /^[[:space:]]*$/) next
    pop     = $(h["Population"])
    rp      = $(h["RP"])
    donorAB = $(h["DONOR_AB"])
    donorD  = $(h["DONOR_D"])
    print pop, rp, donorAB, donorD
  }
' "$POP_MANIFEST" | while IFS=$'\t' read -r POP RP DONOR_AB DONOR_D; do

  # Choose donor by subgenome
  if [[ "$CHR_LETTER" == "D" ]]; then
    DONOR="$DONOR_D"
  else
    DONOR="$DONOR_AB"
  fi

  POP_DIR="${OUT_BASE}/${POP}"
  SUBSET_DIR="${POP_DIR}/subsets"
  FAIL_DIR="${POP_DIR}/failing"
  KEEP_DIR="${POP_DIR}/kept"
  SUM_DIR="${POP_DIR}/summaries"
  DROP_UNIQ="${POP_DIR}/unique_dropped_markers.tsv"
  mkdir -p "$SUBSET_DIR" "$FAIL_DIR" "$KEEP_DIR" "$SUM_DIR"

  SUBSET_FILE="${SUBSET_DIR}/chr${CHR}_${POP}_${RP}_vs_${DONOR}.tsv"
  FAIL_FILE="${FAIL_DIR}/chr${CHR}_${POP}_failing.tsv"
  KEEP_FILE="${KEEP_DIR}/chr${CHR}_${POP}_kept.tsv"
  SUM_FILE="${SUM_DIR}/chr${CHR}_${POP}_summary.tsv"

  # 1) Create parent-only subset
  if ! subset_two_cols "$IN_FILE" "$RP" "$DONOR" > "$SUBSET_FILE" 2>/dev/null; then
    {
      echo -e "chr\tpopulation\trp\tdonor\ttotal\tkept\tdrop_missing\tdrop_het\tdrop_mono\tdropped_total"
      echo -e "chr${CHR}\t${POP}\t${RP}\t${DONOR}\t0\t0\t0\t0\t0\t0"
    } > "$SUM_FILE"
    echo "WARN: Missing parent column(s) for ${POP} on chr${CHR}: RP=${RP}, DONOR=${DONOR}" >&2
    continue
  fi

  # 2) Filter markers (drop if missing/het in either parent or monomorphic between parents)
  {
    echo -e "Marker\tFailReasons" > "$FAIL_FILE"
    echo -e "Marker\tREF\tALT\t${RP}\t${DONOR}" > "$KEEP_FILE"
  }
  awk -F'\t' -v OFS='\t' -v CHR="chr'"$CHR"'" -v RP_NAME="$RP" -v DN_NAME="$DONOR" -v FAIL_FILE="$FAIL_FILE" -v KEEP_FILE="$KEEP_FILE" '
    function norm(gt,   n,a,a1,a2){
      gsub(/\r/,"",gt); gsub(/[ \t]/,"",gt); gsub(/\|/,"/",gt)
      if(gt=="" || gt=="." || toupper(gt)=="NA" || gt=="./.") return "MISS"
      if(index(gt,"/")==0 && length(gt)==2) gt=substr(gt,1,1) "/" substr(gt,2,1)
      n=split(gt,a,"/"); if(n!=2 || a[1]=="." || a[2]==".") return "MISS"
      a1=toupper(a[1]); a2=toupper(a[2])
      if(a1!=a2) return "HET"
      return a1
    }
    NR==1{ next }  # header already written
    {
      total++
      rp = norm($4)
      dn = norm($5)
      reasons=""
      if(rp=="MISS" || dn=="MISS"){
        reasons="MISSING_PARENT"
      } else if(rp=="HET" || dn=="HET"){
        reasons=(reasons=="" ? "HET_PARENT" : reasons";HET_PARENT")
      } else if(rp==dn){
        reasons=(reasons=="" ? "MONOMORPHIC" : reasons";MONOMORPHIC")
      }

      if(reasons!=""){
        print $1, reasons >> FAIL_FILE
        dropped++
        if(index(reasons,"MISSING_PARENT")) drop_miss++
        if(index(reasons,"HET_PARENT"))     drop_het++
        if(index(reasons,"MONOMORPHIC"))    drop_mono++
      } else {
        kept++
        print $1, $2, $3, $4, $5 >> KEEP_FILE
      }
    }
    END{
      printf "chr\tpopulation\trp\tdonor\ttotal\tkept\tdrop_missing\tdrop_het\tdrop_mono\tdropped_total\n" > "'"$SUM_FILE"'"
      printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
             CHR, "'"$POP"'", RP_NAME, DN_NAME, total+0, kept+0, drop_miss+0, drop_het+0, drop_mono+0, dropped+0 >> "'"$SUM_FILE"'"
    }
  ' "$SUBSET_FILE"

  # 3) Maintain per-pop unique dropped markers across chromosomes (safe with flock)
  LOCK_FILE="${DROP_UNIQ}.lock"
  exec 9> "$LOCK_FILE"
  flock 9
  {
    TMP="${DROP_UNIQ}.tmp.$$"
    [[ -s "$DROP_UNIQ" ]] && cat "$DROP_UNIQ" > "$TMP" || : > "$TMP"
    if [[ -s "$FAIL_FILE" ]]; then
      tail -n +2 "$FAIL_FILE" | cut -f1 >> "$TMP"
    fi
    sort -u "$TMP" > "$DROP_UNIQ"
    rm -f "$TMP"
  }
  flock -u 9
  exec 9>&-
  echo "chr${CHR} ${POP}: subset=${SUBSET_FILE}; kept=${KEEP_FILE}; failing=${FAIL_FILE}; summary=${SUM_FILE}"
done
