#!/usr/bin/env bash
#SBATCH --job-name=ai2_pedonly_donorprops
#SBATCH --output=ai2_pedonly_donorprops_%A_%a.out
#SBATCH --error=ai2_pedonly_donorprops_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=1-147

set -euo pipefail

# --- AlphaImpute2 availability (match your environment) ---
module load Python/3.10.4-GCCcore-11.3.0
export PATH="$HOME/.local/bin:$PATH"   # expects ~/.local/bin/AlphaImpute2

# --- Discover all (pop × chr) inputs ---
ROOT="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices/Merged_byPopulation/filtered_parentHomPoly"
LIST="${ROOT}/ai2_inputs.list"
LOCK="${LIST}.lock"
# Build the list once, atomically, under a lock
{
  flock -x 200
  if [[ ! -s "$LIST" ]]; then
    tmp="$(mktemp "${LIST}.XXXX")"
    find "${ROOT}"/Pop*/AI2 -maxdepth 1 -type f -name 'chr*_Pop*.genotypes' \
      | sort > "$tmp"
    mv -f "$tmp" "$LIST"
  fi
} 200>"$LOCK"

# Wait briefly if another task just created it but the FS hasn’t flushed yet
for _ in {1..10}; do
  [[ -s "$LIST" ]] && break
  sleep 0.2
done

IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" || true)
[[ -n "${IN:-}" ]] || { echo "No input for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

# --- Parse pop & chr from the path ---
# IN example: /.../Pop80/AI2/chr1A_Pop80.genotypes
POP_DIR=$(basename "$(dirname "$IN")")             # AI2
POP_ROOT=$(basename "$(dirname "$(dirname "$IN")")")  # Pop80
POP_NUM=${POP_ROOT#Pop}                             # 80
BASE=$(basename "$IN" .genotypes)                   # chr1A_Pop80
CHR=${BASE%%_*}                                     # chr1A
CHR_LETTER=${CHR: -1}                               # A|B|D
SUBKEY=$([[ "$CHR_LETTER" == "D" ]] && echo "D" || echo "AB")
POP_AI2_DIR="$(dirname "$IN")"                      # .../Pop80/AI2
PED_FILE="${POP_AI2_DIR}/pop${POP_NUM}_ped.txt"
[[ -s "$PED_FILE" ]] || { echo "ERROR: missing pedigree $PED_FILE"; exit 2; }

# --- Where AI2 will write outputs for this pop × chr ---
RUN_DIR="${POP_AI2_DIR}/runs_pedOnly_l1/${CHR}"
OUT_PREFIX="${RUN_DIR}/${CHR}"                      # yields ${RUN_DIR}/${CHR}.genotypes
mkdir -p "$RUN_DIR"

# --- Parent mapping (exactly your table) ---
declare -A pop_map
pop_map["30,AB,REC"]="RP2"; pop_map["30,AB,DON"]="WE30"; pop_map["30,D,REC"]="RP2"; pop_map["30,D,DON"]="TA2"
pop_map["32,AB,REC"]="RP2"; pop_map["32,AB,DON"]="WE32"; pop_map["32,D,REC"]="RP2"; pop_map["32,D,DON"]="TA2"
pop_map["34,AB,REC"]="RP2"; pop_map["34,AB,DON"]="WE34"; pop_map["34,D,REC"]="RP2"; pop_map["34,D,DON"]="TA2"
pop_map["41,AB,REC"]="RP1"; pop_map["41,AB,DON"]="WE41"; pop_map["41,D,REC"]="RP1"; pop_map["41,D,DON"]="TA2"
pop_map["42,AB,REC"]="RP1"; pop_map["42,AB,DON"]="WE42"; pop_map["42,D,REC"]="RP1"; pop_map["42,D,DON"]="TA2"
pop_map["46,AB,REC"]="RP1"; pop_map["46,AB,DON"]="WE46"; pop_map["46,D,REC"]="RP1"; pop_map["46,D,DON"]="TA2"
pop_map["80,AB,REC"]="RP1"; pop_map["80,AB,DON"]="WE80"; pop_map["80,D,REC"]="RP1"; pop_map["80,D,DON"]="TA1"
REC="${pop_map["${POP_NUM},${SUBKEY},REC"]:-}"
DON="${pop_map["${POP_NUM},${SUBKEY},DON"]:-}"
if [[ -z "$REC" || -z "$DON" ]]; then
  echo "ERROR: Missing parent mapping for pop ${POP_NUM} on ${SUBKEY}"; exit 2
fi

# --- AlphaImpute2 knobs (ped-only) ---

AlphaImpute2 \
  -genotypes "$IN" \
  -pedigree  "$PED_FILE" \
  -out       "$OUT_PREFIX" \
  -ped_only \
  -error 0 \
  -length 2 \
  -final_peeling_threshold .1 \

# --- Donor allele proportions (computed from AI2 output for this pop×chr) ---
OUT_GENO="${RUN_DIR}/${CHR}.genotypes"   # AI2 output (no header; ID then calls)
[[ -s "$OUT_GENO" ]] || { echo "ERROR: expected AI2 output $OUT_GENO"; exit 3; }
DP_DIR="${POP_AI2_DIR}/donor_props_l1"
mkdir -p "$DP_DIR"
OUT_TSV="${DP_DIR}/donor_props_${CHR}_Pop${POP_NUM}.tsv"

# Normalize whitespace to tabs (in-place via temp)
tmpnorm="$(mktemp)"; trap 'rm -f "$tmpnorm"' EXIT
awk -v OFS='\t' 'BEGIN{FS="[ \t]+"} {$1=$1}1' "$OUT_GENO" | tr -d '\r' > "$tmpnorm" && mv "$tmpnorm" "$OUT_GENO"

# Header
echo -e "Individual\tPopulation\tChr\tSubgenome\tN_informative_markers\tN_used_nonmissing\tDonorProp" > "$OUT_TSV"

# Compute donor proportion only for WIL{POP_NUM}.*
awk -v CHR="$CHR" -v SUBG="$CHR_LETTER" -v POPNUM="$POP_NUM" -v RECID="$REC" -v DONID="$DON" '
  BEGIN { FS=OFS="\t"; haveREC=haveDON=0; }
  {
    id=$1;
    if (id==RECID) { haveREC=1; for (i=2;i<=NF;i++) rec[i]=$i+0; next; }
    if (id==DONID) { haveDON=1; for (i=2;i<=NF;i++) don[i]=$i+0; next; }
    store[NR]=$0; nf[NR]=NF; maxN=NR;
  }
  END {
    if (!haveREC || !haveDON) {
      printf("ERROR: Missing REC(%s) or DON(%s) on %s\n", RECID, DONID, CHR) > "/dev/stderr";
      exit 1;
    }
    # Informative markers: parents both homozygous (0/2) and different
    n_inf=0;
    for (i=2; i in rec; i++) {
      r=rec[i]; d=don[i];
      if ((r==0 || r==2) && (d==0 || d==2) && r!=d) {
        informative[i]=(d==2?1:0); # 1=>DON=2/REC=0; 0=>DON=0/REC=2
        n_inf++;
      }
    }
    for (k=1; k<=maxN; k++) {
      split(store[k], F, FS);
      id=F[1];
      if (id ~ "^WIL" POPNUM "\\.") {
        donor_sum=0; used=0;
        for (i in informative) {
          if (i<=nf[k]) {
            g = F[i]+0;
            if (g==0 || g==1 || g==2) {
              if (informative[i]==1) { donor_sum += g; }      # DON=2
              else                   { donor_sum += (2 - g); }# DON=0
              used++;
            }
          }
        }
        denom = 2*used;
        prop = (denom>0 ? donor_sum/denom : "NA");
        printf("%s\tWIL%d\t%s\t%s\t%d\t%d\t%s\n", id, POPNUM, CHR, SUBG, n_inf, used, prop);
      }
    }
  }
' "$OUT_GENO" >> "$OUT_TSV"

# --- Mark completion for this pop×chr and aggregate per-pop when all 21 finish ---
TOUCH_DIR="${DP_DIR}/.done_Pop${POP_NUM}"
mkdir -p "$TOUCH_DIR"
touch "${TOUCH_DIR}/.done_${CHR}"

# Attempt aggregation once all 21 chromosomes for this population are done.
# Use a simple lock directory to avoid duplicate aggregation.
LOCKDIR="${DP_DIR}/.agg_lock_Pop${POP_NUM}"
EXPECTED=21

# Count .done_*; if all present, try to acquire lock and aggregate
have_done=$(ls -1 "${TOUCH_DIR}"/.done_* 2>/dev/null | wc -l | awk "{print \$1}")
if [[ "$have_done" -ge "$EXPECTED" ]]; then
  if mkdir "$LOCKDIR" 2>/dev/null; then
    echo "[Aggregate] ${POP_ROOT}: all ${EXPECTED} chromosomes done; aggregating donor proportions."
    # Gather only this population's donor_props files
    shopt -s nullglob
    files=( "${DP_DIR}/donor_props_"*"_Pop${POP_NUM}.tsv" )
    if (( ${#files[@]} )); then
      # Per-Individual × Subgenome (weighted by N_used_nonmissing across chromosomes)
      awk -F'\t' 'BEGIN{
        OFS="\t"; first=1;
      }
      FNR==1 { if (first) { first=0; next } else { next } }  # skip headers
      {
        id=$1; pop=$2; subg=$4; used=$6; dp=$7;
        if (dp=="NA" || dp=="") next;
        key=id OFS pop OFS subg;
        w=used+0; val=dp+0.0;
        sum_w[key]+=w; sum_wx[key]+=w*val;
      }
      END{
        print "Individual","Population","Subgenome","DonorProp_weighted","N_used_total";
        for (key in sum_w) {
          split(key,a,OFS);
          id=a[1]; pop=a[2]; subg=a[3];
          wtot=sum_w[key];
          prop=(wtot>0 ? sum_wx[key]/wtot : "NA");
          print id,pop,subg,prop,wtot;
        }
      }' "${files[@]}" \
      | sort -t$'\t' -k2,2 -k3,3 -k1,1 \
      > "${DP_DIR}/per_individual_subgenome_Pop${POP_NUM}.tsv"

      # Population × Subgenome summary (mean & SD across individuals; unweighted across individuals)
      awk -F'\t' 'BEGIN{ OFS="\t"; }
      NR==1 { next }
      {
        pop=$2; subg=$3; dp=$4;
        if (dp=="NA" || dp=="") next;
        key=pop OFS subg;
        n[key] += 1;
        d = dp - mean[key];
        mean[key] += d / n[key];
        d2 = dp - mean[key];
        m2[key] += d * d2;
      }
      END{
        print "Population","Subgenome","N_individuals","MeanDonorProp","SDDonorProp";
        for (key in n) {
          split(key,a,OFS);
          pop=a[1]; subg=a[2];
          count=n[key]; mu=mean[key];
          sd=(count>1 ? sqrt(m2[key]/(count-1)) : 0.0);
          print pop,subg,count,mu,sd;
        }
      }' "${DP_DIR}/per_individual_subgenome_Pop${POP_NUM}.tsv" \
      | sort -t$'\t' -k1,1 -k2,2 \
      > "${DP_DIR}/population_subgenome_summary_Pop${POP_NUM}.tsv"
      echo "[Aggregate] Wrote:"
      echo " - ${DP_DIR}/per_individual_subgenome_Pop${POP_NUM}.tsv"
      echo " - ${DP_DIR}/population_subgenome_summary_Pop${POP_NUM}.tsv"
    else
      echo "[Aggregate] No donor_props files found for ${POP_ROOT}; nothing to aggregate."
    fi
    rmdir "$LOCKDIR" || true
  fi
fi
echo "[DONE] ${POP_ROOT} ${CHR}"
