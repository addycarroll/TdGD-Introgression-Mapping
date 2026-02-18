#!/usr/bin/env bash
#SBATCH --job-name=prep_ai2_min
#SBATCH --output=prep_ai2_min_%A_%a.out
#SBATCH --error=prep_ai2_min_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=1-147

set -euo pipefail

# --- CONFIG ---
BASE_DIR="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices/Merged_byPopulation/filtered_parentHomPoly"

# Build a one-time list of all expected input files (per-pop × per-chr)
LIST="${BASE_DIR}/prep_ai2_input.list"
if [[ ! -s "$LIST" ]]; then
  find "$BASE_DIR" -type f -name 'chr*_Pop*_merged.parentHomPoly.tsv' | sort > "$LIST"
fi

IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" || true)
[[ -n "${IN:-}" ]] || { echo "No input for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

# Derive Pop ID and CHR
# Path like .../Pop80/chr1A_Pop80_merged.tsv
PP=$(basename "$(dirname "$IN")")          # Pop80
PP_NUM=${PP#Pop}                            # 80
BASE=$(basename "$IN" .tsv)                 # chr1A_Pop80_merged
CHR=${BASE%%_*}                             # chr1A
CHR_LETTER=${CHR: -1}                       # A|B|D

OUT_DIR="$(dirname "$IN")/AI2"
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

SORTED="$OUT_DIR/${BASE}_sorted.tsv"
NUMERIC="$OUT_DIR/${BASE}_numeric.tsv"
AI2="$OUT_DIR/${CHR}_Pop${PP_NUM}.genotypes"
SNP="$OUT_DIR/${CHR}_Pop${PP_NUM}.markers.snp"
BAD="$LOG_DIR/${BASE}_unparsable_positions.txt"
DUPREP="$LOG_DIR/${BASE}_duplicate_positions.tsv"
echo ">>> ${CHR} Pop${PP_NUM}"

# --- Parent mapping (from table) ---
# For A/B chromosomes use Donor_AB; for D chromosomes use Donor_D.
rp=""
don_ab=""
don_d=""

case "$PP_NUM" in
  30) rp="RP2"; don_ab="WE30"; don_d="TA2" ;;
  32) rp="RP2"; don_ab="WE32"; don_d="TA2" ;;
  34) rp="RP2"; don_ab="WE34"; don_d="TA2" ;;
  41) rp="RP1"; don_ab="WE41"; don_d="TA2" ;;
  42) rp="RP1"; don_ab="WE42"; don_d="TA2" ;;
  46) rp="RP1"; don_ab="WE46"; don_d="TA2" ;;
  80) rp="RP1"; don_ab="WE80"; don_d="TA1" ;;
  *) echo "ERROR: unknown population $PP_NUM"; exit 2 ;;
esac

if [[ "$CHR_LETTER" == "A" || "$CHR_LETTER" == "B" ]]; then
  donor="$don_ab"
else
  donor="$don_d"
fi

# --- 1) Header renaming: KSTdGD{PP}.{id} -> WIL{PP}.{id}; RP -> rp; DONOR -> donor ---
# Keep Marker, REF, ALT asis.
TMP_RENAMED="$(mktemp)"
awk -F'\t' -v OFS='\t' -v PP="$PP_NUM" -v RP_NEW="$rp" -v DON_NEW="$donor" '
  function rename_sample(x,   a){
    if (x=="RP") return RP_NEW;
    if (x=="DONOR") return DON_NEW;
    if (match(x, "^KSTdGD" PP "\\.([0-9]+)$", a)) return "WIL" PP "." a[1];
    return x;
  }
  NR==1{
    for(i=1;i<=NF;i++) $i = rename_sample($i);
    print; next
  }
  { print }
' "$IN" > "$TMP_RENAMED"

# --- 2) Sort by numeric position parsed from Marker ---
{
  head -n1 "$TMP_RENAMED"
  tail -n +2 "$TMP_RENAMED" \
  | awk -F'\t' -v OFS='\t' -v BADF="$BAD" '
      {
        m=$1
        if (match(m, /[:_]+([0-9]+)(\.[0-9]+)?$/, a)) {
          pos=a[1]
          key=sprintf("%012d", pos) OFS m
        } else {
          key=sprintf("%012d", 999999999999) OFS m
          print m >> BADF
        }
        print key, $0
      }' \
  | sort -t $'\t' -k1,1 -k2,2 \
  | cut -f3-
} > "$SORTED"

# Report duplicate numeric positions
awk -F'\t' '
  NR==1{next}
  {
    m=$1
    if (match(m, /[:_]+([0-9]+)(\.[0-9]+)?$/, a)) pos=a[1]; else next
    c[pos]++
  }
  END{ for(p in c) if(c[p]>1) printf "%s\t%d\n", p, c[p] }
' "$SORTED" | sort -k1,1n > "$DUPREP"

# --- 3) Recode to 0/1/2/9 (REF/ALT-based), dropping REF/ALT from numeric matrix ---
awk -F'\t' -v OFS='\t' '
  function UP(x){ return toupper(x) }
  function code(gt, ref, alt,   g, n, arr, a1, a2, R, T){
    g=gt; gsub(/[[:space:]]+/,"",g); gsub(/\|/,"/",g)
    if (g=="" || g=="9" || g=="." || g=="./." || UP(g)=="NA" || g=="N/N") return 9
    n=split(g, arr, "/"); if (n!=2) return 9
    a1=UP(arr[1]); a2=UP(arr[2]); R=UP(ref); T=UP(alt)
    if (a1==R && a2==R) return 0
    if (a1==T && a2==T) return 2
    if ((a1==R && a2==T) || (a1==T && a2==R)) return 1
    return 9
  }
  NR==1{
    for(i=1;i<=NF;i++) H[$i]=i
    m=H["Marker"]; r=H["REF"]; a=H["ALT"]
    if(!m || !r || !a){ print "ERROR: need Marker/REF/ALT columns" > "/dev/stderr"; exit 2 }
    printf "Marker"
    k=0
    for(i=1;i<=NF;i++){
      if(i!=m && i!=r && i!=a){ IDX[++k]=i; printf OFS $i }
    }
    printf "\n"; next
  }
  {
    ref=$r; alt=$a
    printf "%s", $m
    for(i=1;i<=k;i++) printf OFS "%d", code($(IDX[i]), ref, alt)
    printf "\n"
  }
' "$SORTED" > "$NUMERIC"

# --- 4) SNP map (.markers.snp) ---
{
  echo -e "SNP_ID\tChr\tPos\tA1\tA2"
  awk -F'\t' -v OFS='\t' -v CHR="$CHR" '
    NR==1{
      for(i=1;i<=NF;i++) H[$i]=i
      m=H["Marker"]; r=H["REF"]; a=H["ALT"]; next
    }
    {
      id=$m; ref=$r; alt=$a; pos=""
      if (match(id, /[:_]+([0-9]+)(\.[0-9]+)?$/, z)) pos=z[1]
      print id, CHR, pos, toupper(ref), toupper(alt)
    }
  ' "$SORTED"
} > "$SNP"

# --- 5) Transpose to AI2 (.genotypes), rows = individuals---
awk -f /dev/fd/3 "$NUMERIC" 3<<'AWKX' > "$AI2"
BEGIN{ FS=OFS="\t" }
NR==1{
  n=NF
  for(i=2;i<=n;i++){ id[i]=$i; line[i]=$i }
  next
}
{
  for(i=2;i<=n;i++){ line[i]=line[i] OFS $i }
}
END{
  for(i=2;i<=n;i++) print line[i]
}
AWKX

echo ">>> Done: $IN
    Sorted      : $SORTED
    Numeric     : $NUMERIC
    SNP map     : $SNP
    AI2 geno    : $AI2
    UnparsedPos : $BAD (if any)
    Dup pos rep : $DUPREP"

