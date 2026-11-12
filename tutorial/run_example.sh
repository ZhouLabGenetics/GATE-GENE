#!/usr/bin/env bash
# GATE-GENE end-to-end example on the bundled data. See tutorial/README.md.
#
#   Docker:  docker run --rm -v "$PWD":/gate -w /gate \
#              wzhou88/saigetest:1.5.2 bash tutorial/run_example.sh
#   Source:  pixi run bash tutorial/run_example.sh
#
# Deliberately linear: this script is documentation that runs.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT/extdata" || { echo "cannot cd to $ROOT/extdata"; exit 1; }  # so ./input/... resolves
OUT=./output_tutorial
mkdir -p "$OUT"

# The step scripts live on PATH in the Docker image and in extdata/ in a source
# checkout. Resolve either way so one script covers both routes.
run_step() {
  local s="$1"; shift
  if command -v "$s" >/dev/null 2>&1; then "$s" "$@"; else Rscript "./$s" "$@"; fi
}

fail() { echo; echo "FAILED at: $1"; exit 1; }

echo "=== 1/4  sparse GRM ==========================================="
run_step createSparseGRM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --nThreads=2 \
    --outputPrefix="$OUT/sparseGRM" \
    --numRandomMarkerforSparseKin=1000 \
    --relatednessCutoff=0.125 || fail "createSparseGRM.R"

GRM="$OUT/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
[[ -s "$GRM" ]] || fail "sparse GRM not written ($GRM)"

echo
echo "=== 2/4  step 1: null survival model =========================="
# --traitType=survival selects GATE-GENE. --phenoCol is the EVENT INDICATOR
# (1 = event, 0 = censored); --eventTimeCol is the TIME. They are not
# interchangeable: swapping them runs without error and means nothing.
#
# --eventTimeBinSize is REQUIRED in 1.5.2. Without it step 1 fits the null
# model and then aborts with "object 'minEventTime' not found", because
# minEventTime is assigned only inside the bin-size branch
# (R/SAIGE_fitGLMM_fast.R:1149) but read unconditionally at lines 1575/1592.
# Note this is not a no-op: event times are floored into bins of this width.
run_step step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --phenoFile=./input/pheno_1000samples_survival.txt \
    --phenoCol=casecontrol \
    --eventTimeCol=AgeOfEventFinal \
    --covarColList=X \
    --sampleIDColinphenoFile=IND_ID \
    --traitType=survival \
    --eventTimeBinSize=1 \
    --outputPrefix="$OUT/example_survival" \
    --nThreads=2 \
    --LOCO=FALSE \
    --minMAFforGRM=0.01 \
    --skipModelFitting=FALSE \
    --tol=0.01 \
    --isCovariateOffset=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE || fail "step1_fitNULLGLMM.R"

[[ -s "$OUT/example_survival.rda" ]] || fail "null model .rda not written"
[[ -s "$OUT/example_survival.varianceRatio.txt" ]] || fail "varianceRatio.txt not written"

echo
echo "=== 3/4  step 2: gene-based test ==============================="
# No survival-specific options here: step 2 reads the trait type from the
# fitted model, so this call is identical in form to SAIGE-GENE+.
#
# BGEN, not PLINK: the bundled group file names variants as 1:1:A:C, which
# matches the BGEN. genotype_100markers.bim uses rs1/rs2/... instead, so a
# bed-based run reports "101 markers in 'RegionFile' are not in 'GenoFile'"
# and writes an empty result. Variant IDs must match the genotype file exactly.
run_step step2_SPAtests.R \
    --bgenFile=./input/genotype_100markers.bgen \
    --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
    --sampleFile=./input/samplelist.txt \
    --SAIGEOutputFile="$OUT/example_survival_groupTest_out.txt" \
    --GMMATmodelFile="$OUT/example_survival.rda" \
    --varianceRatioFile="$OUT/example_survival.varianceRatio.txt" \
    --groupFile=./input/group_new_chrposa1a2.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --chrom=1 \
    --AlleleOrder=ref-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --LOCO=FALSE || fail "step2_SPAtests.R"

RES="$OUT/example_survival_groupTest_out.txt"
if [[ ! -s "$RES" ]]; then
  if [[ -e "$RES" ]]; then
    echo "  $RES exists but is 0 bytes -- step 2 matched no variants."
    echo "  Almost always the variant IDs in the group file do not match the genotype"
    echo "  file. Look for \"markers in 'RegionFile' are not in 'GenoFile'\" above."
  fi
  fail "no group-test output"
fi

echo
echo "=== 4/4  checks ==============================================="
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i
          if(!("Pvalue" in h)) { print "  no Pvalue column"; exit 1 }
          printf "  columns (%d): %s\n", NF, $0; next }
  { n++
    p=$(h["Pvalue"])
    if (p!="NA") { np++; if (p+0<0 || p+0>1) { printf "  P out of range: %s\n", p; bad=1 }
                   if (mn==""||p+0<mn) mn=p+0; if (mx==""||p+0>mx) mx=p+0 } }
  END { printf "  rows (set x mask): %d\n", n
        printf "  P values: %d non-NA, range %.3g to %.3g\n", np, mn, mx
        if (n==0)  { print "  header but no data rows"; exit 1 }
        if (np==0) { print "  no usable P values"; exit 1 }
        if (bad)   exit 1 }
' "$RES" || fail "result sanity check"

echo
echo "TUTORIAL PASSED"
echo "results: extdata/$(basename "$OUT")/$(basename "$RES")"
