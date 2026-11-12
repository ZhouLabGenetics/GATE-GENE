#!/usr/bin/env bash
# GATE-GENE: time-to-event (survival) command reference.
#
# This file is the fuller command list, including the sparse-GRM null model and
# the --is_fastTest variant. Every command here has been run against
# wzhou88/saigetest:1.5.2, but only ../tutorial/run_example.sh is self-checking --
# use that one to verify an installation.
#
# Run from this directory (extdata/). Inside the Docker image the four step
# scripts are on PATH, so drop the "Rscript" prefix.
#
# The phenotype file ./input/pheno_1000samples_survival.txt supplies:
#   IND_ID           sample ID
#   casecontrol      event indicator, 1 = event observed, 0 = censored
#   AgeOfEventFinal  event / follow-up time
#   X                a covariate
#
# NOTE: these commands use the BGEN, whose IDs (1:1:A:C) match
# group_new_chrposa1a2.txt. For the bed/bim/fam use ./input/group_new_snpid.txt
# instead, which names the same genes rs1/rs2/... . Crossing the two matches nothing
# and writes a 0-byte result file -- not even a header -- while exiting 0.
#
# NOTE: --eventTimeBinSize is required in 1.5.2. Without it step 1 aborts with
# "object 'minEventTime' not found" after fitting the null model. See
# R/SAIGE_fitGLMM_fast.R:1149 (assigned) vs 1575/1592 (read).
#
# Output goes to ./output_survival/, NOT ./output/. extdata/output/ holds example files
# that are committed to the repository -- including example_survival.rda and
# example_survival.varianceRatio.txt -- which these commands would otherwise overwrite.
set -euo pipefail
mkdir -p output_survival

##### 0. sparse GRM (once per cohort; reused by every phenotype) ##############

Rscript createSparseGRM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --nThreads=4 \
    --outputPrefix=./output_survival/sparseGRM \
    --numRandomMarkerforSparseKin=1000 \
    --relatednessCutoff=0.125

##### 1. step 1: fit the null survival model ##################################
# --traitType=survival is what selects GATE-GENE.
# --phenoCol is the EVENT INDICATOR; --eventTimeCol is the TIME. Swapping them
# runs without error and produces meaningless results.

Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --phenoFile=./input/pheno_1000samples_survival.txt \
    --phenoCol=casecontrol \
    --eventTimeCol=AgeOfEventFinal \
    --covarColList=X \
    --sampleIDColinphenoFile=IND_ID \
    --traitType=survival \
    --eventTimeBinSize=1 \
    --outputPrefix=./output_survival/example_survival \
    --nThreads=4 \
    --LOCO=FALSE \
    --minMAFforGRM=0.01 \
    --skipModelFitting=FALSE \
    --tol=0.01 \
    --isCovariateOffset=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE

# Same fit, using the sparse GRM from step 0 instead of the full GRM.
Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --phenoFile=./input/pheno_1000samples_survival.txt \
    --phenoCol=casecontrol \
    --eventTimeCol=AgeOfEventFinal \
    --covarColList=X \
    --sampleIDColinphenoFile=IND_ID \
    --traitType=survival \
    --eventTimeBinSize=1 \
    --outputPrefix=./output_survival/example_survival_sparseGRM \
    --nThreads=4 \
    --LOCO=FALSE \
    --minMAFforGRM=0.01 \
    --skipModelFitting=FALSE \
    --tol=0.01 \
    --isCovariateOffset=TRUE \
    --IsOverwriteVarianceRatioFile=TRUE \
    --sparseGRMFile=./output_survival/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=./output_survival/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE
# --useSparseGRMforVarRatio=TRUE is required if step 2 will be given the sparse GRM.
# It adds a "sparse" row to the .varianceRatio.txt file. Without it step 2 aborts with
#   "sparse GRM is specified but the variance ratio for sparse GRM was not estimated
#    in Step 1. Please remove --sparseGRMFile and --sparseGRMSampleIDFile"
# --useSparseGRMtoFitNULL=TRUE on its own only fits the null model; it writes no
# "sparse" row. Verified against wzhou88/saigetest:1.5.2.

##### 2. step 2: gene-based (set-based) test ##################################
# No survival-specific options here: step 2 reads the trait type from the
# fitted null model, so this call is identical in form to SAIGE-GENE+.
#
# annotation x maxMAF are crossed: 3 x 3 masks per set, plus the Cauchy
# combination across them.

Rscript step2_SPAtests.R \
    --bgenFile=./input/genotype_100markers.bgen \
    --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
    --sampleFile=./input/samplelist.txt \
    --SAIGEOutputFile=./output_survival/example_survival_groupTest_out.txt \
    --GMMATmodelFile=./output_survival/example_survival.rda \
    --varianceRatioFile=./output_survival/example_survival.varianceRatio.txt \
    --groupFile=./input/group_new_chrposa1a2.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --chrom=1 \
    --AlleleOrder=ref-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --LOCO=FALSE

# With the sparse-GRM null model, and the sparse GRM supplied to step 2 as well.
Rscript step2_SPAtests.R \
    --bgenFile=./input/genotype_100markers.bgen \
    --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
    --sampleFile=./input/samplelist.txt \
    --SAIGEOutputFile=./output_survival/example_survival_sparseGRM_groupTest_out.txt \
    --GMMATmodelFile=./output_survival/example_survival_sparseGRM.rda \
    --varianceRatioFile=./output_survival/example_survival_sparseGRM.varianceRatio.txt \
    --sparseGRMFile=./output_survival/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=./output_survival/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=./input/group_new_chrposa1a2.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --chrom=1 \
    --AlleleOrder=ref-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --LOCO=FALSE \
    --is_fastTest=TRUE

##### output ##################################################################
# ./output_survival/example_survival_groupTest_out.txt          one row per set x mask
#   Region Group max_MAF Pvalue Pvalue_Burden Pvalue_SKAT BETA_Burden SE_Burden
#   MAC Number_rare Number_ultra_rare
#   (verified on a survival run in wzhou88/saigetest:1.5.2. Unlike a binary-trait
#    run, the survival output carries no MAC_case / MAC_control columns.)
# One row per set x mask, plus a "Cauchy" row per set combining the masks.
# ./output_survival/example_survival_groupTest_out.txt.singleAssoc.txt   single-variant results
# ./output_survival/example_survival_groupTest_out.txt.index              resume/progress index
#
# BETA_Burden from a survival fit is NOT on the log-odds scale, so a binary-trait
# interpretation does not carry over. Confirm the exact scale before reporting.
