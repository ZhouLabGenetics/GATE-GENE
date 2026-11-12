# GATE-GENE tutorial

An end-to-end time-to-event analysis on the bundled example data: sparse GRM, null
model, gene-based test — all four stages, with every command self-checked. No biobank access and no data application are needed. It is the
fastest way to confirm an installation works and to see what every input and output
looks like.

## Run it

With Docker, from the repository root:

```bash
docker run --rm -v "$PWD":/gate -w /gate \
  wzhou88/saigetest:1.5.2 bash tutorial/run_example.sh
```

With a pixi source install:

```bash
pixi run bash tutorial/run_example.sh
```

Output goes to `extdata/output_tutorial/`. On success the last lines read:

```
  columns (11): Region	Group	max_MAF	Pvalue	Pvalue_Burden	Pvalue_SKAT	BETA_Burden	SE_Burden	MAC	Number_rare	Number_ultra_rare
  rows (set x mask): 14
  P values: 14 non-NA, range 0.172 to 0.96

TUTORIAL PASSED
```

The image is `linux/amd64` only. On Apple Silicon it runs under emulation — Docker picks
the platform itself and warns that it does not match the host — and takes a few minutes;
native Linux is considerably faster.

## What it runs

| Step | Script | Produces |
|---|---|---|
| 1 | `createSparseGRM.R` | `sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx` and its `.sampleIDs.txt` |
| 2 | `step1_fitNULLGLMM.R` | `example_survival.rda`, `example_survival.varianceRatio.txt` |
| 3 | `step2_SPAtests.R` | `example_survival_groupTest_out.txt`, `.txt.singleAssoc.txt`, `.txt.index` |
| 4 | checks | column, row and P-value validation |

`step1_fitNULLGLMM.R` fits with the **full** GRM estimated from `--plinkFile`, not with
the sparse GRM built a moment earlier. That is deliberate — one run then exercises
`createSparseGRM.R` *and* the no-sparse-GRM path. To use the sparse GRM instead, pass
`--sparseGRMFile`, `--sparseGRMSampleIDFile` and `--useSparseGRMtoFitNULL=TRUE`.

The phenotype file `extdata/input/pheno_1000samples_survival.txt` supplies 1,000
samples with:

| Column | Role |
|---|---|
| `IND_ID` | sample ID (`--sampleIDColinphenoFile`) |
| `casecontrol` | **event indicator**, 1 = event observed, 0 = censored (`--phenoCol`) |
| `AgeOfEventFinal` | **event / follow-up time** (`--eventTimeCol`) |
| `X` | a covariate |

## Three things this example exists to teach

**1. `--phenoCol` and `--eventTimeCol` are not interchangeable.** One is the 0/1
indicator, the other the time. Swapping them runs to completion and produces
meaningless results.

**2. `--eventTimeBinSize` is required in 1.5.2.** Omit it and step 1 fits the null model,
then aborts with `object 'minEventTime' not found`. The variable is assigned only inside
that option's branch (`R/SAIGE_fitGLMM_fast.R:1149`) but read unconditionally at lines
1575 and 1592. It is not a no-op — event times are floored into bins of the given width,
so choose it deliberately rather than copying `1` blindly.

**3. Variant IDs must match the genotype file exactly.** This example uses the BGEN,
whose IDs (`1:1:A:C`) match `group_new_chrposa1a2.txt`. The bundled `group_new_snpid.txt`
is the same genes named `rs1`, `rs2`, … , which is what the bed/bim/fam needs. Pairing
either group file with the other genotype file finds nothing — and
**step 2 still exits 0 and prints "Analysis done!"** — and the set-based output file is
created **completely empty, 0 bytes, without even a header row**. (The companion
`.singleAssoc.txt` does get a header, and nothing else.) The only sign is two lines
buried in the log — `101  markers in 'RegionFile' are not in 'GenoFile'.` and, three
lines later, `The chunk is empty`. Step 3's existence check is what catches this.

## Reading the results

One row per set x mask, plus a `Cauchy` row per set combining its masks.

| Column | Meaning |
|---|---|
| `Region` | set name from the group file |
| `Group` | annotation class tested |
| `max_MAF` | MAF cutoff for this mask |
| `Pvalue` | SKAT-O P value (the headline) |
| `Pvalue_Burden`, `Pvalue_SKAT` | component tests |
| `BETA_Burden`, `SE_Burden` | burden effect and its standard error |
| `MAC` | minor allele count in the mask |
| `Number_rare`, `Number_ultra_rare` | variants contributing; ultra-rare ones are collapsed |

Unlike a binary-trait run, the survival output carries **no** `MAC_case` / `MAC_control`
columns.

Only masks with variants appear: in this example the `max_MAF=0.0001` masks are absent
because no variant is that rare, and `GENE3` is absent because its only annotation
(`intergenic`) is not among those requested.

A reference copy of the expected result is in `expected/`. **The run is not
deterministic**, so do not expect an exact match: `createSparseGRM.R` samples
`--numRandomMarkerforSparseKin` markers at random and step 1 samples
`--numRandomMarkerforVarianceRatio` markers for the variance ratio, and neither script
exposes a seed. Two runs of this tutorial on the same machine and image gave P values
differing by up to 1.8e-07. That is why the step-4 check is structural rather than an
exact comparison.

## If it fails

| Symptom | Cause |
|---|---|
| `object 'minEventTime' not found` | `--eventTimeBinSize` not set; see above |
| `FAILED at: no group-test output` | variant IDs in the group file do not match the genotype file — the result file is written 0 bytes |
| `Rscript: command not found` | expected inside the image — only the four step scripts are on `PATH`; the checks use `awk` for this reason |
| a set missing from the output | its annotations are not among those in `--annotation_in_groupTest`, or none of its variants pass the MAF cutoffs |
