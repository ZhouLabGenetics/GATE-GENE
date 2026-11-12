# GATE-GENE

Here we maintain GATE-GENE, gene- and region-based rare-variant association testing for
**time-to-event (survival)** phenotypes in large, related biobank cohorts. Please find
<https://zhoulabgenetics.github.io/GATE-GENE-doc/> for documentation.

Results from GATE-GENE are browsable at **<https://gategene.org>**, which hosts the
PheWAS and multi-state trajectory analyses in the UK Biobank and the *All of Us*
Research Program.

GATE-GENE is not a separate program: it is the time-to-event mode of the SAIGE R package
in this repository, selected with `--traitType=survival`. The two-step workflow, the
executables and the input formats are SAIGE's; what changes is that step 1 fits a
survival null model. For SAIGE and SAIGE-GENE+ themselves, see
<https://saigegit.github.io/SAIGE-doc/>.

The method

- tests a time to event with right censoring, rather than a binary case/control proxy
  that discards when the event happened
- accounts for sample relatedness through a frailty model with a genetic relationship
  matrix (GRM)
- allows for model fitting with either a full or a sparse GRM
- handles low event counts and heavy censoring through a saddlepoint approximation
- is computationally efficient for large data sets: the null model is fitted once per
  phenotype and reused for every gene
- performs single-variant association tests
- performs conditional association analysis (PLINK/PGEN input; see the documentation
  for the BGEN limitation)

Set-based tests for rare variants

- perform BURDEN, SKAT and SKAT-O tests
- allow for tests on multiple minor allele frequency cutoffs and functional annotations
- allow for specifying weights for markers in the set-based tests, defaulting to
  beta(MAF, 1, 25)
- combine results across annotation and MAF masks with the Cauchy combination

The package takes genotype file input in the following formats

PLINK (bed, bim, fam), PGEN, BGEN, VCF, BCF, SAV

## Installation

```bash
docker pull wzhou88/saigetest:1.5.2
docker run --rm wzhou88/saigetest:1.5.2 step1_fitNULLGLMM.R --help
```

Singularity, a source install with pixi, and pinning the image by digest for published
runs are covered under
[Installation](https://zhoulabgenetics.github.io/GATE-GENE-doc/docs/installation/).

## Getting started

A worked example that runs end to end on the bundled test data, with the expected
output committed for comparison, is in [`tutorial/`](tutorial/) and on the
[Example](https://zhoulabgenetics.github.io/GATE-GENE-doc/docs/example/) page.

```bash
docker run --rm -v "$PWD":/app -w /app \
  wzhou88/saigetest:1.5.2 bash tutorial/run_example.sh
```

## Citation

If you use GATE-GENE, please cite the GATE-GENE preprint.

## About

Developed in the [Zhou Lab](https://github.com/ZhouLabGenetics).

Please report problems through
[GitHub issues](https://github.com/ZhouLabGenetics/GATE-GENE/issues), including the
version, how you installed it, the full command and the complete log.
