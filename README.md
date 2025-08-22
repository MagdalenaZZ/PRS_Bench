# PRS_Bench
Package for constructing PRS in a few different ways

# PRS-Bench

Benchmark polygenic risk score (PRS) construction with multiple methods using **Nextflow** + **Docker/Conda** and **Python/R** helpers. We pretend GWAS has been run and start from summary statistics.

# PRS-Bench


Benchmark different PRS construction methods under a unified, reproducible pipeline.


## Methods Implemented
- **C+T (PLINK2)**: LD clumping + P-value thresholding; scored with PLINK2 `--score`.
- **PRS-CS (stub)**: Wrap PRS-CS; currently a pass-through stub; replace with real call.
- **LDpred2 (stub)**: R skeleton using bigsnpr/bigstatsr; currently pass-through.
- **Toy pure-Python C+T**: For tiny matrices when PLINK/VCF not available.


## Why stubs?
To keep this repo runnable without large LD refs or heavy installs, PRS-CS and LDpred2 are stubs you can swap for the real implementations once you provide references and packages.


## Data Contracts
- **sumstats.tsv** columns: `SNP, CHR, BP, A1, A2, BETA, SE, P`
- **phenotypes.tsv** columns: `sample_id, pheno[, covariates...]`
- **Genotypes**: PLINK `--bfile prefix` or `--vcf path`.


## Running
See Quickstart in the top-level README of this canvas.


## Evaluation
Outputs `scores.tsv` (sample_id, PRS, pheno, covariates) and metrics (`metrics.json`) plus a calibration plot.


## Bias & Fairness Hooks
- Provide ancestry-specific LD references; evaluate within and across ancestry groups.
- Export stratified metrics by group (TODO in `evaluate_prs.py`).


## Roadmap
- Replace stubs with real PRS-CS and LDpred2 runners.
- Add cross-validation of hyperparameters for C+T.
- Add covariate-adjusted logistic/linear models.
- Stratified performance reports and plots.
- Container images pinned with digests.

---

## Features

* Methods: **Clumping+Thresholding (C+T)** via PLINK2, **PRS-CS**, **LDpred2**, and a minimalist **pure-Python C+T** fallback for tiny toy data.
* End-to-end pipeline with **Nextflow** (local or cloud), containerized tools.
* Clear inputs/outputs, metrics (AUC, R², calibration), and multi-ancestry-aware evaluation splits.
* Reproducible with **Docker** or **Conda**.

---

## Repo layout

```
prs-bench/
├── README.md
├── nextflow.config
├── main.nf
├── conf/
│   ├── docker.config
│   ├── conda.config
│   └── profiles.config
├── envs/
│   ├── base-conda.yml
│   ├── r-ldpred2.yml
│   └── prscs.yml
├── docker/
│   ├── Dockerfile.plink
│   ├── Dockerfile.r-ldpred2
│   └── Dockerfile.prscs
├── bin/
│   ├── score_prs.py
│   ├── evaluate_prs.py
│   ├── toy_ct_prs.py
│   ├── wrap_prscs.py
│   └── ldpred2.R
├── resources/
│   ├── toy/
│   │   ├── genotypes.vcf.gz (toy; optional)
│   │   ├── phenotypes.tsv
│   │   └── sumstats.tsv
│   └── ldref/
│       ├── eur.ldblk.hdf5 (placeholder)
│       └── afr.ldblk.hdf5 (placeholder)
└── outputs/ (created at runtime)
```

> **Note**: Large LD reference files are placeholders; document where to download (e.g., 1000G-based LD blocks). Users can swap references by ancestry.

---

## Inputs

* **Summary statistics** (`sumstats.tsv`):

  * Required columns: `SNP, CHR, BP, A1, A2, BETA, SE, P`.
* **Genotypes**: PLINK bed/bim/fam **or** VCF; per-sample IDs should match phenotype file.
* **Phenotypes** (`phenotypes.tsv`): columns: `sample_id`, `pheno` (0/1 for case-control or continuous), optional covariates: `age, sex, PCs...`.

---

## Quickstart

### 1) With Docker + Nextflow

```bash
# Install nextflow (>=23)
NXF_VER=23.10.0
curl -s https://get.nextflow.io | bash

# Run with docker profile
./nextflow run main.nf -profile docker \
  --sumstats resources/toy/sumstats.tsv \
  --geno_prefix /path/to/genos # prefix of PLINK bed/bim/fam, or --vcf /path/to/genotypes.vcf.gz \
  --pheno resources/toy/phenotypes.tsv \
  --methods ct,prscs,ldpred2 \
  --outdir outputs/toy
```

### 2) With Conda (no Docker)

```bash
mamba env create -f envs/base-conda.yml
mamba env create -f envs/r-ldpred2.yml
mamba env create -f envs/prscs.yml
conda activate prs-bench
nextflow run main.nf -profile conda [...]
```

