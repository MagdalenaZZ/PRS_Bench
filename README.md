# PRS_Bench
Package for constructing PRS in a few different ways

# PRS-Bench

Benchmark polygenic risk score (PRS) construction with multiple methods using **Nextflow** + **Docker/Conda** and **Python/R** helpers. We pretend GWAS has been run and start from summary statistics.

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

---

## `nextflow.config`

```groovy
params.outdir = params.outdir ?: 'outputs/run'
params.methods = (params.methods ?: 'ct').split(',')
params.sumstats = params.sumstats
params.geno_prefix = params.geno_prefix
params.vcf = params.vcf
params.pheno = params.pheno
params.ldref = params.ldref ?: 'resources/ldref'

process {
  executor = 'local'
  errorStrategy = 'retry'
  maxRetries = 2
  withName: CT_* { cpus = 4; memory = '8 GB' }
  withName: PRSCS_* { cpus = 8; memory = '16 GB' }
  withName: LDPRED2_* { cpus = 8; memory = '16 GB' }
}

profiles {
  docker { includeConfig 'conf/docker.config' }
  conda  { includeConfig 'conf/conda.config' }
}
```

---

## `conf/docker.config`

```groovy
process.container = 'ubuntu:22.04'
params.plink_image   = 'biocontainers/plink2:v2.00-ava7'
params.r_image       = 'ghcr.io/rocker-org/r-ubuntu:22.04'
params.py_image      = 'python:3.11-slim'
```

---

## `conf/conda.config`

```groovy
conda.enabled = true
process.conda = 'envs/base-conda.yml'
```

---

## `envs/base-conda.yml`

```yaml
name: prs-bench
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.11
  - numpy
  - pandas
  - scikit-learn
  - matplotlib
  - plink2
  - bcftools
  - tabix
  - r-base
  - r-data.table
  - r-optparse
```

---

## `envs/r-ldpred2.yml`

```yaml
name: r-ldpred2
channels: [conda-forge, bioconda, defaults]
dependencies:
  - r-base
  - r-bigsnpr
  - r-bigstatsr
  - r-data.table
  - r-optparse
```

---

## `envs/prscs.yml`

```yaml
name: prscs
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.10
  - numpy
  - scipy
  - pandas
  - pip
  - pip:
      - prscs # or local clone if needed
```

---

## `main.nf`

```groovy
nextflow.enable.dsl=2

Channel
  .fromPath(params.sumstats)
  .set { CH_SUMSTATS }

// Choose genotype source: PLINK prefix or VCF

def GENO_CHANNEL = params.vcf ? Channel.fromPath(params.vcf) : Channel.fromPath(params.geno_prefix)

def METHODS = params.methods

process CT_CLUMP {
  tag { file(params.sumstats).baseName }
  publishDir "${params.outdir}/ct", mode: 'copy'
  container params.plink_image
  input:
    path sumstats from CH_SUMSTATS
    val geno from GENO_CHANNEL
  output:
    path 'ct.snplist'
  script:
    def sum = "$sumstats"
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    # Extract SNP list by P-value threshold and LD clump (C+T)
    awk 'NR==1 || $8 < 5e-8 {print \$1}' ${sum} > p5e-8.snps
    plink2 ${genoArg} --clump ${sum} --clump-field P --clump-p1 5e-8 --clump-p2 1e-2 \
      --clump-r2 0.1 --clump-kb 250 --out clumped
    cut -f3 clumped.clumps | tail -n +2 > ct.snplist
    """
}

process CT_SCORE {
  publishDir "${params.outdir}/ct", mode: 'copy'
  container params.plink_image
  input:
    path sumstats from CH_SUMSTATS
    val geno from GENO_CHANNEL
    path snplist from CT_CLUMP.out
  output:
    path 'ct.profile'
  script:
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    awk 'NR==1 || FNR==NR {a[$1]=1; next} a[$1]' ${snplist} ${sumstats} > ct.sumstats
    # PLINK2 scoring expects: SNP A1 BETA
    awk 'NR>1{print $1"\t"$4"\t"$6}' ct.sumstats > weights.txt
    plink2 ${genoArg} --score weights.txt 1 2 3 header-read cols=scoresums \
       --out ct
    mv ct.sscore ct.profile
    """
}

process PRSCS_RUN {
  when: METHODS.contains('prscs')
  publishDir "${params.outdir}/prscs", mode: 'copy'
  container params.py_image
  input:
    path sumstats from CH_SUMSTATS
  output:
    path 'prscs.weights.tsv'
  script:
    """
    python3 /bin/wrap_prscs.py \
      --sumstats ${sumstats} \
      --ldref ${params.ldref} \
      --out prscs.weights.tsv
    """
}

process LDPRED2_RUN {
  when: METHODS.contains('ldpred2')
  publishDir "${params.outdir}/ldpred2", mode: 'copy'
  container params.r_image
  input:
    path sumstats from CH_SUMSTATS
  output:
    path 'ldpred2.weights.tsv'
  script:
    """
    Rscript /bin/ldpred2.R \
      --sumstats ${sumstats} \
      --ldref ${params.ldref} \
      --out ldpred2.weights.tsv
    """
}

process SCORE_ALL {
  publishDir "${params.outdir}/scores", mode: 'copy'
  container params.py_image
  input:
    path pheno from Channel.fromPath(params.pheno)
    tuple path(weights), val(method) from (
      Channel.fromPath("${params.outdir}/ct/weights.txt").map{[it,'ct']} \
      | (METHODS.contains('prscs') ? Channel.fromPath("${params.outdir}/prscs/prscs.weights.tsv").map{[it,'prscs']} : Channel.empty()) \
      | (METHODS.contains('ldpred2') ? Channel.fromPath("${params.outdir}/ldpred2/ldpred2.weights.tsv").map{[it,'ldpred2']} : Channel.empty())
    )
    val geno from GENO_CHANNEL
  output:
    path 'scores.tsv'
  script:
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    python3 /bin/score_prs.py \
      --weights ${weights} --method ${method} \
      ${params.vcf ? '--vcf ' + geno : '--bfile ' + geno} \
      --pheno ${pheno} --out scores.tsv
    """
}

process EVALUATE {
  publishDir "${params.outdir}/eval", mode: 'copy'
  container params.py_image
  input:
    path scores from SCORE_ALL.out
  output:
    path 'metrics.json'
    path 'calibration_plot.png'
  script:
    """
    python3 /bin/evaluate_prs.py --scores ${scores} \
      --out_metrics metrics.json --out_plot calibration_plot.png
    """
}
```

---

## `bin/score_prs.py`

```python
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, subprocess, tempfile, os, json

def read_weights(path):
    # Accept: three columns [SNP, A1, BETA] or more
    df = pd.read_csv(path, sep='\t', header=None)
    if df.shape[1] >= 3:
        df = df.iloc[:, [0,1,2]]
        df.columns = ['SNP','A1','BETA']
    else:
        raise ValueError('weights must have at least 3 columns: SNP, A1, BETA')
    return df

def score_with_plink(bfile=None, vcf=None, weights=None, out='scores.tmp'):
    wfile = 'weights.tmp.txt'
    weights[['SNP','A1','BETA']].to_csv(wfile, sep='\t', index=False, header=True)
    if bfile:
        cmd = ['plink2','--bfile', bfile,'--score', wfile,'1','2','3','header-read','cols=scoresums','--out','prs']
    else:
        cmd = ['plink2','--vcf', vcf,'--score', wfile,'1','2','3','header-read','cols=scoresums','--out','prs']
    subprocess.run(cmd, check=True)
    # plink2 outputs .sscore
    sscore = pd.read_csv('prs.sscore', sep='\t')
    sscore = sscore.rename(columns={'IID':'sample_id','SCORE1_SUM':'PRS'})
    sscore[['sample_id','PRS']].to_csv(out, sep='\t', index=False)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--weights', required=True)
    ap.add_argument('--method', required=True)
    ap.add_argument('--bfile')
    ap.add_argument('--vcf')
    ap.add_argument('--pheno', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    W = read_weights(args.weights)
    out_scores = score_with_plink(bfile=args.bfile, vcf=args.vcf, weights=W)
    S = pd.read_csv(out_scores, sep='\t')
    P = pd.read_csv(args.pheno, sep='\t')
    M = S.merge(P, on='sample_id', how='inner')
    M.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
```

---

## `bin/evaluate_prs.py`

```python
#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, json
from sklearn.metrics import roc_auc_score, brier_score_loss, r2_score
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--scores', required=True)
ap.add_argument('--out_metrics', required=True)
ap.add_argument('--out_plot', required=True)
args = ap.parse_args()

df = pd.read_csv(args.scores, sep='\t')
metrics = {}
if 'pheno' in df.columns:
    if set(np.unique(df['pheno'])) <= {0,1}:
        # binary
        y = df['pheno'].values
        x = (df['PRS'] - df['PRS'].mean())/df['PRS'].std()
        # naive logistic-like score using PRS only
        # report AUC and Brier
        metrics['auc'] = float(roc_auc_score(y, x))
        metrics['brier'] = float(brier_score_loss(y, (x-x.min())/(x.max()-x.min()+1e-9)))
        # calibration curve (quantile bins)
        df['q'] = pd.qcut(x, 20, duplicates='drop')
        cal = df.groupby('q').agg(obs=('pheno','mean'), prs=('PRS','mean')).reset_index()
        plt.figure()
        plt.scatter(cal['prs'], cal['obs'])
        plt.xlabel('Mean PRS in bin')
        plt.ylabel('Observed case rate')
        plt.title('Calibration (by PRS quantile)')
        plt.savefig(args.out_plot, dpi=160)
    else:
        # continuous trait
        y = df['pheno'].values
        x = df['PRS'].values
        metrics['r2'] = float(r2_score(y, x))
        plt.figure()
        plt.scatter(x, y, s=6, alpha=0.6)
        plt.xlabel('PRS')
        plt.ylabel('Phenotype')
        plt.title('PRS vs Trait')
        plt.savefig(args.out_plot, dpi=160)

with open(args.out_metrics, 'w') as f:
    json.dump(metrics, f, indent=2)
```

---

## `bin/toy_ct_prs.py` (pure-Python, small toy VCF or matrix)

```python
#!/usr/bin/env python3
"""
Minimal clumping+thresholding for tiny toy data using a genotype matrix (CSV) and summary stats.
This is NOT production-grade; it ignores LD except for a window/r2 heuristic and assumes phased additive coding.
"""
import argparse, pandas as pd, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--geno_csv', required=True, help='samples x SNPs matrix, columns=SNP ids')
ap.add_argument('--sumstats', required=True)
ap.add_argument('--p_thresh', type=float, default=5e-8)
ap.add_argument('--kb', type=int, default=250)
ap.add_argument('--r2', type=float, default=0.1)
ap.add_argument('--out_weights', required=True)
args = ap.parse_args()

G = pd.read_csv(args.geno_csv)
G = G.set_index(G.columns[0]) if G.columns[0] == 'sample_id' else G
S = pd.read_csv(args.sumstats, sep='\t')
S = S.sort_values('P')

selected = []
for _, row in S.iterrows():
    if row['P'] > args.p_thresh:
        break
    snp = row['SNP']
    if snp not in G.columns:
        continue
    keep = True
    for k in selected:
        # naive LD proxy: correlation between genotype vectors
        r = np.corrcoef(G[snp].values, G[k].values)[0,1]
        if r**2 > args.r2:
            keep = False; break
    if keep:
        selected.append(snp)

W = S[S['SNP'].isin(selected)][['SNP','A1','BETA']]
W.to_csv(args.out_weights, sep='\t', index=False)
print(f'Selected {len(W)} SNPs out of {len(S)} at P<{args.p_thresh}')
```

---

## `bin/wrap_prscs.py` (stub wrapper)

```python
#!/usr/bin/env python3
import argparse, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--sumstats', required=True)
ap.add_argument('--ldref', required=True)
ap.add_argument('--out', required=True)
args = ap.parse_args()

# In real runs, call PRS-CS python package / script here and convert to 3-col weights.
# For now, pass-through BETA as-is (placeholder) so pipeline runs end-to-end on toy data.
df = pd.read_csv(args.sumstats, sep='\t')
out = df[['SNP','A1','BETA']].copy()
out.to_csv(args.out, sep='\t', index=False)
print('PRS-CS stub: copied GWAS betas (replace with real PRS-CS call).')
```

---

## `bin/ldpred2.R` (skeleton)

```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table)
})

option_list <- list(
  make_option(c('--sumstats'), type='character'),
  make_option(c('--ldref'), type='character'),
  make_option(c('--out'), type='character')
)
opt <- parse_args(OptionParser(option_list=option_list))

ss <- fread(opt$sumstats)
# TODO: implement full LDpred2 using bigsnpr/bigstatsr with external LD ref.
# For toy flow, write passthrough weights to keep pipeline running.
out <- ss[, .(SNP, A1, BETA)]
fwrite(out, opt$out, sep='\t')
cat('LDpred2 stub: copied GWAS betas (replace with real LDpred2 implementation).\n')
```

---

## `README.md` (full)

```markdown
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
```

---

## Toy data format examples

### `resources/toy/sumstats.tsv`

```tsv
SNP	CHR	BP	A1	A2	BETA	SE	P
rs1	1	10100	A	G	0.12	0.05	1e-8
rs2	1	13500	C	T	-0.08	0.04	2e-6
rs3	2	70500	G	A	0.05	0.02	1e-3
```

### `resources/toy/phenotypes.tsv`

```tsv
sample_id	pheno
id1	1
id2	0
id3	1
```

---

## Notes on real implementations

* **PRS-CS**: supply ancestry-matched LD reference (e.g., 1000G EUR/AFR/...) and call the official script with appropriate `--phi` or auto mode; convert outputs to 3-column weights for scoring.
* **LDpred2**: use `bigsnpr` with external LD blocks; write beta-adjusted weights.
* **C+T grid**: run multi-threshold grid (e.g., p ∈ {5e-8, 1e-5, 1e-3, 0.05}, r² ∈ {0.1, 0.2}, kb ∈ {250, 500}) and pick best on validation set.

## Reproducibility

* Pin tool versions in Dockerfiles/Conda; capture `nextflow.log`; set seeds where applicable.

## License

MIT


