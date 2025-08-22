nextflow.enable.dsl=2

// Channels (fail fast if a path is wrong)
Channel
  .fromPath(params.sumstats, checkIfExists: true)
  .set { CH_SUMSTATS }

def GENO_CHANNEL = (
  params.vcf
    ? Channel.fromPath(params.vcf, checkIfExists: true)
    : Channel.fromPath(params.geno_prefix, checkIfExists: true)
).broadcast()

def METHODS = (params.methods instanceof String) ? params.methods.split(',') : params.methods


process CT_CLUMP {
  publishDir "${params.outdir}/ct", mode: 'copy'
  input:
    path sumstats from CH_SUMSTATS
    path geno     from GENO_CHANNEL      // <-- was `val geno`
  output:
    path 'ct.snplist' into CH_CT_SNPS
  script:
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    plink2 ${genoArg} --clump ${sumstats} --clump-field P \
      --clump-p1 5e-8 --clump-p2 1e-2 --clump-r2 0.1 --clump-kb 250 \
      --out clumped
    cut -f3 clumped.clumps | tail -n +2 > ct.snplist
    """
}

process CT_SCORE {
  publishDir "${params.outdir}/ct", mode: 'copy'
  input:
    path sumstats from CH_SUMSTATS
    path geno     from GENO_CHANNEL      // <-- was `val geno`
    path snplist  from CH_CT_SNPS
  output:
    path 'ct.profile'
    path 'weights.txt' into CH_WEIGHTS_CT   // <-- expose weights to next step
  script:
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    awk 'NR==FNR{a[\$1]=1; next} NR>1 && a[\$1]{print \$1"\\t"\$4"\\t"\$6}' \
        ${snplist} ${sumstats} > weights.txt
    plink2 ${genoArg} --score weights.txt 1 2 3 header-read cols=scoresums --out ct
    mv ct.sscore ct.profile
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

process PRSCS_RUN {
  // Only runs if you include "prscs" in --methods; otherwise it is skipped.
  when: METHODS.contains('prscs')

  publishDir "${params.outdir}/prscs", mode: 'copy'
  container params.py_image   // remove this line if you don't define params.py_image

  input:
    path sumstats from CH_SUMSTATS

  output:
    path 'prscs.weights.tsv'

  script:
    """
    # PRS-CS stub wrapper: converts GWAS betas -> 3-col weights
    wrap_prscs.py \
      --sumstats ${sumstats} \
      --ldref ${params.ldref} \
      --out prscs.weights.tsv
    """
}



process SCORE_ALL {
  publishDir "${params.outdir}/scores", mode: 'copy'
  input:
    path pheno from Channel.fromPath(params.pheno, checkIfExists: true)
    tuple path(weights), val(method) from (
      CH_WEIGHTS_CT.map{ w -> tuple(w, 'ct') }
      // add PRS-CS/LDpred2 here later when you enable them
    )
    path geno from GENO_CHANNEL
  output:
    path 'scores.tsv'
  script:
    def genoArg = params.vcf ? "--vcf ${geno}" : "--bfile ${geno}"
    """
    python3 /bin/score_prs.py \\
      --weights ${weights} --method ${method} \\
      ${params.vcf ? '--vcf ' + geno : '--bfile ' + geno} \\
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
