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
