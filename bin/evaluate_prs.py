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

