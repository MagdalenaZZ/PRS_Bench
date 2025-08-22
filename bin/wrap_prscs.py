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


