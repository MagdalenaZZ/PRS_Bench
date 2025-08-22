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

