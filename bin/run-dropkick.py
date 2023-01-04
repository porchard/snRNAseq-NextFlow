#!/usr/bin/env python
# coding: utf-8

import sys
import scanpy
import dropkick as dk
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pickle

INDIR, PREFIX = sys.argv[1:]

x = scanpy.read_10x_mtx(INDIR)

# TODO: may wish to select min_genes based on the overall distribution or the value for e.g. the 100kth barcode.
adata = dk.recipe_dropkick(x, X_final="raw_counts", mito_names='^mt-|^MT-|^Mt-')

qc_plt = dk.qc_summary(adata)
qc_plt.savefig(f'{PREFIX}dropkick-qc.png', dpi=300)
qc_plt.clf()

adata_model = dk.dropkick(adata, n_jobs=5, mito_names='^mt-|^MT-|^Mt-')


with open(f'{PREFIX}adata.pickle', 'wb') as fh:
    pickle.dump(adata, fh)
with open(f'{PREFIX}adata_model.pickle', 'wb') as fh:
    pickle.dump(adata_model, fh)


score_plt = dk.score_plot(adata)
score_plt.savefig(f'{PREFIX}dropkick-scoreplot.png', dpi=300)
plt.clf()

adata.obs['dropkick_score'].to_frame().rename_axis(index='barcode').to_csv(f'{PREFIX}dropkick-score.tsv', sep='\t')

fig, ax = plt.subplots()
adata.obs['dropkick_score'].hist(bins=50)
fig.tight_layout()
fig.savefig(f'{PREFIX}dropkick-score-distribution.png', dpi=300)
fig.clf()