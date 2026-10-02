#!/usr/bin/env python
# coding: utf-8

import os
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Plot whitelist matching stats (output of quantify-whitelist-matching.py) across libraries.')
parser.add_argument('--tsv', nargs='+', help='Paths to the TSV files (library name = filename minus .whitelist-matching.tsv)')
parser.add_argument('--out', help='Output figure file')
args = parser.parse_args()

TSVS = args.tsv
OUT = args.out

FRACTION_OF_PRIMARY = ['cr_match', 'cb_match', 'cr_equals_cb']
XLABELS = {
    'total_reads': 'Total alignments (millions)',
    'primary_read': 'Primary alignments (millions)',
    'cr_match': 'CR on whitelist\n(fraction of primary alignments)',
    'cb_match': 'CB on whitelist\n(fraction of primary alignments)',
    'cr_equals_cb': 'CR == CB\n(fraction of primary alignments)',
}

dfs = []
for f in TSVS:
    tmp = pd.read_csv(f, sep='\t')
    tmp['library'] = os.path.basename(f).removesuffix('.whitelist-matching.tsv')
    dfs.append(tmp)
df = pd.concat(dfs).sort_values('library').reset_index(drop=True)

columns = [c for c in df.columns if c != 'library']
for c in FRACTION_OF_PRIMARY:
    df[c] = df[c] / df['primary_read']

fig, axs = plt.subplots(ncols=len(columns), figsize=(3 * len(columns), 0.4 * len(df) + 1.5), sharey=True)
for ax, c in zip(axs, columns):
    sns.barplot(data=df, x=c, y='library', color='#4c72b0', ax=ax)
    ax.set_ylabel('Library' if ax is axs[0] else '')
    ax.set_xlabel(XLABELS.get(c, c))
    if c in FRACTION_OF_PRIMARY:
        ax.set_xlim(0, 1.2)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.bar_label(ax.containers[0], fmt='%.3f', padding=2, fontsize='small')
    else:
        ax.bar_label(ax.containers[0], labels=[f'{v / 1e6:.1f}M' for v in df[c]], padding=2, fontsize='small')
        ax.set_xlim(0, df[c].max() * 1.3)
        ax.xaxis.set_major_formatter(lambda x, pos: f'{x / 1e6:g}')
fig.tight_layout()
fig.savefig(OUT, bbox_inches='tight', dpi=300)
