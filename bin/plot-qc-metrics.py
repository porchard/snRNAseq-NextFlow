#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from skimage.filters import threshold_multiotsu
import numpy as np
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--prefix', default='qc.', help='Prefix for output files (default: "qc.")')
parser.add_argument('rna_metrics', help='Path to file of RNA metrics')
parser.add_argument('dropkick', help='Path to file of dropkick output')
args = parser.parse_args()

RNAQC_METRIC_FILE = args.rna_metrics
DROPKICK = args.dropkick
PREFIX = args.prefix

@ticker.FuncFormatter
def read_count_formatter(x, pos):
    if x >= 1e9:
        return '{}B'.format(x/1e9)
    if x >= 1e6:
        return '{}M'.format(x/1e6)
    if x >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


qc = pd.read_csv(RNAQC_METRIC_FILE, sep='\t')
dropkick = pd.read_csv(DROPKICK, sep='\t')

bulk_stats_no_barcodes = qc.loc[qc.barcode=='-',['total_reads', 'uniquely_mapped_reads']].sum().rename('value').reset_index().assign(grp='Reads w/o\nwhitelisted barcode')
bulk_stats_whitelisted_barcodes = qc.loc[qc.barcode!='-',['total_reads', 'uniquely_mapped_reads']].sum().rename('value').reset_index().assign(grp='Reads w/\nwhitelisted barcode')
bulk_stats = pd.concat([bulk_stats_no_barcodes, bulk_stats_whitelisted_barcodes]).rename(columns={'index': 'stat'})
bulk_stats.stat = bulk_stats.stat.map({'total_reads': 'Total reads', 'uniquely_mapped_reads': 'Uniquely mapped reads'})

cumulative = qc[(qc.barcode!='-') & (qc.umis>0)].sort_values('umis', ascending=False)
cumulative['rnk'] = range(1, len(cumulative)+1)

qc = qc[qc.barcode!='-'].merge(dropkick, how='left')
qc.dropkick_score = qc.dropkick_score.fillna(0)

# try to infer UMI threshold
def estimate_threshold(x):
    # do on logscale
    values = np.log10(x).values
    values = values.reshape((len(values),1))
    thresholds = threshold_multiotsu(image=values, classes=3, nbins=256)
    # convert back to linear scale
    thresholds = [pow(10, i) for i in thresholds]
    UMI_THRESHOLD = round(thresholds[1])
    return UMI_THRESHOLD


MAX_EXPECTED_NUMBER_NUCLEI = int(10e3)
LOWERBOUNDS = [1, 5, 10, 50, 100, 500]
for i in LOWERBOUNDS:
    UMI_THRESHOLD = estimate_threshold(qc[(qc.barcode!='-') & (qc.umis>=i)].umis.astype(int))
    NUMBER_MEETING_UMI_THRESHOLD = (qc.umis>=UMI_THRESHOLD).sum()
    if NUMBER_MEETING_UMI_THRESHOLD <= MAX_EXPECTED_NUMBER_NUCLEI:
        break

if NUMBER_MEETING_UMI_THRESHOLD > MAX_EXPECTED_NUMBER_NUCLEI:
    # just fall back to 500
    UMI_THRESHOLD = 500
    NUMBER_MEETING_UMI_THRESHOLD = (qc.umis>=UMI_THRESHOLD).sum()


suggested_thresholds = pd.DataFrame({'metric': ['min_UMI'], 'threshold': UMI_THRESHOLD})
suggested_thresholds.to_csv(f'{PREFIX}suggested-thresholds.tsv', sep='\t', index=False)


fig, axs = plt.subplots(ncols=4, figsize=(5.5*4, 5.5))

# read counts
ax = axs[0]
sns.barplot(x='stat', y='value', hue='grp', ax=ax, data=bulk_stats)
ax.legend().set_title('')
ax.set_ylim(0, bulk_stats.value.max()*1.5)
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('')
ax.set_ylabel('Reads')
ax.grid(True)

# UMI rank plot
ax = axs[1]
sns.scatterplot(x='rnk', y='umis', data=cumulative, edgecolor=None, ax=ax)
ax.axhline(UMI_THRESHOLD, color='red', ls='--', label='Recommended\nmin. UMIs = {:,}'.format(UMI_THRESHOLD))
ax.axvline(NUMBER_MEETING_UMI_THRESHOLD, color='blue', ls='--', label='# nuclei meeting\nUMI threshold = {:,}'.format(NUMBER_MEETING_UMI_THRESHOLD))
ax.set_xscale('log')
ax.set_yscale('log')
ax.yaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Barcode rank')
ax.set_ylabel('# UMIs')
ax.grid(True)
ax.legend()

# UMIs vs fraction mitochondrial
ax = axs[2]
sns.scatterplot(x='umis', y='fraction_mitochondrial', hue='dropkick_score', data=qc, alpha=0.01, ax=ax)
ax.set_xscale('log')
ax.set_xlim(left=1)
ax.axvline(UMI_THRESHOLD, color='red', ls='--')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('# UMIs')
ax.set_ylabel('Fraction mitochondrial')
ax.grid(True)

ax = axs[3]
sns.scatterplot(x='umis', y='fraction_mitochondrial', hue='dropkick_score', data=qc, alpha=0.01, ax=ax)
ax.set_xscale('log')
ax.set_xlim(left=min([100, UMI_THRESHOLD*0.8]))
ax.axvline(UMI_THRESHOLD, color='red', ls='--')
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('# UMIs')
ax.set_ylabel('Fraction mitochondrial')
ax.grid(True)

fig.tight_layout()
fig.savefig(f'{PREFIX}metrics.png')
fig.clf()