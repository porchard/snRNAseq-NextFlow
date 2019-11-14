#!/usr/bin/env python
# coding: utf-8

import sys
import os
from matplotlib import pyplot as plt
import seaborn as sns
from QoRTs import QoRTs
import glob
import re
import pandas as pd


LIST_OF_QORTS_DIRS = sys.argv[1]
OUT_PREFIX = sys.argv[2]


QORTS_DIRS = []
with open(LIST_OF_QORTS_DIRS, 'r') as f:
    QORTS_DIRS = [line.rstrip() for line in f]

q = [QoRTs(directory=d) for d in QORTS_DIRS]


def dir_to_nucleus_info(d):
    library, genome, barcode = re.match('^(.*)-(.*)-(.*)$', os.path.basename(d)).groups()
    return locals()

for i in range(0, len(QORTS_DIRS)):
    q[i].library = '{library}-{barcode}'.format(**dir_to_nucleus_info(QORTS_DIRS[i]))

# will want to extract:
# chromosome counts
# mapping location
# strandedness
# gene counts
# tmp.number_genes_with_nonzero_counts
# read counts

# chromosome_counts
chromosome_counts = pd.concat([pd.DataFrame([chrom, count] for chrom, count in i.chrom_count.items()).assign(library=i.library) for i in q])
chromosome_counts.columns = ['chrom', 'count', 'library']
chromosome_counts.to_csv(f'{OUT_PREFIX}chrom_counts.txt', index = False)

# mapping location
mapping_location = pd.concat([pd.DataFrame([loc, count] for loc, count in i.mapping_location.items()).assign(library=i.library) for i in q])
mapping_location.columns = ['loc', 'count', 'library']
mapping_location.to_csv(f'{OUT_PREFIX}mapping_location.txt', index = False)

# strandedness
strandedness = pd.concat([pd.DataFrame([loc, count] for loc, count in i.strand_test.items()).assign(library=i.library) for i in q])
strandedness.columns = ['loc', 'count', 'library']
strandedness.to_csv(f'{OUT_PREFIX}strandedness.txt', index = False)

# read pair counts
read_pair_counts = pd.DataFrame([[i.total_read_pairs, i.quality_read_pairs, i.library] for i in q])
read_pair_counts.columns = ['total_read_pairs', 'quality_read_pairs', 'library']
read_pair_counts.to_csv(f'{OUT_PREFIX}read_pair_counts.txt', index = False)

# gene counts
gene_counts = pd.concat([pd.DataFrame([gene, count] for gene, count in i.gene_count.items() if count != 0).assign(library=i.library) for i in q])
gene_counts.columns = ['gene', 'count', 'library']
gene_counts.to_csv(f'{OUT_PREFIX}gene_counts.txt', index = False)

# genes detected
genes_detected = pd.DataFrame([[i.number_genes_with_nonzero_counts, i.number_genes_with_zero_counts, i.library] for i in q])
genes_detected.columns = ['genes_with_nonzero_counts', 'genes_with_zero_counts', 'library']
genes_detected.to_csv(f'{OUT_PREFIX}genes_detected.txt', index = False)

gene_body_coverage = pd.concat([pd.DataFrame([block, pct['total'], pct['bottom_half'], pct['upper_mid_quartile'], pct['75_to_90'], pct['high']] for block, pct in i.gene_body_coverage_pct.items()).assign(library=i.library) for i in q])
gene_body_coverage.columns = ['block', 'total', 'bottom_half', 'upper_mid_quartile', '75_to_90', 'high', 'library']
gene_body_coverage.to_csv(f'{OUT_PREFIX}gene_body_coverage.txt', index = False)
