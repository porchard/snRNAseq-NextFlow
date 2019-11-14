#!/usr/bin/env python

import sys
import logging
import copy
import os
import json
import gzip
import numpy
import argparse
import re
import pysam
import random
from collections import Counter

parser = argparse.ArgumentParser(description='Make a file of nucleus, gene, count from a deduped bam')
parser.add_argument('bam', help = 'Input bam file')
parser.add_argument('--nucleus-tag', dest = 'nucleus_tag', default = 'CB', help = 'Tag that indicates nucleus (default: CB)')
parser.add_argument('--gene-tag', dest = 'gene_tag', default = 'XT', help = 'Tag that indicates which gene the read has been assigned to (default: XT)')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

logging.info(f'Reading {args.bam}')

counts = dict() # cell -> gene -> count
count = 0

with pysam.AlignmentFile(args.bam, 'rb') as f:
    for read in f.fetch(until_eof=True):
        count += 1
        if count % 1000000 == 0:
            logging.info(f'Processed {count} reads')
        if not read.has_tag(args.nucleus_tag):
            logging.warning(f'Skipping read {read.query_name} (no tag {args.nucleus_tag})')
            continue
        nucleus = read.get_tag(args.nucleus_tag)
        if nucleus not in counts:
            counts[nucleus] = {}
        gene = read.get_tag(args.gene_tag) if read.has_tag(args.gene_tag) else None
        if not read.has_tag(args.gene_tag):
            logging.warning(f'Skipping read {read.query_name} (no tag {args.gene_tag})')
            continue
        gene = read.get_tag(args.gene_tag)
        if gene not in counts[nucleus]:
            counts[nucleus][gene] = 0
        counts[nucleus][gene] += 1

logging.info('Finished reading through bam file')

for nucleus, stats in counts.items():
    for gene, count in stats.items():
        print(f'{nucleus}\t{gene}\t{count}')

logging.info('Done.')

