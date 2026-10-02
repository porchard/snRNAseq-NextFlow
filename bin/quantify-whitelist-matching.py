#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
import logging

import pandas as pd
import pysam

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')



parser = argparse.ArgumentParser(description='Gather stats from a snRNA-seq bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file (output by starsolo).')
parser.add_argument('whitelist', type = str,  help = 'Barcode whitelist..')
args = parser.parse_args()


logging.info('Reading whitelist')
whitelist = set(pd.read_csv(args.whitelist, header=None)[0].to_list())

total_reads = 0
primary_reads = 0
cr_match = 0 # primary only
cb_match = 0 # primary only
cr_equals_cb = 0

logging.info('Reading bam file')
with pysam.AlignmentFile(args.bam, 'rb') as f:
    for read in f.fetch(until_eof=True):
        total_reads += 1
        if total_reads % 1000000 == 0:
            logging.info('Processed {} reads'.format(total_reads))
        if read.is_secondary or read.is_supplementary:
            continue
        primary_reads += 1
        cr = read.get_tag('CR') if read.has_tag('CR') else None
        cb = read.get_tag('CB') if read.has_tag('CB') else None
        if cr in whitelist:
            cr_match += 1
        if cb in whitelist:
            cb_match += 1
        if cr is not None and cb is not None and cr == cb:
            cr_equals_cb += 1
            
logging.info('Finished reading bam file.')

df = pd.DataFrame([[total_reads, primary_reads, cr_match, cb_match, cr_equals_cb]], columns=['total_reads', 'primary_read', 'cr_match', 'cb_match', 'cr_equals_cb'])
df.to_csv(sys.stdout, sep='\t', index=False)

logging.info('Done')

