#!/usr/bin/env python

import os
import sys
import pybedtools as bt
import pysam
import logging
import argparse

parser = argparse.ArgumentParser(description='Given a bam file, create scaled and stranded bigwigs.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file, sorted by read name.')
parser.add_argument('prefix', type = str, help = 'Prefix for output bedgraph files ({prefix}.fwd.bdg; {prefix}.rev.bdg)')
parser.add_argument('--scale-per-reads', dest='scale_per_reads', type = int, default = 1000000, help = 'Scale per this many reads (default: 1000000, i.e. signal per million reads)')
parser.add_argument('--no-scale', dest='no_scale', default = False, action='store_true', help = 'Do not scale the bedgraph file (in this case, --scale-per-reads is ignored).')
parser.add_argument('--rev-strand-negative', dest='rev_strand_negative', default = False, action='store_true', help = 'Make the reverse strand values negative.')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


# Given a bam file (probably RNA-seq):
# 1. Determine read counts
# 2. Calculate scaling factor
# 3. Use bedtools genomecov to create bedgraph files


read_count = 0

if not args.no_scale:
    logging.info('Counting reads in bam file')

    with pysam.AlignmentFile(args.bam, 'rb') as f:
        for read in f.fetch(until_eof=True):
            read_count += 1
            if read_count % 1000000 == 0:
                logging.info('Counted {} reads so far'.format(read_count))

    logging.info('Counted {} reads'.format(read_count))

scale_factor = 1 if args.no_scale else args.scale_per_reads / read_count
logging.info('Creating bedgraphs with scaling factor {}'.format(scale_factor))
logging.info('Creating bedgraph for + strand')
bt.BedTool(args.bam).genome_coverage(bg=True, split=True, strand='+', scale=scale_factor).saveas('{}.fwd.bdg'.format(args.prefix))
logging.info('Creating bedgraph for - strand')
sign = 1 if not args.rev_strand_negative else -1
bt.BedTool(args.bam).genome_coverage(bg=True, split=True, strand='-', scale=sign*scale_factor).saveas('{}.rev.bdg'.format(args.prefix))
logging.info('Done')
