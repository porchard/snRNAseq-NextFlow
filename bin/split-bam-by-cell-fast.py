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
import copy

parser = argparse.ArgumentParser()
parser.add_argument('--max-reads-in-memory', dest = 'max_reads', type = int, default = 1000000, help = 'Hold up to this many aligned reads in memory before writing to disk.')
parser.add_argument('--min-reads-to-output', dest = 'min_cell_reads', type = int, default = 1, help = 'Ignore cells with fewer than this many reads.')
parser.add_argument('--keep-barcodes', dest = 'keep', type = str, default = '', help = 'Keep the barcodes listed in this file. In this case, --min-reads-to-output will be ignored')
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('--prefix', default = '', help = 'Prefix of output files {PREFIX}{cell}{SUFFIX}.bam')
parser.add_argument('--suffix', default = '', help = 'Suffix of output files {PREFIX}{cell}{SUFFIX}.bam')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')
logging.info('Splitting cells from file {}'.format(args.bam_in))


keep_barcodes = set()
if args.keep != '':
    with open(args.keep, 'r') as f:
        for line in f:
            keep_barcodes.add(line.rstrip())


def readgroup_to_cell (readgroup):
    flowcell, lane, cell = readgroup.split('.')
    return cell


def readname_to_cell(readname):
    name, cell, umi = readname.split('_')
    return cell


def make_sam_name (cell):
    return '{}{}{}.sam'.format(args.prefix, cell, args.suffix)


def append_sam(cell, header, reads):
    sam = make_sam_name(cell)
    if not os.path.exists(sam):
        with pysam.AlignmentFile(sam, 'w', header = header) as f:
            pass
    with open(sam, 'a+') as f:
        # just append the reads...
       for read in reads:
           f.write(read + '\n')
    return True


# fetch the header
logging.info('Reading header from {}'.format(args.bam_in))
f = pysam.AlignmentFile(args.bam_in, 'rb')
if not isinstance(f.header, dict):
    header = f.header.as_dict()
f.close()
logging.info('Finished reading header from {}'.format(args.bam_in))


keep = set()
if args.keep != '':
    logging.info(f'Reading {args.keep}; --min-reads-to-output argument being ignored.')
    with open(args.keep, 'r') as f:
        for line in f:
            keep.add(line.rstrip())
    logging.info('Read {} barcodes to keep'.format(len(keep)))
else:
    logging.info('Counting reads per barcode')
    # get the read counts for each library
    cell_read_counts = Counter()
    with pysam.AlignmentFile(args.bam_in, 'rb') as f:
        total_read_count = 0 # number of reads encountered
        for read in f.fetch(until_eof=True):
            total_read_count += 1
            if total_read_count % 1000000 == 0:
                logging.info('Read through {} reads'.format(total_read_count))
            cell = readname_to_cell(read.query_name)
            cell_read_counts[cell] += 1
            
    logging.info('Identified {} libraries with a total of {} reads'.format(len(cell_read_counts), sum(cell_read_counts.values())))
    # remove cells with fewer than the threshold number of reads...
    keep = set([cell for cell, count in cell_read_counts.items() if count >= args.min_cell_reads])
    drop = [cell for cell, count in cell_read_counts.items() if count < args.min_cell_reads]
    logging.info('Keeping {} libraries with at least {} reads'.format(len(keep), args.min_cell_reads))

# filter through the file, writing text files and counting the reads per cell
logging.info('Collecting reads from libraries to keep...')
reads = dict()
with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    total_read_count = 0 # number of reads encountered
    for read in f.fetch(until_eof=True):
        total_read_count += 1
        if total_read_count % 1000000 == 0:
            logging.info('Read through {} reads'.format(total_read_count))
        #cell = readgroup_to_cell(read.get_tag('RG'))
        cell = readname_to_cell(read.query_name)
        if cell not in keep:
            continue
        if cell not in reads:
            reads[cell] = []
        reads[cell].append(read.to_string())
        if total_read_count % args.max_reads == 0:
            logging.info('Writing {} reads to sam files...'.format(args.max_reads))
            for cell in reads:
                append_sam(cell, header, reads[cell])
            reads = {}

# now do one final write...
logging.info('Writing {} reads to sam files...'.format(args.max_reads))
for cell in reads:
    append_sam(cell, header, reads[cell])
reads = dict()

logging.info('Finished.')

