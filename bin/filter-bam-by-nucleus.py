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

# samtools view or samtools split are quite inefficient 
# store all reads in a dict, and then write them separately...?
# or could store e.g. 1000000 reads at a time, and iteratively write the files...
# that would be less memory intensive and probably just as fast

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Input bam file')
parser.add_argument('nuclei', help = 'List of nuclei to keep')
parser.add_argument('bam_out', help = 'Output bam file')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

logging.info('Loading list of nuclei to keep')
keep = set()

def get_nucleus(read):
    return read.query_name.split('_')[1]

with open(args.nuclei, 'r') as f:
    for line in f:
        keep.add(line.rstrip())

logging.info('Writing new bam')

with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    header = f.header
    with pysam.AlignmentFile(args.bam_out, 'wb', header=header) as f_new:
        for read in f.fetch(until_eof=True):
            if get_nucleus(read) in keep:
                f_new.write(read)

logging.info('Done.')

