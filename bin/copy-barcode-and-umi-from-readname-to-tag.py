#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import re
import pysam
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out', help = 'Bam file to write')
args = parser.parse_args()

logging.info(f'Writing new bam file ({args.bam_out})')
with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    with pysam.AlignmentFile(args.bam_out, 'wb', template = f) as f_new:
        count = 0
        for read in f.fetch(until_eof=True):
            count += 1
            if count % 1000000 == 0:
                logging.info(f'Processed {count} reads so far')
            readname, barcode, umi = read.query_name.split('_')
            read.set_tag('CB', barcode)
            read.set_tag('LB', barcode)
            read.set_tag('UR', umi)
            f_new.write(read)

logging.info('Done')
