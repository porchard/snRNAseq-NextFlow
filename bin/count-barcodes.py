#!/usr/bin/env python
import os
import sys
import argparse
import gzip
from Bio import SeqIO
import re
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description="Count barcodes from an index fastq files")
parser.add_argument('fastq_files', action='store', nargs='+', type=str, help = 'Path to the fastq files from which barcodes should be counted')
args = parser.parse_args()

barcode_counts = dict()

for fastq_file in args.fastq_files:
    logging.info(f'Processing fastq files {fastq_file}')
    count = 0
    with gzip.open(fastq_file, 'rt') as fastq:
        for record in SeqIO.parse(fastq, 'fastq'):
            count += 1
            if count % 1000000 == 0:
                logging.info(f'Processed {count} reads')
            barcode = str(record.seq)
            if barcode not in barcode_counts:
                barcode_counts[barcode] = 0
            barcode_counts[barcode] += 1

for barcode, count in barcode_counts.items():
    print(f'{barcode}\t{count}')

logging.info('Done')

