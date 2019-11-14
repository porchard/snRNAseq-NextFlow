#!/usr/bin/env python
import os
import sys
import argparse
import gzip
import re
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description="Filter a fastq file to reads from good_barcodes")
parser.add_argument('--min-counts', dest = 'min_counts', action='store', type=int, default = 1, help = 'Whitelist will consist of barcodes in the good_barcodes list that appear at least this many times.')
parser.add_argument('good_barcodes', action='store', type=str, help = 'Path to 10X whitelist')
parser.add_argument('barcode_count_files', action='store', nargs='+', type=str, help = 'Path to the barcode count files')
args = parser.parse_args()

logging.info('Reading 10X whitelist')
barcodes = dict()
with open(args.good_barcodes, 'r') as f:
    for line in f:
        line = line.rstrip()
        barcodes[line] = 0

logging.info('Finished reading 10X whitelist')

for barcode_count_file in args.barcode_count_files:
    logging.info(f'Processing file {barcode_count_file}')
    with open(barcode_count_file, 'r') as f:
        for line in f:
            barcode, count = line.rstrip().split()
            count = int(count)
            if barcode in barcodes:
                barcodes[barcode] += count

keep = [barcode for barcode, count in barcodes.items() if count >= args.min_counts]

logging.info('Keeping {} barcodes for whitelist'.format(len(keep)))
logging.info('Outputting whitelist')

for i in keep:
    print(i)

logging.info('Done')

