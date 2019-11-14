#!/usr/bin/env python
import os
import sys
import re
import logging
import argparse 

parser = argparse.ArgumentParser(description='Assign barcodes to genomes.', add_help = True)
parser.add_argument('--genomes', type = str, nargs = '+', help = 'List of genomes, in the same order as the respective count files')
parser.add_argument('--counts', type = str, nargs = '+', help = 'List of count files')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

COUNT_FILES = dict(zip(args.genomes, args.counts))

counts = dict() # barcode --> genome --> count
genomes = list(COUNT_FILES.keys())

logging.info('Identified genomes: {}'.format(', '.join(list(genomes))))

for g, count_file in COUNT_FILES.items():
    logging.info(f'Reading barcode counts from file {count_file}')
    with open(count_file, 'r') as f:
        for line in f:
            count, barcode = line.rstrip().split()
            count = int(count)
            if barcode not in counts:
                counts[barcode] = {genome: 0 for genome in genomes}
            counts[barcode][g] += count

logging.info('Assigning barcodes to genomes (threshold: 0.6)')
THRESHOLD = 0.6

assignments = {barcode: None for barcode in counts}

for barcode, c in counts.items():
    total = sum(c.values())
    for genome, i in c.items():
        if float(i) / total >= THRESHOLD:
            assignments[barcode] = genome
            continue

for barcode, assignment in assignments.items():
    print(f'{barcode}\t{assignment}')
