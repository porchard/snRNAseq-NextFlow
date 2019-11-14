#!/usr/bin/env python
import os
import sys
import argparse
import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description="Fastq file from which to cut out UMI and barcode sequences.")
parser.add_argument('fastq_file', action='store', type=str, help = 'Path to the fastq file from which UMI and barcode sequences should be extracted.')
parser.add_argument('barcode_fastq_file', action='store', type=str, help = 'File of barcodes to write.')
parser.add_argument('umi_fastq_file', action='store', type=str, help = 'File of UMIs to write.')
args = parser.parse_args()

BARCODE_START = 0
BARCODE_LENGTH = 16
BARCODE_END = BARCODE_START + BARCODE_LENGTH
UMI_START = 16
UMI_LENGTH = 10
UMI_END = UMI_START + UMI_LENGTH


with gzip.open(args.fastq_file, 'rt') as fastq_in:
    with gzip.open(args.barcode_fastq_file, 'wt') as barcode_fastq:
        with gzip.open(args.umi_fastq_file, 'wt') as umi_fastq:
            count = 0
            for title, seq, qual in FastqGeneralIterator(fastq_in):
                barcode_fastq.write('@{}\n{}\n+\n{}\n'.format(title, seq[BARCODE_START:BARCODE_END], qual[BARCODE_START:BARCODE_END]))
                umi_fastq.write('@{}\n{}\n+\n{}\n'.format(title, seq[UMI_START:UMI_END], qual[UMI_START:UMI_END]))
                count += 1
                if count % 1000000 == 0:
                    logging.info(f'Processed {count} reads')

logging.info('Done')

