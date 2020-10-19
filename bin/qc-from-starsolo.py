#!/usr/bin/env python

import os
import sys
import pysam
import json
import argparse
import logging
import numpy as np
from scipy.io import mmread

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Gather stats from a snRNA-seq bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file (output by starsolo).')
parser.add_argument('count_matrix', type = str,  help = 'Count matrix, in market matrix format.')
parser.add_argument('barcodes_tsv', type = str,  help = 'List of the barcodes to interpret market matrix format.')
parser.add_argument('--cell-tag', dest='cell_tag', type = str, default = 'CB', help = 'Tag denoting the cell/nucleus (default: CB)')
parser.add_argument('--min-reads', dest='min_reads', type = int, default = 0, help = 'Suppress output for cells with fewer than this many reads (default: 0).')
args = parser.parse_args()

CELL_TAG = args.cell_tag
MIN_READS_FOR_OUTPUT = args.min_reads


class Cell:

    def __init__(self, barcode):
        self.barcode = barcode
        self.supplementary_alignments = 0
        self.secondary_alignments = 0
        self.total_reads = 0 # excludes secondary/supplementary alignments
        self.uniquely_mapped_reads = 0
        self.chromosome_read_counts = dict() # chrom -> count


    def record_alignment(self, read):
        if read.is_secondary or read.is_supplementary:
            if read.is_secondary:
                self.secondary_alignments += 1
            if read.is_supplementary:
                self.supplementary_alignments += 1
            return 0
        self.total_reads += 1
        if read.mapping_quality == 255:
            self.uniquely_mapped_reads += 1
            chrom = read.reference_name
            if chrom not in self.chromosome_read_counts:
                self.chromosome_read_counts[chrom] = 0
            self.chromosome_read_counts[chrom] += 1
        return 0


    def gather_metrics(self):
        metrics = dict()
        metrics['barcode'] = self.barcode
        metrics['total_reads'] = self.total_reads
        metrics['uniquely_mapped_reads'] = self.uniquely_mapped_reads
        metrics['secondary_alignments'] = self.secondary_alignments
        metrics['supplementary_alignments'] = self.supplementary_alignments
        metrics['fraction_mitochondrial'] = self.chromosome_read_counts['chrM'] / self.uniquely_mapped_reads if 'chrM' in self.chromosome_read_counts else 0
        return metrics



def get_umis_per_barcode(mtx, barcodes_tsv):
    barcodes = []
    with open(barcodes_tsv, 'r') as f:
        for line in f:
            barcodes.append(line.rstrip())
    
    mat = mmread(mtx)
    umis_per_barcode = mat.sum(axis=0)
    return dict(zip(barcodes, umis_per_barcode.A.flatten()))




if __name__ == '__main__':
    cells = dict()
    no_cell_tag = 0
    total_reads = 0

    logging.info('Using count matrix to get final number of UMIs per barcode')
    umis_per_barcode = get_umis_per_barcode(args.count_matrix, args.barcodes_tsv)
    
    logging.info('Reading bam file')
    with pysam.AlignmentFile(args.bam, 'rb') as f:
        for read in f.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {} reads'.format(total_reads))
            barcode = 'no_barcode' if not read.has_tag(CELL_TAG) else read.get_tag(CELL_TAG)
            if barcode not in cells:
                cells[barcode] = Cell(barcode)
            cells[barcode].record_alignment(read)

    logging.info('Finished reading bam file.')
    logging.info('Outputting metrics')

    print_metrics = ['barcode', 'total_reads', 'uniquely_mapped_reads', 'secondary_alignments', 'supplementary_alignments', 'umis', 'fraction_mitochondrial']
    print('\t'.join(print_metrics))
    for cell in cells.values():
        metrics = cell.gather_metrics()
        metrics['umis'] = umis_per_barcode[metrics['barcode']] if metrics['barcode'] in umis_per_barcode else None
        to_print = [str(metrics[i]) for i in print_metrics]
        print('\t'.join(to_print))

    logging.info('Done')
