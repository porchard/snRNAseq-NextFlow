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


def mean(d):
    """Return the mean of values represented as {value: count}.

    >>> mean({1: 2, 3: 1})
    1.6666666666666667
    >>> mean({5: 10})
    5.0
    >>> mean({0: 1, 10: 1})
    5.0
    """
    total = sum(v * c for v, c in d.items())
    n = sum(d.values())
    return total / n


def median(d):
    """Return the median of values represented as {value: count}.

    >>> median({1: 2, 3: 1})
    1
    >>> median({1: 1, 2: 1, 3: 1, 4: 1})
    2.5
    >>> median({5: 10})
    5
    >>> median({1: 3, 2: 1, 3: 1})
    1
    """
    sorted_values = sorted(d.keys())
    n = sum(d.values())
    cumulative = 0
    mid = n / 2

    if n % 2 == 1:
        target = (n + 1) // 2
        for v in sorted_values:
            cumulative += d[v]
            if cumulative >= target:
                return v
    else:
        lower = None
        for v in sorted_values:
            cumulative += d[v]
            if lower is None and cumulative >= mid:
                lower = v
            if cumulative >= mid + 1:
                return (lower + v) / 2 if lower != v else v


parser = argparse.ArgumentParser(description='Gather stats from a snRNA-seq bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file (output by starsolo).')
parser.add_argument('genefull_exonoverintron_count_matrix', type = str,  help = 'STARsolo GeneFull_ExonOverIntron count matrix.')
parser.add_argument('gene_count_matrix', type = str,  help = 'STARsolo Gene count matrix.')
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
        # everything from here on down is primary alignments only
        self.primary_alignments = 0
        self.mapped = 0
        self.uniquely_mapped = 0
        self.chromosome_read_counts = dict() # chrom -> count
        self.mapq = dict() # mapq -> count
        self.assigned_to_gene = 0


    def record_alignment(self, read):
        if read.is_secondary or read.is_supplementary:
            if read.is_secondary:
                self.secondary_alignments += 1
            if read.is_supplementary:
                self.supplementary_alignments += 1
            return 0
        self.primary_alignments += 1
        if not read.is_unmapped:
            self.mapped += 1
        if read.mapping_quality == 255:
            self.uniquely_mapped += 1
            chrom = read.reference_name
            if chrom not in self.chromosome_read_counts:
                self.chromosome_read_counts[chrom] = 0
            self.chromosome_read_counts[chrom] += 1
        if read.mapping_quality not in self.mapq:
            self.mapq[read.mapping_quality] = 0
        self.mapq[read.mapping_quality] += 1
        if read.has_tag('GX') and read.get_tag('GX') != '-':
            self.assigned_to_gene += 1
        return 0


    def gather_metrics(self):
        metrics = dict()
        metrics['barcode'] = self.barcode
        metrics['secondary_alignments'] = self.secondary_alignments
        metrics['supplementary_alignments'] = self.supplementary_alignments
        metrics['primary_alignments'] = self.primary_alignments
        metrics['mapped_primary_alignments'] = self.mapped
        metrics['fraction_primary_alignments_mapped'] = self.mapped / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['uniquely_mapped_primary_alignments'] = self.uniquely_mapped
        metrics['fraction_primary_alignments_uniquely_mapped'] = self.uniquely_mapped / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['primary_alignments_assigned_to_gene'] = self.assigned_to_gene
        metrics['fraction_primary_alignments_assigned_to_gene'] = self.assigned_to_gene / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['mean_mapq'] = mean(self.mapq)
        metrics['median_mapq'] = median(self.mapq)
        metrics['primary_alignments_with_mapq_0'] = self.mapq.get(0, 0)
        metrics['fraction_primary_alignments_with_mapq_0'] = self.mapq.get(0, 0) / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['fraction_mitochondrial'] = self.chromosome_read_counts['chrM'] / self.uniquely_mapped if 'chrM' in self.chromosome_read_counts else 0
        return metrics



def get_umis_per_barcode(mtx):
    barcodes = []
    barcodes_tsv = os.path.join(os.path.dirname(mtx), 'barcodes.tsv')
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
    genefull_exonoverintron_umis_per_barcode = get_umis_per_barcode(args.genefull_exonoverintron_count_matrix)
    gene_umis_per_barcode = get_umis_per_barcode(args.gene_count_matrix)
    assert(list(gene_umis_per_barcode.keys()) == list(genefull_exonoverintron_umis_per_barcode.keys()))
    exon_full_gene_body_ratio = {barcode: gene_umis_per_barcode[barcode] / genefull_exonoverintron_umis_per_barcode[barcode] if genefull_exonoverintron_umis_per_barcode[barcode] > 0 else 0 for barcode in gene_umis_per_barcode}

    
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

    print_metrics = [
        'barcode',
        'umis',
        'exon_to_full_gene_body_ratio',
        'secondary_alignments',
        'supplementary_alignments',
        'primary_alignments',
        'mapped_primary_alignments',
        'fraction_primary_alignments_mapped',
        'uniquely_mapped_primary_alignments',
        'fraction_primary_alignments_uniquely_mapped',
        'primary_alignments_assigned_to_gene',
        'fraction_primary_alignments_assigned_to_gene',
        'mean_mapq',
        'median_mapq',
        'primary_alignments_with_mapq_0',
        'fraction_primary_alignments_with_mapq_0',
        'fraction_mitochondrial'
    ]

    print('\t'.join(print_metrics))
    for cell in cells.values():
        metrics = cell.gather_metrics()
        metrics['umis'] = genefull_exonoverintron_umis_per_barcode[metrics['barcode']] if metrics['barcode'] in genefull_exonoverintron_umis_per_barcode else 'NA'
        metrics['exon_to_full_gene_body_ratio'] = exon_full_gene_body_ratio[metrics['barcode']] if metrics['barcode'] in exon_full_gene_body_ratio else 'NA'
        to_print = [str(metrics[i]) for i in print_metrics]
        print('\t'.join(to_print))

    logging.info('Done')
