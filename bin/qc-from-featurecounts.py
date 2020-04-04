#!/usr/bin/env python

import os
import sys
import pysam
import json
import argparse
import logging
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser(description='Gather stats from a snRNA-seq bam file.', add_help = True)
parser.add_argument('bam', type = str,  help = 'BAM file, sorted by read name.')
parser.add_argument('--cell-tag', dest='cell_tag', type = str, default = 'CB', help = 'Tag denoting the cell/nucleus (default: CB)')
parser.add_argument('--gene-tag', dest='gene_tag', type = str, default = 'GX', help = 'Tag denoting the gene assignment (default: GX)')
parser.add_argument('--umi-tag', dest='umi_tag', type = str, default = 'UB', help = 'Tag denoting the UMI (default: UB)')
parser.add_argument('--min-reads', dest='min_reads', type = int, default = 0, help = 'Suppress output for cells with fewer than this many reads (default: 0).')
args = parser.parse_args()

CELL_TAG = args.cell_tag
GENE_TAG = args.gene_tag
UMI_TAG = args.umi_tag
MIN_READS_FOR_OUTPUT = args.min_reads
GENE_STATUS_TAG = None


class Cell:

    def __init__(self, barcode):
        self.barcode = barcode
        self.umis = dict() # UMI --> gene --> count
        self.total_reads = 0
        self.reads_no_gene = 0
        self.reads_no_umi = 0
        self.uniquely_mapped_reads = 0
        self.reads_passing_filters = 0
        self.chromosome_read_counts = dict() # chrom -> count


    def record_alignment(self, read):
        self.total_reads += 1
        if not read.has_tag(GENE_TAG):
            self.reads_no_gene += 1
        if not read.has_tag(UMI_TAG):
            self.reads_no_umi += 1
        if read.mapping_quality == 255:
            self.uniquely_mapped_reads += 1
            chrom = read.reference_name
            if chrom not in self.chromosome_read_counts:
                self.chromosome_read_counts[chrom] = 0
            self.chromosome_read_counts[chrom] += 1
        if read.has_tag(UMI_TAG) and read.has_tag(GENE_TAG) and read.mapping_quality == 255 and not read.is_secondary:
            self.reads_passing_filters += 1
            umi = read.get_tag(UMI_TAG)
            gene = read.get_tag(GENE_TAG)
            if umi not in self.umis:
                self.umis[umi] = dict()
            if gene not in self.umis[umi]:
                self.umis[umi][gene] = 0
            self.umis[umi][gene] += 1


    def number_genes(self):
        all_genes = set()
        for umi, genes in self.umis.items():
            for gene in genes:
                all_genes.add(gene)
        return len(all_genes)


    def number_umis(self):
        return len(self.umis)


    def handle_multigeneumi(self):
        # each UMI is only allowed to be assigned to one gene. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview
        # therefore, find and fix such cases as described in the link
        new_umis = dict()
        for umi, genes in self.umis.items():
            if len(genes) > 1:
                # find the gene with the max counts...
                max_value = max(genes.values())
                genes_with_max_value = list()
                for gene, count in genes.items():
                    if count == max_value:
                        genes_with_max_value.append(gene)
                if len(genes_with_max_value) > 1:
                    new_umis[umi] = dict()
                else:
                    new_umis[umi] = {genes_with_max_value[0]: max_value}
            else:
                new_umis[umi] = genes
        self.umis = new_umis



    def gene_counts(self):
        counts = dict()  # gene --> umi counts
        # self.handle_multigeneumi()
        for umi, genes in self.umis.items():
            for gene in genes:
                if gene not in counts:
                    counts[gene] = 0
                counts[gene] += 1
        return counts


    def cumulative_gene_dist(self):
        x = self.gene_counts()
        s = sum(x.values())
        number_genes = 0
        cumulative = dict()
        for gene, count in sorted(x.items(), key=lambda y: y[1], reverse=True):
            number_genes += 1
            cumulative[number_genes] = float(count)/s if number_genes == 1 else cumulative[number_genes-1] + (float(count)/s)
        return cumulative


    def estimate_fraction_duplicate(self):
        reads_assigned_to_gene = self.total_reads - self.reads_no_gene
        unique_reads_assigned_to_gene = sum(self.gene_counts().values())
        return (reads_assigned_to_gene - unique_reads_assigned_to_gene) / reads_assigned_to_gene if reads_assigned_to_gene > 0 else np.nan


    def gather_metrics(self):
        metrics = dict()
        metrics['barcode'] = self.barcode
        metrics['total_reads'] = self.total_reads
        metrics['uniquely_mapped_reads'] = self.uniquely_mapped_reads
        metrics['reads_passing_filters'] = self.reads_no_umi
        metrics['chromosome_read_counts'] = self.chromosome_read_counts
        metrics['reads_no_umi'] = self.reads_no_umi
        metrics['reads_no_gene'] = self.reads_no_gene
        metrics['number_genes'] = self.number_genes()
        metrics['number_umis'] = self.number_umis()
        metrics['estimated_fraction_duplicate'] = self.estimate_fraction_duplicate()
        metrics['gene_counts'] = {gene: count for gene, count in self.gene_counts().items()}
        metrics['cumulative_gene_dist'] = {gene_number: fraction for gene_number, fraction in sorted(self.cumulative_gene_dist().items())}
        return metrics



if __name__ == '__main__':
    cells = dict()
    no_cell_tag = 0
    total_reads = 0

    with pysam.AlignmentFile(args.bam, 'rb') as f:
        for read in f.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {} reads'.format(total_reads))
            if not read.has_tag(CELL_TAG):
                no_cell_tag += 1
                continue
            barcode = read.get_tag(CELL_TAG)
            if barcode not in cells:
                cells[barcode] = Cell(barcode)
            cells[barcode].record_alignment(read)

    logging.info('Finished reading bam file.')
    JSON = []
    for barcode, cell in cells.items():
        if cell.total_reads >= MIN_READS_FOR_OUTPUT:
            logging.info('Preparing output for cell with barcode {}'.format(barcode))
            JSON.append(cell.gather_metrics())
    print(json.dumps(JSON, indent=4, sort_keys=True))
    logging.info('Encountered {} reads (of {}) with no cell tag'.format(no_cell_tag, total_reads))
