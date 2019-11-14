#!/usr/bin/env python
# coding: utf-8

import sys
import os
import logging
import gzip
import csv
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')


class QoRTs():

    def __init__(self, directory = None, d = None, library = None):
        self.library = library
        self.read_length = None
        self.chrom_count = dict() # chromosome --> count
        self.gene_count = dict() # gene --> count
        self.insert_sizes = dict() # insert_size --> count
        self.quality_read_pairs = None
        self.total_read_pairs = None
        self.dropped_not_proper_pair = None
        self.dropped_fail_vendor_qc = None
        self.dropped_marked_not_valid = None
        self.dropped_chrom_mismatch = None
        self.dropped_pairs_strands_mismatch = None
        self.dropped_ignored_chromosome = None
        self.dropped_not_unique_alignment = None
        self.kept_not_unique_alignment = None
        self.min_read_length = None
        self.max_read_length = None
        self.max_legal_phred_score = None
        self.is_single_end = None
        self.mapping_location = dict()
        self.number_genes_with_zero_counts = None
        self.number_genes_with_nonzero_counts = None
        self.splice_loci = dict()
        self.splice_events = dict()
        self.strand_test = dict()
        self.chromosomes_covered = None
        self.read_length = None
        self.warning = None
        self.error = None
        self.fraction_strandedness = None
        self.gene_body_coverage_pct = dict() # QUANTILE -> total/bottom_half/upper_mid_quartile/75_to_90/high
        
        if directory is not None:
            with gzip.open(os.path.join(directory, 'QC.chromCount.txt.gz'), 'rt') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for line in reader:
                    chrom, count = line['CHROM'], int(line['CT'])
                    self.chrom_count[chrom] = count

            with gzip.open(os.path.join(directory, 'QC.geneCounts.txt.gz'), 'rt') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for line in reader:
                    gene, count = line['GENEID'], int(line['COUNT'])
                    self.gene_count[gene] = count

            if os.path.isfile(os.path.join(directory, 'QC.insert.size.txt.gz')):
                with gzip.open(os.path.join(directory, 'QC.insert.size.txt.gz'), 'rt') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for line in reader:
                        insert_size, count = int(line['InsertSize']), int(line['Ct'])
                        self.insert_sizes[insert_size] = count

            if os.path.isfile(os.path.join(directory, 'QC.geneBodyCoverage.byExpr.avgPct.txt.gz')):
                with gzip.open(os.path.join(directory, 'QC.geneBodyCoverage.byExpr.avgPct.txt.gz'), 'rt') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for line in reader:
                        self.gene_body_coverage_pct[line['QUANTILE']] = {
                                'total': float(line['TOTAL']),
                                'bottom_half': float(line['1.bottomHalf']),
                                'upper_mid_quartile': float(line['2.upperMidQuartile']),
                                '75_to_90': float(line['3.75to90']),
                                'high': float(line['4.high'])
                            }

            with open(os.path.join(directory, 'QC.summary.txt'), 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for line in reader:
                    field = line['FIELD']
                    count = line['COUNT']
                    if field == 'READ_PAIR_OK':
                        self.quality_read_pairs = int(count)
                    elif field == 'TOTAL_READ_PAIRS':
                        self.total_read_pairs = int(count)
                    elif field == 'DROPPED_NOT_PROPER_PAIR':
                        self.dropped_not_proper_pair = int(count)
                    elif field == 'DROPPED_READ_FAILS_VENDOR_QC':
                        self.dropped_fail_vendor_qc = int(count)
                    elif field == 'DROPPED_MARKED_NOT_VALID':
                        self.dropped_marked_not_valid = int(count)
                    elif field == 'DROPPED_CHROMS_MISMATCH':
                        self.dropped_chrom_mismatch = int(count)
                    elif field == 'DROPPED_PAIR_STRANDS_MISMATCH':
                        self.dropped_pairs_strands_mismatch = int(count)
                    elif field == 'DROPPED_IGNORED_CHROMOSOME':
                        self.dropped_ignored_chromosome = int(count)
                    elif field == 'DROPPED_NOT_UNIQUE_ALIGNMENT':
                        self.dropped_not_unique_alignment = int(count)
                    elif field == 'KEPT_NOT_UNIQUE_ALIGNMENT':
                        self.kept_not_unique_alignment = int(count)
                    elif field == 'minObservedReadLength':
                        self.min_read_length = int(count)
                    elif field == 'maxObservedReadLength':
                        self.max_read_length = int(count)
                    elif field == 'maxLegalPhredScore':
                        self.max_legal_phred_score = int(count)
                    elif field == 'IS_SINGLE_END':
                        self.is_single_end = bool(int(count))
                    elif 'ReadPairs_' in field:
                        m = re.match('ReadPairs_(.*)', field)
                        category = m.group(1)
                        self.mapping_location[category] = int(count)
                    elif field == 'Genes_WithZeroCounts':
                        self.number_genes_with_zero_counts = int(count)
                    elif field == 'Genes_WithNonzeroCounts':
                        self.number_genes_with_nonzero_counts = int(count)
                    elif 'SpliceLoci_' in field:
                        m = re.match('SpliceLoci_(.*)', field)
                        category = m.group(1)
                        self.splice_loci[category] = int(count)
                    elif 'SpliceEvents_' in field:
                        m = re.match('SpliceEvents_(.*)', field)
                        category = m.group(1)
                        self.splice_events[category] = int(count)
                    elif 'StrandTest_' in field:
                        if field == 'StrandTest_STRANDEDNESS_MATCHES_INFERRED':
                            continue
                        m = re.match('StrandTest_(.*)', field)
                        category = m.group(1)
                        self.strand_test[category] = int(count)
                    elif field == 'NumberOfChromosomesCovered':
                        self.chromosomes_covered = int(count)
                    elif field == 'READ_LENGTH':
                        self.read_length = int(count)
                    elif field == 'COMPLETED_WITHOUT_WARNING':
                        self.warning = not bool(int(count))
                    elif field == 'COMPLETED_WITHOUT_ERROR':
                        self.error = not bool(int(count))

            if self.strand_test['frFirstStrand'] == 0 and self.strand_test['frSecondStrand'] == 0:
                self.fraction_strandedness = 0
            else:
                self.fraction_strandedness = float(max([self.strand_test['frFirstStrand'], self.strand_test['frSecondStrand']]))/ (self.strand_test['frFirstStrand'] + self.strand_test['frSecondStrand'])
        elif d is not None:
            for key, val in d.items():
                setattr(self, key, val)
        else:
            pass
    
    def as_dict(self):
        attributes = [x for x in dir(self) if not x.startswith('__') and not callable(getattr(self, x))]
        a = dict()
        for key in attributes:
            val = getattr(self, key)
            a[key] = val

        return a
