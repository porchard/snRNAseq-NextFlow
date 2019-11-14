#!/usr/bin/env python

import argparse
import sys
import os
import json

# Given a library json file, a json file containing other
# config information, and a path to the desired results
# directory, create a config file that can be used
# for the Snakemake ATAC-seq pipeline


def parse_arguments():
    parser = argparse.ArgumentParser(prog='python make_atacseq_config.py',
    description = """
    Given a library json file, a json file containing other
    config information, and a path to the desired results
    directory, create a config file that can be used
    for the Snakemake ATAC-seq pipeline
    """
    )

    parser.add_argument('reference_json', action='store', help='Path to the JSON file encoding e.g. BWA index path (see README).')
    parser.add_argument('library_json', action='store', nargs = '+', help='Path to the JSON file encoding library information.')
    parser.add_argument('-r', '--results', action='store', dest='results', help='Base directory in which the ATAC-seq analysis results should be stored (default: current working directory).')

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()

    atacseq_config = {}
    atacseq_config['libraries'] = {}

    # now we can add in the information containing paths to generic data such as BWA indices, etc.
    with open(args.reference_json, 'r') as r:
        generic_config = json.load(r)
        for k, v in generic_config.items():
            atacseq_config[k] = v

    # read in the library data
    for library_json in args.library_json:
        with open(library_json, 'r') as s:
            atacseq_config['libraries'].update(json.load(s))

    atacseq_config['results'] = args.results or os.getcwd()

    print(json.dumps(atacseq_config, sort_keys = True, indent=4))
