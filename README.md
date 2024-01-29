# NextFlow pipeline for 10X snRNA-seq data

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library).


## Configuration
Paths to reference files must be included in the `nextflow.config` file -- check that file and change paths accordingly. These include:

1. STAR indices (compatible with STAR v. 2.7.9a)
2. GTF files

When launching the pipeline, as shown in the `nextflow` command below, you'll also need to set the following:

1. The location of the results directory (e.g., `--results /path/to/results`)
2. The 10X Chromium chemistry version ('V2', 'V3', 'GEX5', or 'multiome'; e.g., `--chemistry multiome``). Note for 5' GEX ('GEX5'), read 1 is currently assumed to include cDNA in addition to the CB, UMI, and adapter; i.e. it's length should be > 39 bp.
3. The location of the barcode whitelist (e.g., `--barcode-whitelist /path/to/737K-arc-v1.txt`). The whitelists are within this directory but must be unzipped before use. For Chromium 3' GEX v3, this is the 3M-february-2018.txt.gz file; for 3' GEX v2 or 5' GEX v1 and v2, this is the 737K-august-2016.txt.gz file; for multiome, this is 737K-arc-v1.txt.gz)

Lastly, you'll need to include information about each RNA-seq library, including the genome to which it should be mapped, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json. For each readgroup, the '1' fastq file corresponds to the sequencing read including the UMI and the nucleus index (and, for 5' GEX, some cDNA); the '2' fastq file refers to the sequencing read representing only cDNA. Also, note that the 'genome' attribute is given as a list (because I will be adding the ability to map to multiple genomes, in the case that nuclei from multiple species are mixed together).

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -resume -params-file library-config.json --barcode-whitelist /path/to/barcode-whitelist.txt --chemistry multiome --results /path/to/results /path/to/main.nf
```

## Output
* `cellbender/*`: Cellbender results
* `dropkick/*`: Results of running dropkick per-library (QC plots, as well as per-barcode dropkick score that can be used for selecting quality barcodes for downstream analysis)
* `multiqc/fastq/*`: multiqc summaries of fastqc results
* `multiqc/star/*`: multiqc summaries of STAR logs
* `prune/*`: filtered bam files (duplicates NOT removed)
* `qc/*`: Per-barcode QC metrics and QC metric plots
* `starsolo/*`: starsolo output. Count matrices derived using a variety of counting methods (see STAR manual) are in `starsolo/{library}/{library}.Solo.out/*`