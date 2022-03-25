# NextFlow pipeline for 10X snATAC-seq data

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library).


## Configuration
Paths to reference files must be included in the nextflow.config file -- check that file and change paths accordingly. These include:

1. STAR indices (compatible with STAR v. 2.7.9a)
2. GTF files
3. Barcode whitelist (for Chromium v3, that is the 3M-february-2018.txt file; for v2, that is the 737K-august-2016.txt file; for multiome, that is 737K-arc-v1.txt)

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results') as well as the 10X Chromium chemistry version ('V2', 'V3', or 'multiome')

Lastly, you'll need to include information about each RNA-seq library, including the genome to which it should be mapped, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json. For each readgroup, the '1' fastq file corresponds to the sequencing read including the UMI and the nucleus index; the '2' fastq file refers to the sequencing read representing the actual transcript. Also, note that the 'genome' attribute is given as a list (because I will be adding the ability to map to multiple genomes, in the case that nuclei from multiple species are mixed together).

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -params-file library-config.json --chemistry multiome --results /path/to/results /path/to/main.nf
```

## Output
* `dropkick/*`: Results of running dropkick per-library (QC plots, as well as per-barcode dropkick score that can be used for selecting quality barcodes for downstream analysis)
* `multiqc/fastq/*`: multiqc summaries of fastqc results
* `multiqc/star/*`: multiqc summaries of STAR logs
* `prune/*`: filtered bam files (duplicates NOT removed)
* `qc/*`: Per-barcode QC metrics and QC metric plots
* `starsolo/*`: starsolo output. Count matrices derived using a variety of counting methods (see STAR manual) are in `starsolo/{library}/{library}.Solo.out/*`