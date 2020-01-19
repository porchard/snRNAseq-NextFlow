# NextFlow pipeline for 10X snATAC-seq data

## Dependencies
If you have Singularity installed, you can use the config provided here ('Singularity') to build a container with all the dependencies.

Otherwise, you'll need to have the following installed:
1. bedtools
2. biopython
3. STAR
4. fastqc
5. samtools
6. pysam
7. pybedtools

I've used this pipeline with NextFlow v. 19.04.1

## Configuration
Paths to various generic files (e.g., STAR indices) must be included in the nextflow.config file -- check that file and change paths accordingly. These include:

1. STAR indices
2. Barcode whitelist (for Chromium v3, that is the 3M-february-2018.txt file; for v2, that is the 737K-august-2016.txt file)

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results') as well as the 10X Chromium chemistry version ('V2' or 'V3')

Lastly, you'll need to include information about each RNA-seq library, including the genome to which it should be mapped, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json. For each readgroup, the '1' fastq file corresponds to the sequencing read including the UMI and the nucleus index; the '2' fastq file refers to the sequencing read representing the actual transcript. Also, note that the 'genome' attribute is given as a list (because I will be adding the ability to map to multiple genomes, in the case that nuclei from multiple species are mixed together).

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -with-singularity /path/to/Singularity.simg -params-file library-config.json --results /path/to/results /path/to/main.nf
```
