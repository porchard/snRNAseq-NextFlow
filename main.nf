#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IONICE = 'ionice -c2 -n7'

def get_star_index (genome) {
	params.star_index[genome]
}

def get_gtf (genome) {
	return(get_star_index(genome) + '/annotation.gtf')
}

def get_chrom_sizes (genome) {
	return(get_star_index(genome) + '/chrNameLength.txt')
}
	
def get_genome (library) {
	params.libraries[library].genome
}

def library_to_readgroups (library) {
	params.libraries[library].readgroups.keySet()
}

def library_and_readgroup_to_fastqs (library, readgroup) {
	params.libraries[library].readgroups[readgroup]
}


process starsolo {

    publishDir "${params.results}/starsolo/${library}-${genome}", mode: 'rellink', overwrite: true
    memory '75 GB'
    cpus 10
    tag "${library}-${genome}"

    input:
    tuple val(library), val(genome), path(barcode_fastq), path(insert_fastq)

    output:
    tuple val(library), path("Aligned.sortedByCoord.out.bam"), path("Log.final.out"), path("Log.out"), path("Log.progress.out"), path("SJ.out.tab"), path("Solo.out")
    tuple val(library), val(genome), path("Aligned.sortedByCoord.out.bam"), path("Solo.out"), emit: for_qc
    tuple val(library), val(genome), path("Aligned.sortedByCoord.out.bam"), emit: for_prune

    script:
    soloUMIlen = params.chemistry == 'V2' ? 10 : 12

    """
    ${IONICE} STAR --soloBarcodeReadLength 0 --runThreadN 10 --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --soloType Droplet --soloUMIlen $soloUMIlen --soloFeatures Transcript3p Gene GeneFull SJ Velocyto --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCellFilter None --soloCBwhitelist ${params['barcode-whitelist']} --readFilesIn ${insert_fastq.join(',')} ${barcode_fastq.join(',')}
    """

}

process prune {

    publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
    maxForks 10
    tag "${library}-${genome}"

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple path("${library}-${genome}.before-dedup.bam"), path("${library}-${genome}.before-dedup.bam.bai")

    """
    ${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam > ${library}-${genome}.before-dedup.bam && samtools index ${library}-${genome}.before-dedup.bam
    """

}

process fastqc {

    publishDir "${params.results}/fastqc", mode: 'rellink', overwrite: true
    maxForks 6
    tag "${library} ${readgroup}"

    input:
    tuple val(library), val(readgroup), path(fastq)

    output:
    tuple path(outfile_1), path(outfile_2)

    script:
    outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
    outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

    """
    fastqc $fastq
    """

}


process qc {

    memory '25 GB'
    publishDir "${params.results}/qc"
    tag "${library}-${genome}"

    input:
    tuple val(library), val(genome), path("star.bam"), path(solo_out)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.qc.txt")


    """
    qc-from-starsolo.py star.bam ${solo_out}/GeneFull/raw/matrix.mtx ${solo_out}/GeneFull/raw/barcodes.tsv > ${library}-${genome}.qc.txt
    """

}

process plot_qc {

    memory '15 GB'
    publishDir "${params.results}/qc"
    tag "${library}-${genome}"

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.metrics.png")


    """
    plot-qc-metrics.py --prefix ${library}-${genome}. $metrics
    """

}


workflow {

    libraries = params.libraries.keySet()

    fastq_in = []
    fastqc_in = []

    for (library in libraries) {
        for (readgroup in library_to_readgroups(library)) {
            fastqs = library_and_readgroup_to_fastqs(library, readgroup)
            insert_read = fastqs['2']
            barcode_read = fastqs['1']
            fastqc_in << [library, readgroup, file(insert_read)]
            fastqc_in << [library, readgroup, file(barcode_read)]
            for (genome in get_genome(library)) {
                fastq_in << [library, genome, file(barcode_read), file(insert_read)]
            }
        }
    }

    fastqc(Channel.from(fastqc_in))
    star_in = Channel.from(fastq_in).groupTuple(by: [0,1])
    starsolo_out = starsolo(star_in)
    prune(starsolo_out.for_prune)
    qc(starsolo_out.for_qc) | plot_qc

}