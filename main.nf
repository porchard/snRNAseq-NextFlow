#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IONICE = 'ionice -c2 -n7'

def get_star_index (genome) {
	return(params.star_index[genome])
}

def get_gtf (genome) {
	return(params.gtf[genome])
}

def get_chrom_sizes (genome) {
	return(get_star_index(genome) + '/chrNameLength.txt')
}
	
def get_genome (library) {
	return(params.libraries[library].genome)
}

def library_to_readgroups (library) {
	return(params.libraries[library].readgroups.keySet())
}

def library_and_readgroup_to_fastqs (library, readgroup) {
	return(params.libraries[library].readgroups[readgroup])
}


process starsolo {

    publishDir "${params.results}/starsolo/${library}-${genome}", mode: 'rellink', overwrite: true
    memory { 70.GB + (30.GB * task.attempt) }
    cpus 10
    tag "${library}-${genome}"
    container 'library://porchard/default/star:2.7.10a'
    maxRetries 3

    input:
    tuple val(library), val(genome), path(barcode_fastq), path(insert_fastq)

    output:
    tuple val(library), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), path("${library}-${genome}.Log.final.out"), path("${library}-${genome}.Log.out"), path("${library}-${genome}.Log.progress.out"), path("${library}-${genome}.SJ.out.tab"), path("${library}-${genome}.Solo.out")
    tuple val(library), val(genome), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), path("${library}-${genome}.Solo.out"), emit: for_qc
    tuple val(library), val(genome), path("${library}-${genome}.Solo.out"), emit: for_dropkick
    tuple val(library), val(genome), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), emit: for_prune
    path("${library}-${genome}.Log.final.out"), emit: for_multiqc

    script:
    soloUMIlen = params.chemistry == 'V2' ? 10 : 12

    """
    ${IONICE} STAR --soloBarcodeReadLength 0 --runThreadN 10 --outFileNamePrefix ${library}-${genome}. --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM ${30000000000 + (10000000000 * task.attempt)} --outBAMsortingBinsN 200 --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --soloType Droplet --soloUMIlen $soloUMIlen --soloFeatures Transcript3p Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS SJ Velocyto --soloMultiMappers Uniform PropUnique EM Rescue --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCellFilter None --soloCBwhitelist ${params['barcode-whitelist']} --readFilesIn ${insert_fastq.join(',')} ${barcode_fastq.join(',')}
    """

}


process star_multiqc {

    publishDir "${params.results}/multiqc/star", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    memory '4 GB'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process prune {

    publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
    maxForks 10
    tag "${library}-${genome}"
    container 'library://porchard/default/general:20220107'
    memory '2 GB'

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
    container 'library://porchard/default/general:20220107'
    memory '4 GB'

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


process fastq_multiqc {

    publishDir "${params.results}/multiqc/fastq", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    memory '4 GB'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process qc {

    memory '25 GB'
    publishDir "${params.results}/qc"
    tag "${library}-${genome}"
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(library), val(genome), path("star.bam"), path(solo_out)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.qc.txt")


    """
    qc-from-starsolo.py star.bam ${solo_out}/GeneFull/raw/matrix.mtx ${solo_out}/GeneFull/raw/barcodes.tsv > ${library}-${genome}.qc.txt
    """

}


process dropkick {

    memory '150 GB'
    cpus 5
    publishDir "${params.results}/dropkick"
    tag "${library}-${genome}"
    container 'library://porchard/default/dropkick:20220225'
    errorStrategy 'ignore'

    input:
    tuple val(library), val(genome), path(solo_out)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.dropkick-score.tsv"), emit: dk_score
    path("*.png")


    """
    mkdir -p dropkick-in
    cp ${solo_out}/GeneFull/raw/matrix.mtx dropkick-in/matrix.mtx
    cp ${solo_out}/GeneFull/raw/features.tsv dropkick-in/genes.tsv
    cp ${solo_out}/GeneFull/raw/barcodes.tsv dropkick-in/barcodes.tsv
    run-dropkick.py dropkick-in/ ${library}-${genome}.
    """

}


process plot_qc {

    memory '15 GB'
    publishDir "${params.results}/qc"
    tag "${library}-${genome}"
    container 'library://porchard/default/dropkick:20220225'

    input:
    tuple val(library), val(genome), path(metrics), path(dk)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.metrics.png"), path("${library}-${genome}.suggested-thresholds.tsv")


    """
    plot-qc-metrics.py --prefix ${library}-${genome}. $metrics $dk
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

    fastqc(Channel.from(fastqc_in)).flatten().toSortedList() | fastq_multiqc
    star_in = Channel.from(fastq_in).groupTuple(by: [0,1])
    starsolo_out = starsolo(star_in)
    star_multiqc(starsolo_out.for_multiqc.toSortedList())
    prune(starsolo_out.for_prune)
    qc(starsolo_out.for_qc).combine(dropkick(starsolo_out.for_dropkick).dk_score, by: [0, 1]) | plot_qc

}
