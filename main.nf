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

    publishDir "${params.results}/starsolo/${library}-${genome}"
    memory { 70.GB + (30.GB * task.attempt) }
    cpus 10
    tag "${library}-${genome}"
    container 'library://porchard/default/star:2.7.10a'
    maxRetries 3
    time '48h'

    input:
    tuple val(library), val(genome), path(barcode_fastq), path(insert_fastq)

    output:
    tuple val(library), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), path("${library}-${genome}.Log.final.out"), path("${library}-${genome}.Log.out"), path("${library}-${genome}.Log.progress.out"), path("${library}-${genome}.SJ.out.tab"), path("${library}-${genome}.Solo.out")
    tuple val(library), val(genome), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), path("${library}-${genome}.Solo.out"), emit: for_qc
    tuple val(library), val(genome), path("${library}-${genome}.Solo.out"), emit: solo_out
    tuple val(library), val(genome), path("${library}-${genome}.Aligned.sortedByCoord.out.bam"), emit: for_prune
    path("${library}-${genome}.Log.final.out"), emit: for_multiqc

    script:
    soloUMIlen = ['V2', 'GEX5'].contains(params.chemistry) ? 10 : 12
    clip5pNbases = params.chemistry == 'GEX5' ? '39 0' : '0'
    soloBarcodeMate = params.chemistry == 'GEX5' ? 1 : 0
    fastq_1 = params.chemistry == 'GEX5' ? barcode_fastq : insert_fastq
    fastq_2 = params.chemistry == 'GEX5' ? insert_fastq : barcode_fastq

    """
    ${IONICE} STAR --soloBarcodeReadLength 0 --runThreadN 10 --outFileNamePrefix ${library}-${genome}. --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM ${30000000000 + (10000000000 * task.attempt)} --outBAMsortingBinsN 200 --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --soloType CB_UMI_Simple --soloUMIlen $soloUMIlen --soloFeatures GeneFull_ExonOverIntron GeneFull GeneFull_Ex50pAS Gene SJ --soloMultiMappers Uniform PropUnique EM Rescue --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCellFilter None --soloCBwhitelist ${params['barcode-whitelist']} --clip5pNbases ${clip5pNbases} --soloBarcodeMate ${soloBarcodeMate} --readFilesIn ${fastq_1.join(',')} ${fastq_2.join(',')}
    find ${library}-${genome}.Solo.out -type d -exec chmod 755 {} +
    find ${library}-${genome}.Solo.out -type f -exec chmod 644 {} +
    """

}


process star_multiqc {

    publishDir "${params.results}/multiqc/star"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '2h'

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

    publishDir "${params.results}/prune"
    tag "${library}-${genome}"
    container 'library://porchard/default/general:20220107'
    memory '2 GB'
    cpus 1
    time '5h'

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple path("${library}-${genome}.before-dedup.bam"), path("${library}-${genome}.before-dedup.bam.bai")

    script:
    flags = params.chemistry == 'GEX5' ? '-f 3 -F 256 -F 2048' : '-F 4 -F 256 -F 2048'

    """
    ${IONICE} samtools view -h -b -q 255 ${flags} $bam > ${library}-${genome}.before-dedup.bam && samtools index ${library}-${genome}.before-dedup.bam
    """

}


process fastqc {

    publishDir "${params.results}/fastqc"
    tag "${library} ${readgroup}"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '5h'

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

    publishDir "${params.results}/multiqc/fastq"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '5h'

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
    cpus 1
    time '5h'

    input:
    tuple val(library), val(genome), path("star.bam"), path(solo_out)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.qc.txt")


    """
    qc-from-starsolo.py star.bam ${solo_out}/GeneFull_ExonOverIntron/raw/matrix.mtx ${solo_out}/GeneFull_ExonOverIntron/raw/barcodes.tsv > ${library}-${genome}.qc.txt
    """

}


process plot_qc {

    memory '15 GB'
    publishDir "${params.results}/qc"
    tag "${library}-${genome}"
    container 'library://porchard/default/dropkick:20220225'
    cpus 1
    time '5h'

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.metrics.png"), path("${library}-${genome}.suggested-thresholds.tsv")


    """
    plot-qc-metrics.py --prefix ${library}-${genome}. $metrics
    """

}


process interactive_barcode_rank_plot {

    memory '15 GB'
    publishDir "${params.results}/interactive-barcode-rank-plots"
    tag "${library}-${genome}"
    container "docker://porchard/plotly:20230705"
    cpus 1
    time '3h'

    input:
    tuple val(library), val(genome), path(solo_out)

    output:
    path("${library}-${genome}.barcode-rank-plot.html")


    """
    interactive-barcode-rank-plot.py ${solo_out}/GeneFull_ExonOverIntron/raw/matrix.mtx ${library}-${genome}.barcode-rank-plot.html
    """

}


process cellbender {

    cpus 1
    memory '40 GB'
    publishDir "${params.results}/cellbender"
    container 'docker://porchard/cellbender:0.3.0'
    time '72h'
    errorStrategy 'ignore'

    input:
    tuple val(library), val(genome), path(solo_out)

    output:
    path("${library}*")
    path("${library}*.h5"), emit: h5_files

    """
    cp ${solo_out}/GeneFull_ExonOverIntron/raw/matrix.mtx matrix.mtx
    cp ${solo_out}/GeneFull_ExonOverIntron/raw/features.tsv genes.tsv
    cp ${solo_out}/GeneFull_ExonOverIntron/raw/barcodes.tsv barcodes.tsv

    cellbender remove-background --cuda --epochs 150 --fpr 0.01 0.05 0.1 --input . --output ./${library}-${genome}.cellbender.h5
    cp .command.log ${library}-${genome}.log
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
    qc(starsolo_out.for_qc) | plot_qc
    interactive_barcode_rank_plot(starsolo_out.solo_out)
    cellbender(starsolo_out.solo_out)

}
