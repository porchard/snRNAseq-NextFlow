#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

libraries = params.libraries.keySet()

get_star_index = {
	genome ->
	params.star_index[genome]
}

get_gtf = {
	genome ->
	return(get_star_index(genome) + '/annotation.gtf')
}

get_chrom_sizes = {
	genome ->
	return(get_star_index(genome) + '/chrNameLength.txt')
}
	
get_genome = {
	library ->
	params.libraries[library].genome
}

library_to_readgroups = {
	library ->
	params.libraries[library].readgroups.keySet()
}


library_and_readgroup_to_fastqs = {
	library, readgroup ->
	params.libraries[library].readgroups[readgroup]
}


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

star_in = Channel.from(fastq_in).groupTuple(by: [0,1 ])

process starsolo {

	publishDir "${params.results}/starsolo/${library}-${genome}", mode: 'rellink', overwrite: true
	memory '75 GB'
	cpus 10

	input:
	set val(library), val(genome), file(barcode_fastq), file(insert_fastq) from star_in	

	output:
	set val(library), file("Aligned.sortedByCoord.out.bam"), file("Log.final.out"), file("Log.out"), file("Log.progress.out"), file("SJ.out.tab"), file("Solo.out")
	set val(library), val(genome), file("Aligned.sortedByCoord.out.bam") into feature_counts_in
	set val(library), val(genome), file("Aligned.sortedByCoord.out.bam") into prune_in

	script:
	soloUMIlen = params.chemistry == 'V2' ? 10 : 12

	"""
	${IONICE} STAR --soloBarcodeReadLength 0 --runThreadN 10 --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c --outSAMattributes NH HI nM AS CR CY CB UR UY UB sM GX GN --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --soloType Droplet --soloUMIlen $soloUMIlen --soloFeatures Transcript3p Gene GeneFull SJ Velocyto --soloUMIfiltering MultiGeneUMI --soloCBmatchWLtype 1MM_multi_pseudocounts --soloCellFilter None --soloCBwhitelist ${params['barcode-whitelist']} --readFilesIn ${insert_fastq.join(',')} ${barcode_fastq.join(',')}
	"""

}

process prune {

	publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
	maxForks 10

	input:
	set val(library), val(genome), file(bam) from prune_in

	output:
	set file("${library}-${genome}.before-dedup.bam"), file("${library}-${genome}.before-dedup.bam.bai")

	"""
	${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam > ${library}-${genome}.before-dedup.bam
	samtools index ${library}-${genome}.before-dedup.bam
	"""

}

process fastqc {

	publishDir "${params.results}/fastqc", mode: 'rellink', overwrite: true
	maxForks 6
	
	input:
	set val(library), val(readgroup), file(fastq) from Channel.from(fastqc_in)

	output:
	set file(outfile_1), file(outfile_2)

	script:
	outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
    	outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

    	"""
    	fastqc $fastq
    	"""

}


process feature_counts {

	cpus 20
	publishDir "${params.results}/feature-counts"

	input:
	set val(library), val(genome), file("star.bam") from feature_counts_in

	output:
	set val(library), val(genome), file("${library}-${genome}.featureCounts.bam") into feature_counts_out
	

	"""
	~/sw/subread-1.6.4-Linux-x86_64/bin/featureCounts -a ${get_gtf(genome)} -T 20 -t transcript -o ${library}-${genome} star.bam -R BAM -O -s 1; mv star.bam.featureCounts.bam ${library}-${genome}.featureCounts.bam
	"""

}


process qc {

	publishDir "${params.results}/qc"
	memory '75 GB'

	input:
	set val(library), val(genome), file(bam) from feature_counts_out

	output:
	set val(library), val(genome), file("${library}-${genome}.qc.json.gz")

	"""
	qc-from-featurecounts.py --cell-tag CB --gene-tag XT --umi-tag UB --min-reads 100 $bam | gzip -c > ${library}-${genome}.qc.json.gz
	"""

}
