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
		fastq_in << [library, readgroup, file(barcode_read)]
		fastqc_in << [library, readgroup, file(insert_read)]
	}
}

Channel.from(fastq_in).into{fastq_in_chan_1; fastq_in_chan_2}

process nucleus_barcodes {
	publishDir "${params.results}/nucleus-barcode-and-umi", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 1

	input:
	set val(library), val(readgroup), file(fastq) from fastq_in_chan_1

	output:
	set val(library), val(readgroup), file("${library}-${readgroup}.index.fastq.gz") into index_chan

	"""
	make-nucleus-barcode-and-umi.py $fastq ${library}-${readgroup}.index.fastq.gz ${library}-${readgroup}.umi.fastq.gz
	"""
}

fastqc_in_chan = Channel.from(fastqc_in)

process fastqc {

	publishDir "${params.results}/fastqc", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2
	
	input:
	set val(library), val(readgroup), file(fastq) from fastqc_in_chan

	output:
	set file(outfile_1), file(outfile_2)

	script:
	outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
        outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

        """
        fastqc $fastq
        """

}

process count_raw_barcodes {
	
	publishDir "${params.results}/alevin-whitelist", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2

	input:
	set val(library), val(readgroup), file(fastq) from index_chan

	output:
	set val(library), file("${library}-${readgroup}.raw_barcode_counts.txt") into count_out_chan

	"""
	count-barcodes.py $fastq > ${library}-${readgroup}.raw_barcode_counts.txt
	"""
}

whitelist_in_chan = count_out_chan.groupTuple()

process alevin_whitelist {

	publishDir "${params.results}/alevin", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2

	input:
	set val(library), file(counts) from whitelist_in_chan

	output:
	set val(library), file("${library}-whitelist.txt") into alevin_chan

	"""
	zcat ${params['barcode-whitelist']} > whitelist.txt
	make-alevin-whitelist.py --min-counts 2 whitelist.txt ${counts.join(' ')} > ${library}-whitelist.txt
	rm whitelist.txt
	"""
}

get_fastq_files = {
	library, read_number ->
	fastq_files = []
	for (readgroup in params.libraries[library]['readgroups'].keySet()) {
		fastq_files << file(params.libraries[library]['readgroups'][readgroup][read_number])
	}
	return fastq_files
}
alevin_in_chan = alevin_chan.map({it -> [it[0], it[1], get_fastq_files(it[0], '1'), get_fastq_files(it[0], '2')]})

process extract_and_correct_barcodes {
	publishDir "${params.results}/alevin", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2
	
	input:
	set val(library), file("${library}-whitelist.txt"), file(first_fastqs), file(second_fastqs) from alevin_in_chan

	output:
	set val(library), file("${library}.1.fastq.gz") into map_chan
	
	"""
	salmon alevin -l ISR -1 ${first_fastqs.join(' ')} -2 ${second_fastqs.join(' ')} -p 5 --chromiumV3 --noQuant --dumpfq --index ${params.salmon_index[get_genome(library)[0]]} --output ${params.results + '/' + library} --tgMap ${params.tgMap[get_genome(library)[0]]} --whitelist ${library}-whitelist.txt | gzip -c > ${library}.1.fastq.gz
	"""	
}


// add genome information here...
tmp = []
for (library in libraries) {
	for (genome in get_genome(library)) {
		tmp << [library, genome]
	}
}

map_in_chan = Channel.from(tmp).combine(map_chan, by: 0)

process map {
	cpus 5
	memory '50 GB'
	publishDir "${params.results}/star/${library}-${genome}", mode: 'rellink', overwrite: true
	validExitStatus 0,141

	input:
	set val(library), val(genome), file("${library}.1.fastq.gz") from map_in_chan

	output:
	set val(library), val(genome), file("Aligned.sortedByCoord.out.bam") into map_out_chan

	"""
	${IONICE} STAR --runThreadN 5 --genomeLoad NoSharedMemory --runRNGseed 789727 --readFilesCommand gunzip -c --genomeDir ${get_star_index(genome)} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --sjdbGTFfile ${get_gtf(genome)} --readFilesIn ${library}.1.fastq.gz
	"""
}


map_out_chan.into{filter_genome_in; count_in_chan}

process count_barcode_per_library {
	cpus 1
	memory '20 GB'
	errorStrategy 'retry'
	maxRetries 10

	publishDir "${params.results}/counts-per-barcode", mode: 'rellink', overwrite: true

	input:
	set val(library), val(genome), file(bam) from count_in_chan

	output:
	set val(library), val(genome), file("${library}-${genome}.counts.txt") into counts_threshold_chan

	"""
	${IONICE} samtools view -F 4 -F 256 -F 2048 -q 255 $bam | cut -f1 | perl -pe 's/.*_(.*)_(.*)\$/\$1/' > barcodes.txt
	count.py barcodes.txt > ${library}-${genome}.counts.txt
	"""
}

counts_threshold_chan.into{counts_threshold_chan_1; counts_threshold_chan_2}
assign_genome_in = counts_threshold_chan_1.groupTuple().transpose().groupTuple() // [library, genome, file] --> [library: [[genome, file], [genome, file]]] -> [library, genome, genome], [library, file, file] -> [library: [[genome, genome], [file, file]]]

// assign a genome to each nucleus. In the case that the library is not chimeric, all nuclei will be assigned to the known correct genome
process assign_genome {

        errorStrategy 'retry'
        maxRetries 10

        publishDir "${params.results}/assign-genome", mode: 'rellink', overwrite: true

        input:
        set val(library), val(genomes), file(x) from assign_genome_in

        output:
        set val(library), file("${library}.genomes.txt") into assign_genome_out

        """
        assign-genome.py --genomes ${genomes.join(' ')} --counts ${library}-*.counts.txt > ${library}.genomes.txt
        """
}

filter_genome_in_chan = filter_genome_in.combine(assign_genome_out, by: 0)
process filter_genome {

        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.results}/filter-genome", mode: 'rellink', overwrite: true

        input:
        set val(library), val(genome), file(bam), file(genome_assignments) from filter_genome_in_chan

        output:
        set val(library), val(genome), file("${library}-${genome}.bam") into filter_genome_out

        """
	zcat ${params['barcode-whitelist']} > whitelist.txt
        grep $genome $genome_assignments | cut -f1 | cat - whitelist.txt | sort | uniq -d > keep.txt
	filter-bam-by-nucleus.py $bam keep.txt ${library}-${genome}.bam
	rm whitelist.txt
        """
}

// add UMI/nucleus  barcode tags
process add_tags {

        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.results}/add-nucleus-and-umi-tags", mode: 'rellink', overwrite: true

        input:
        set val(library), val(genome), file("in.bam") from filter_genome_out

        output:
        set val(library), val(genome), file("${library}-${genome}.bam") into add_tags_out

        """
        copy-barcode-and-umi-from-readname-to-tag.py in.bam ${library}-${genome}.bam
        """
}

add_tags_out.into{feature_counts_in; split_in}

process feature_counts {

	cpus 20 

	input:
	set val(library), val(genome), file("star.bam") from feature_counts_in

	output:
	set val(library), val(genome), file("${library}-${genome}.featureCounts.unsorted.bam") into feature_counts_sort_chan

	"""
	~/sw/subread-1.6.4-Linux-x86_64/bin/featureCounts -a ${get_gtf(genome)} -T 20 -t transcript -o ${library}-${genome} star.bam -R BAM -O -Q 255 --primary -s 1; mv star.bam.featureCounts.bam ${library}-${genome}.featureCounts.unsorted.bam
	"""
}

process sort_feature_counts {

	cpus 10 
	memory '45 GB'
	publishDir "${params.results}/feature-counts", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2

	input:
	set val(library), val(genome), file("${library}-${genome}.featureCounts.unsorted.bam") from feature_counts_sort_chan

	output:
	set val(library), val(genome), file("${library}-${genome}.featureCounts.bam") into sort_feature_counts_out_chan

	"""
	samtools sort -m 4g -o ${library}-${genome}.featureCounts.bam -O BAM -T ${library}.sort -@ 9 ${library}-${genome}.featureCounts.unsorted.bam 
	"""
}

process split {
	publishDir "${params.results}/split-nuclei", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2
	memory { 90.GB * task.attempt }
	cache 'deep'

	input:
	set val(library), val(genome), file(bam) from split_in

	output:
	file("${library}-${genome}-*.sam") optional true into split_out_chan mode flatten

	"""
	split-bam-by-cell-fast.py --prefix ${library}-${genome}- --max-reads-in-memory 100000000 --min-reads-to-output ${params.threshold_to_split} $bam
	"""
}

parse_bam_name = {
        (full, library, genome, barcode) = (it =~ /.*\/(.*)-(.*)-(.*).sam/)[0]
        return [library, genome, barcode, file(it)]
}

split_out_chan.map({it -> parse_bam_name(it)}).into{qorts_in_chan; prune_chan}

flatten_gtf_in = []
for (library in libraries) {
	for (genome in get_genome(library)) {
		flatten_gtf_in << genome
	}
}

flatten_gtf_in_chan = Channel.from(flatten_gtf_in).unique().map({it -> [it, file(get_gtf(it))]})

process flatten_gtf {
	errorStrategy 'retry'
	maxRetries 2
	
input:
	set val(genome), file(gtf) from flatten_gtf_in_chan

	output:
	set val(genome), file("${genome}.flattened.gff") into flatten_gtf_out_chan

	"""
	java -Xmx4g -jar \$QORTS_JAR makeFlatGff --stranded $gtf ${genome}.flattened.gff
	"""

}

qorts_in_chan = qorts_in_chan.map({it -> [it[1], it[0], it[2], it[3]]}).combine(flatten_gtf_out_chan, by: 0)


process qorts {
	publishDir "${params.results}/qorts", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 5

	input:
	set val(genome), val(library), val(barcode), file(nucleus_sam), file(flattened_gff) from qorts_in_chan

	output:
	file("${library}-${genome}-${barcode}") into qorts_out_chan

	"""
	${IONICE} java -Xmx4g -jar \$QORTS_JAR QC --stranded --singleEnded --flatgff $flattened_gff --stranded_fr_secondstrand --runFunctions chromCounts,writeGeneCounts,GeneCalcs,StrandCheck,writeGeneBody --title ${library}-${genome}-${barcode} --chromSizes ${get_chrom_sizes(genome)} ${nucleus_sam} ${get_gtf(genome)} ${library}-${genome}-${barcode}
	"""
}

extract_qorts_in_chan = qorts_out_chan.collate(1000)

process extract_qorts {
	errorStrategy 'retry'
	maxRetries 2

	input:
	file(dirs) from extract_qorts_in_chan

	output:
	file("metrics.*") into extract_qorts_out_chan

	"""
	echo ${dirs.join(' ')} | perl -pe 's/ /\\n/' > collect.txt
	extract-qorts-metrics.py collect.txt metrics.
	"""
}

parse_metric_name = {
        (full, metric_name) = (it =~ /.*\/metrics.(.*).txt/)[0]
        return metric_name
}
merge_extract_qorts_in_chan = extract_qorts_out_chan.flatten().map({it -> [parse_metric_name(it), file(it)]}).groupTuple()

process merge_extract_qorts_chan {

	publishDir "${params.results}/qorts-metrics-extracted", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 2

	input:
	set val(metric_name), file("metrics.*.txt") from merge_extract_qorts_in_chan

	output:
	file("${metric_name}.txt")

	"""
	cat metrics*.txt | awk 'NR==1' > ${metric_name}.txt
	cat metrics*.txt | grep -v library >> ${metric_name}.txt
	"""
}


process prune {
	publishDir "${params.results}/prune", mode: 'rellink', overwrite: true

	input:
	set val(library), val(genome), file(bam) from sort_feature_counts_out_chan

	output:
	set val(library), val(genome), file("${library}-${genome}.pruned.bam"), file("${library}-${genome}.pruned.bam.bai") into dedup_in_chan

	"""
	${IONICE} samtools view -h -b -q 255 -F 4 -F 256 -F 2048 $bam > ${library}-${genome}.pruned.bam
	samtools index ${library}-${genome}.pruned.bam
	"""
}

process dedup {
	publishDir "${params.results}/dedup", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 3
	memory { 30.GB * task.attempt }

	input:
	set val(library), val(genome), file(bam), file(index) from dedup_in_chan

	output:
	set val(library), val(genome), file("${library}-${genome}.dedup.bam"), file("${library}-${genome}.dedup.log") into count_matrix_in_chan

	"""
	umi_tools dedup --random-seed=174983 --per-gene --per-cell --gene-tag XT --assigned-status-tag XS --output-stats {params.metrics} -I $bam -S ${library}-${genome}.dedup.bam -L ${library}-${genome}.dedup.log --extract-umi-method tag --umi-tag UR --cell-tag CB --method directional
	"""
}

process count_matrix {
	publishDir "${params.results}/counts", mode: 'rellink', overwrite: true
	errorStrategy 'retry'
	maxRetries 3
	memory { 30.GB * task.attempt }
	maxForks 3

	input:
	set val(library), val(genome), file(bam), file(log) from count_matrix_in_chan

	output:
	set val(library), val(genome), file("${library}-${genome}.counts.txt")

	"""
	make-count-matrix.py $bam | perl -pe 's/^/${library}\t/' > ${library}-${genome}.counts.txt
	"""
}
