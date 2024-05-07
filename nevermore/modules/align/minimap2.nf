process minimap2_align {
	container "docker://registry.git.embl.de/schudoma/minimap2-docker:latest"
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)
	val(do_name_sort)

	output:
	tuple val(sample), path("${sample.id}/${sample.id}.sam"), emit: sam

	script:
	def reads = (sample.is_paired) ? "${sample.id}_R1.fastq.gz ${sample.id}_R2.fastq.gz" : "${sample.id}_R1.fastq.gz"
	def threads = task.cpus.intdiv(2)
	def mm_options = "--sam-hit-only -t ${threads} -x sr --secondary=yes -a"

	def sort_cmd = (do_name_sort) ? "samtools collate -@ ${threads} -o ${sample.id}.bam - tmp/collated_bam" : "samtools sort -@ ${threads} -o ${sample.id}.bam -"


	"""
	mkdir -p ${sample.id}/
	minimap2 ${mm_options} --split-prefix ${sample.id}_split ${reference} ${reads} | ${sort_cmd} > ${sample.id}/${sample.id}.bam
	"""
}