process samtools_coverage {
	container "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"
	tag "${sample.id}"	

	input:
	tuple val(sample), path(bam)

	script:
	"""
	mkdir -p ${sample.id}/samtools_coverage/

	samtools coverage ${bam} > ${sample.id}/samtools_coverage/${sample.id}.coverage.txt
	"""
}