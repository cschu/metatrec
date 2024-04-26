process stringtie {
	container "docker://quay.io/biocontainers/stringtie:2.2.2--h43eeafb_0"

	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), path("${sample.id}/stringtie/${sample.id}.stringtie-transcripts.gtf")

	script:
	"""
	mkdir -p ${sample.id}/stringtie/

	stringtie -o ${sample.id}/stringtie/${sample.id}.stringtie-transcripts.gtf ${bam}
	"""
}