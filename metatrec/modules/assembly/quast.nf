process quast {
	container "quay.io/biocontainers/quast:4.1--py27_0"
	publishDir "${params.output_dir}"

	input:
	tuple val(sample), path(assembly)

	script:
	"""

	quast -o ${sample.id} ${assembly}
	"""


}