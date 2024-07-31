process quast {
	container "quay.io/biocontainers/quast:4.1--py27_0"
	publishDir "${params.output_dir}"

	input:
	tuple val(sample), path(assembly)

	script:
	"""

	if [[ ! -f ${sample.id}.${assembly} ]]; then
		ln -s ${assembly} ${sample.id}.${assembly}
	fi

	quast -t ${task.cpus} -o ${sample.id} ${sample.id}.${assembly}
	"""


}