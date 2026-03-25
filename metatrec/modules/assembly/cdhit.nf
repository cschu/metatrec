process cd_hit_est {
	container "quay.io/biocontainers/cd-hit:4.8.1--h5ca1c30_13"
	tag "${sample.id}"
	cpus 8
	memory {64.GB * task.attempt}
	time {4.h * task.attempt}
	
	input:
	tuple val(sample), path(fastas)

	output:
	tuple val(sample), path("${sample.id}/${sample.id}.reduced.ffn.gz"), emit: contigs
	
	script:
	def mem = task.memory.toMega()

	"""
	mkdir -p ${sample.id}

	zcat ${fastas} > seqs.ffn

	cd-hit-est -T ${task.cpus} -M ${mem} -i seqs.ffn -o ${sample.id}.cdhit -c 0.95 -n 11 -s 0.90 -d 0 -g 1 -G 1

	gzip -c ${sample.id}.cdhit > ${sample.id}/${sample.id}.reduced.ffn.gz

	rm -fv seqs.ffn
	"""

}