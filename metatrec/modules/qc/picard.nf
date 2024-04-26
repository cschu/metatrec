process picard_insert_size {
	container "docker://quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
	tag "${sample.id}"

	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), "${sample.id}/picard/${sample.id}.imetrics.txt", emit: isize_metrics
	tuple val(sample), "${sample.id}/picard/${sample.id}.ihist.pdf", emit: isize_hist

	script:
	def picard_params = "MINIMUM_PCT=0.05"
	"""
	mkdir -p ${sample.id}/picard/

	picard CollectInsertSizeMetrics ${picard_params} INPUT=${bam} OUTPUT=${sample.id}/picard/${sample.id}.imetrics.txt HISTOGRAM_FILE=${sample.id}/picard/${sample.id}.ihist.pdf
	"""
}



// rule picard_insert_size:
//     input:
//         bam = "{sample}.{aligner}.bam"
//     output:
//         metrics = "{sample}.{aligner}.insert_size_metrics.txt",
//         histogram = "{sample}.{aligner}.insert_size_histogram.pdf"
//     params:
//         program_call = config["program_calls"]["picard_insert_size"],
//         program_params = config["parameters"]["picard"]["collect_insertsize_metrics"]
//     resources:
//         mem_mb = HPC_CONFIG.get_memory("picard_insert_size")
//     log:
//         "{sample}.{aligner}.picard_insert.log"
//     shell:
//         "(({params.program_call} INPUT={input.bam} OUTPUT={output.metrics} HISTOGRAM_FILE={output.histogram} {params.program_params}" + \
//         " && touch {output.metrics} {output.histogram})" + \
//         " || (if [[ ! -e {input.bam}.bai ]]; then touch {output.metrics} {output.histogram}; fi)) &> {log}"



// 		collect_insertsize_metrics: "MINIMUM_PCT=0.05"
// 		picard CollectInsertSizeMetrics