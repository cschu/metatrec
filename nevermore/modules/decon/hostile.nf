
params.hostile = [:]
params.hostile.aligner = "bowtie2"

process hostile {
    container "docker://quay.io/biocontainers/hostile:1.1.0--pyhdfd78af_0"

    input:
    tuple val(sample), path(fastqs)
	path(db)

    script:

    def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )

    def r2_input = ""
    if (r2_files.size() != 0) {
        r2_input = "--fastq2 ${r2_files[0]}"
    }

    """
    mkdir -p no_host/${sample.id}

    hostile clean --fastq1 ${r1_files[0]} ${r2_input} --aligner ${params.hostile.aligner} --index ${db} --threads ${task.cpus} --out-dir no_host/${sample.id} --force
    """
}