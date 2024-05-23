params.kmer_steps = "25,29,33,37,41,45,49,53,57,61,65,69,73,77,81,85,89,93,97"


process metaT_megahit {
	container "docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_4"
	label "megahit"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/metaT_megahit/${stage}/${sample.library_source}/${sample.id}/${sample.id}.${stage}.*.fasta"), emit: contigs

	script:
	def mem_gb = task.memory.toGiga()
	def mem = task.memory.toBytes()

	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	if (r1_files.size() != 0) {
		input_files += "-1 ${r1_files.join(' ')}"
	}
	if (r2_files.size() != 0) {
		input_files += " -2 ${r2_files.join(' ')}"
	}
	if (orphans.size() != 0) {
		input_files += " -r ${orphans.join(' ')}"
	}

	def kmer_params = "--k-list ${params.kmer_steps}" //--k-min ${params.mink} --k-max ${params.maxk} --k-step ${params.stepk}"
	def outdir = "assemblies/metaT_megahit/${stage}/${sample.library_source}/${sample.id}"
	
	"""
	mkdir -p ${outdir}/
	megahit -t ${task.cpus} --cpu-only -m ${mem} ${input_files} ${kmer_params} --bubble-level 0 --mem-flag 1
	cp -v megahit_out/final.contigs.fa ${outdir}/${sample.id}.${stage}.transcripts.fasta
	"""
	
}