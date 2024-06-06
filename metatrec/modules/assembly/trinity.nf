process metaT_trinity {
	container "docker://quay.io/biocontainers/trinity:2.15.1--pl5321hdcf5f25_4"
	label "trinity"

	input:
	tuple val(sample), path(fastqs)
	val(stage)

	output:
	tuple val(sample), path("assemblies/metaT_trinity/${stage}/${sample.library_source}/${sample.id}/*.fasta"), emit: contigs

	script:
	def mem_gb = task.memory.toGiga()
	def mem = task.memory.toBytes()

	def input_files = ""
	// we cannot auto-detect SE vs. PE-orphan!
	def r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
	def r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
	def orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

	def append_orphans = ""
	def left = ""
	def right = ""

	def input_string = ""

	def make_left = ""
	def make_right = ""
	if (r1_files.size() != 0 && r2_files.size() != 0) {

		input_string = "--left left.fastq --right right.fastq"

		// make_left += "gzip -dc ${r1_files[0]} | awk 'NR % 4 == 1 { \$1=\$1\"/1\" } { print \$0; }' >> left.fastq\n"
		make_left += "gzip -dc ${r1_files[0]} >> left.fastq\n"
			// printf(\">P%s/1\\n%s\\n\", int(NR/4), \$1); }' >> left.fastq\n"

		if (orphans.size() != 0) {
			// make_left += "gzip -dc ${orphans[0]} | awk 'NR % 4 == 1 { \$1=\$1\"/1\" } { print \$0; }' >> left.fastq\n" 
			// make_left += "gzip -dc ${orphans[0]} >> left.fastq\n"
			make_left += "gzip -dc ${orphans[0]} | awk 'NR % 4 == 1 { \$0=\"@orphan\"int(NR/4)/1; } { print \$0; }' >> left.fastq\n"
			// awk 'NR % 4 == 2 { printf(\">O%s/1\\n%s\\n\", int(NR/4), \$1); }' >> left.fastq"
		}

		// make_right += "gzip -dc ${r2_files[0]} | awk 'NR % 4 == 1 { \$1=\$\1\"/2\" } { print \$0; }' >> right.fastq\n"
		make_right += "gzip -dc ${r2_files[0]} >> right.fastq\n"
		// | awk 'NR % 4 == 2 { printf(\">P%s/2\\n%s\\n\", int(NR/4), \$1); }' >> right.fastq\n"		

	} else if (r1_files.size() != 0) {

		input_string = "--single ${r1_files[0]}"
		
	} else if (r1_files.size() != 0) {

		input_string = "--single ${r2_files[0]}"
		
	} else if (r2_files.size() != 0) {

		input_string = "--single ${orphans[0]}"

	}
	def outdir = "assemblies/metaT_trinity/${stage}/${sample.library_source}/${sample.id}"
	
	"""
	mkdir -p ${outdir}/ trinity/

	${make_left}
	${make_right}

	Trinity --seqType fq --max_memory ${mem_gb}G ${input_string} --CPU ${task.cpus} --output trinity/

	mv -v trinity/* ${outdir}/

	"""
	// cp -v megahit_out/final.contigs.fa ${outdir}/${sample.id}.${stage}.transcripts.fasta
	
}