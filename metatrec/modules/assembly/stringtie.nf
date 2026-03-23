process stringtie {
	container "quay.io/biocontainers/stringtie:2.2.2--h43eeafb_0"
	tag "${sample.id}"

	input:
	tuple val(sample), path(bam)

	output:
	tuple val(sample), path("${sample.id}/stringtie/${sample.id}.stringtie-transcripts.gtf"), emit: gtf

	script:
	"""
	mkdir -p ${sample.id}/stringtie/

	stringtie -o ${sample.id}/stringtie/${sample.id}.stringtie-transcripts.gtf ${bam}
	"""
}

process extract_stringtie_transcripts {
	container "quay.io/biocontainers/bedtools:2.31.1--h13024bc_3"
	tag "${sample.id}"

	input:
	tuple val(sample), path(gtf), path(fasta)

	output:
	tuple val(sample), path("${sample.id}/stringtie/${sample.id}.stringtie.ffn"), emit: transcripts

	script:
	"""
	mkdir -p ${sample.id}/stringtie/

	if [[ ${fasta} == *.gz ]]; then
		zcat ${fasta} > genome.fa
	else
		ln -s ${fasta} genome.fa
	fi

	awk -v OFS='\\t' '/^#/ { print \$0; next } \$3 == "transcript" { print \$0; }' ${gtf} > transcripts.gtf

	bedtools getfasta -fi genome.fa -bed transcripts.gtf -fullHeader -fo ${sample.id}/stringtie/${sample.id}.stringtie.ffn

	rm -vf genome.fa transcripts.gtf
	"""

	output:

}