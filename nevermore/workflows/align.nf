#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bwa_mem_align } from "../modules/align/bwa"
include { minimap2_align } from "../modules/align/minimap2"
include { merge_and_sort } from "../modules/align/helpers"

def asset_dir = "${projectDir}/nevermore/assets"
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


workflow nevermore_align {

	take:
		fastq_ch

	main:		
		/*	align merged single-read and paired-end sets against reference */

		minimap2_align(
			fastq_ch,
			params.reference,
			params.do_name_sort
		)

		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = minimap2_align.out.bam
			.map { sample, bam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, bam)
			}
			.groupTuple(sort: true)

		merge_and_sort(aligned_ch, true)

	emit:
		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}
