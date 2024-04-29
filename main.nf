#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_align } from "./nevermore/workflows/align"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/collate"
include { kallisto_index; kallisto_quant} from "./nevermore/modules/profilers/kallisto"
include { qc_bbmerge_insert_size } from "./nevermore/modules/qc/bbmerge"
include { hisat2_build; hisat2_align } from "./nevermore/modules/align/hisat2"
include { merge_and_sort } from "./nevermore/modules/align/helpers"
include { stringtie } from "./metatrec/modules/assembly/stringtie"
include { picard_insert_size } from "./metatrec/modules/qc/picard"
include { samtools_coverage} from "./metatrec/modules/qc/samtools"


if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


params.ignore_dirs = ""




workflow align_to_reference {

	take:
		fastq_ch

	main:

		hisat2_align(fastq_ch)
		
		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = hisat2_align.out.bam
			.map { sample, bam ->
				def meta = sample.clone()
				meta.id = meta.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(meta, bam)
			}
			.groupTuple(sort: true)

		aligned_ch.dump(pretty: true, tag: "aligned_ch")

		merge_and_sort(aligned_ch, (params.do_name_sort != null && params.do_name_sort))

	emit:
		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}




workflow {

	fastq_input(
		Channel.fromPath(input_dir + "/*", type: "dir")
			.filter { !params.ignore_dirs.split(",").contains(it.name) },
		Channel.of(null)
	)

	annotation_ch = Channel.fromPath(params.annotation_input_dir + "/**.{fna,ffn}.gz")
		.map { file ->
			// SAMEA112553567_METAG_H5WNWDSXC.UDI027-1.psa_megahit.prodigal.fna.gz
			// return tuple(file.name.replaceAll(/_.*/, ""), file.name.replaceAll(/\.psa_megahit.prodigal.fna.gz$/, "").replaceAll(/^.+_METAG_/, ""), file)
			return tuple(file.name.replaceAll(/\.psa_megahit.prodigal.fna.gz$/, ""), file)
			// return tuple(file.getParent().getName(), file)
		}
		.map { sample_id, file -> 
			def meta = [:]
			meta.id = sample_id
			meta.sample_id = sample_id.replaceAll(/_.*/, "")
			// meta.library_source = "metaT"
			return tuple(meta, file)
		}

	assembly_ch = Channel.fromPath(params.assembly_input_dir + "/**.fa.gz")
		.map { file ->
			//SAMEA112489502_METAG_H5WNWDSXC.UDI049-1-assembled.fa.gz
			return tuple(file.name.replaceAll(/-assembled.fa.gz$/, ""), file)
		}
		.map { sample_id, file ->
			def meta = [:]
			meta.id = sample_id
			meta.sample_id = sample_id.replaceAll(/_.*/, "")
			return tuple(meta, file)
		}
	
	annotation_ch.dump(pretty: true, tag: "annotation_ch")

	fastq_ch = fastq_input.out.fastqs
		.map { sample, file -> 
			def meta = sample.clone()
			meta.library_source = "metaT"
			return tuple(meta, file)
		}
	
	qc_bbmerge_insert_size(fastq_ch)

	fastq_ch.dump(pretty: true, tag: "fastq_ch")

	kallisto_index(
		annotation_ch			
	)
	kallisto_index.out.index.dump(pretty: true, tag: "kallisto_index")

	hisat2_build(
		assembly_ch
	)

	// Apr-23 10:49:27.850 [Actor Thread 4] INFO  nextflow.extension.DumpOp - [DUMP: annotation_ch] [
	// 	{
	// 		"id": "SAMEA112551184_METAG_H5WNWDSXC.UDI026-1",
	// 		"sample_id": "SAMEA112551184"
	// 	},
	// 	"/g/scb/bork/schudoma/tasks/metatrec.1892.20240422/annotation_input/SAMEA112551184_METAG_H5WNWDSXC.UDI026-1.psa_megahit.prodigal.fna.gz"
	// ]


	kallisto_quant_input_ch = fastq_ch
		.map { sample, fastqs -> return tuple(sample.id, sample, fastqs) }
		.combine(
			kallisto_index.out.index
				.map { sample, index -> return tuple(sample.sample_id, sample, index) },
			by: 0
		)
		.map { sample_id, sample_fq, fastqs, sample_ix, index  ->
			def meta = sample_fq.clone()
			meta.id = sample_ix.id
			meta.sample_id = sample_ix.sample_id
			return tuple(meta, fastqs, index)
		}

	kallisto_quant_input_ch.dump(pretty: true, tag: "kallisto_quant_input_ch")
	
	kallisto_quant(kallisto_quant_input_ch)
	
	nevermore_main(fastq_ch)

	
	hisat2_input_chx = nevermore_main.out.fastqs
		.map { sample, fastqs -> return tuple(sample.id.replaceAll(/\.singles$/, ""), sample, fastqs) }
		.combine(
			hisat2_build.out.index
				.map { sample, index -> return tuple(sample.sample_id, sample, index) },
			by: 0
		)
	hisat2_input_chx.dump(pretty: true, tag: "hisat2_input_chx")

	hisat2_input_ch = hisat2_input_chx
		.map { sample_id, sample_fq, fastqs, sample_ix, index  ->
			def meta = sample_ix.clone()
			// meta.id = sample_ix.id
			if (sample_fq.id.endsWith(".singles")) {
				meta.id += ".singles"
			}
			// meta.id = sample_fq.id
			meta.sample_id = sample_ix.sample_id
			return tuple(meta, fastqs, index)
		}
	hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")
	
	align_to_reference(hisat2_input_ch)

	
	stringtie(align_to_reference.out.alignments)
	picard_insert_size(align_to_reference.out.alignments)
	samtools_coverage(align_to_reference.out.alignments)

	counts_ch = nevermore_main.out.readcounts
	// counts_ch = counts_ch.concat(
	// 		align_to_reference.out.aln_counts
	// 			.map { sample, file -> return file }
	// 			.collect()
	// 	)

	if (do_preprocessing && params.run_qa) {
		collate_stats(counts_ch.collect())		
	}

}
