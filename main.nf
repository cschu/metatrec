#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_align } from "./nevermore/workflows/align"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/stats"
include { kallisto_index; kallisto_quant} from "./nevermore/modules/profilers/kallisto"
include { qc_bbmerge_insert_size } from "./nevermore/modules/qc/bbmerge"
include { hisat2_build; hisat2_align } from "./nevermore/modules/align/hisat2"
include { merge_and_sort } from "./nevermore/modules/align/helpers"
include { stringtie } from "./metatrec/modules/assembly/stringtie"
include { picard_insert_size } from "./metatrec/modules/qc/picard"
include { samtools_coverage} from "./metatrec/modules/qc/samtools"
include { bowtie2_build; bowtie2_align } from "./nevermore/modules/align/bowtie2"
// include { motus; motus_merge } from "./nevermore/modules/profilers/motus"
include { metaT_megahit; bwa_index; bwa2assembly } from "./metatrec/modules/assembly/megahit"
include { metaT_trinity } from "./metatrec/modules/assembly/trinity"
include { quast } from "./metatrec/modules/assembly/quast"

include { align_to_reference } from "./metatrec/workflows/refalign"
include { handle_input } from "./metatrec/workflows/input"


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
params.do_name_sort = false


workflow {

	handle_input()

	samples_ch = handle_input.out.samples
	// sample: [ meta, source, reads, row.contigs, row.genes ]
	fastq_ch = samples_ch.map { meta, source, reads, contigs, genes -> [ meta, reads ] }
	fastq_ch.dump(pretty: true, tag: "fastq_ch")

	genes_ch = samples_ch.map { meta, source, reads, contigs, genes ->  [ meta, genes ] }

	nevermore_main(fastq_ch)

	if (!params.preprocessing_only) {	

		qc_bbmerge_insert_size(fastq_ch)		

		kallisto_index(genes_ch)
		kallisto_index.out.index.dump(pretty: true, tag: "kallisto_index")

		kallisto_quant_input_ch = fastq_ch
			.map { sample, fastqs -> [ sample.id, sample, fastqs ] }
			.combine(
				kallisto_index.out.index
					.map { sample, index -> return [ sample.id, sample, index ] },
				by: 0
			)
			.map { sample_id, sample_fq, fastqs, sample_ix, index  ->
				def meta = sample_fq.clone()
				meta.id = sample_ix.id
				meta.sample_id = sample_ix.sample_id
				return [ meta, fastqs, index ]
			}

		kallisto_quant_input_ch.dump(pretty: true, tag: "kallisto_quant_input_ch")
		
		kallisto_quant(kallisto_quant_input_ch)

		nevermore_main.out.fastqs.dump(pretty: true, tag: "preprocessed_fastqs")
		
		prep_samples_ch = samples_ch
			.map { sample -> [ sample[0].id, sample ] }
			.combine(
				nevermore_main.out.fastqs
					.map { sample, reads -> [ sample.id, sample, reads ] },
				by: 0
			)
			.map { sample_id, sample_raw, sample_prep, reads ->
				return [ sample_prep, sample_raw[1], reads, sample_raw[3], sample_raw[4] ]
			}
			// .map { sample_id, meta, source, raw_reads, contigs, genes, reads ->
			// 	return [ meta, source, reads, contigs, genes ]
			// }
		

		// [DUMP: preprocessed_fastqs] [
		// 	{
		// 		"id": "SAMEA112490973",
		// 		"library_source": "metaT",
		// 		"is_paired": "true",
		// 		"library": "paired",
		// 		"merged": "true"
		// 	},
		// 	[
		// 		"/scratch/schudoma/WORK/metatrec/5c/0fdc4cd4742cbc8a3210d60843fd7b/qc_reads/SAMEA112490973/SAMEA112490973_R1.fastq.gz",
		// 		"/scratch/schudoma/WORK/metatrec/5c/0fdc4cd4742cbc8a3210d60843fd7b/qc_reads/SAMEA112490973/SAMEA112490973_R2.fastq.gz"
		// 	]
		// ]

		// [DUMP: prep_samples_ch] [
		// 	"SAMEA112490973",
		// 	[
		// 		{
		// 			"id": "SAMEA112490973",
		// 			"library_source": "metaT",
		// 			"is_paired": "true",
		// 			"library": "paired"
		// 		},
		// 		"prokaryote",
		// 		[
		// 			"/scratch/schudoma/WORK/metatrec/a9/a819fdc64ccd29d9198415ec053d50/fastq/SAMEA112490973/SAMEA112490973_R1.fastq.gz",
		// 			"/scratch/schudoma/WORK/metatrec/a9/a819fdc64ccd29d9198415ec053d50/fastq/SAMEA112490973/SAMEA112490973_R2.fastq.gz"
		// 		],
		// 		"/g/bork6/schudoma/data_processing/trec/amr/input/assemblies/SAMEA112490973-assembled.fa.gz",
		// 		"/g/bork6/schudoma/data_processing/trec/amr/input/genes/SAMEA112490973.psa_megahit.prodigal.fna.gz"
		// 	],
		// 	[
		// 		"/scratch/schudoma/WORK/metatrec/5c/0fdc4cd4742cbc8a3210d60843fd7b/qc_reads/SAMEA112490973/SAMEA112490973_R1.fastq.gz",
		// 		"/scratch/schudoma/WORK/metatrec/5c/0fdc4cd4742cbc8a3210d60843fd7b/qc_reads/SAMEA112490973/SAMEA112490973_R2.fastq.gz"
		// 	]
		// ]
		
		prep_samples_ch.dump(pretty: true, tag: "prep_samples_ch")

		// align_to_reference(hisat2_input_ch.mix(bowtie2_input_ch))
		align_to_reference(prep_samples_ch)
		
		stringtie(align_to_reference.out.alignments)
		picard_insert_size(
			align_to_reference.out.alignments
				.filter { !it[0].id.endsWith("singles") && !it[0].id.endsWith("singles.b") }
		)
		samtools_coverage(align_to_reference.out.alignments)

		counts_ch = nevermore_main.out.readcounts
		counts_ch = counts_ch.mix(
			align_to_reference.out.aln_counts
				.map { sample, file -> file }
				.collect()
		)

		downstream_fq_ch = prep_samples_ch.map { meta, source, reads, contigs, genes -> [ meta, reads ] }
		
		// motus(nevermore_main.out.fastqs, params.motus_db)
		// motus_merge(
		// 	motus.out.motus_profile
		// 		.map { sample, profile -> return profile }
		// 		.collect(),
		// 	params.motus_db
		// )

		assembly_input_ch = downstream_fq_ch
			.map { sample, fastqs -> 
				def meta = sample.clone()
				meta.id = meta.id.replaceAll(/\.(singles|orphans)$/, "")
				return [ meta.id, meta.sample_id, fastqs ]
			}
			.groupTuple(by: [0, 1], size: 2, remainder: true)
			.map { sample_protocol_id, sample_id, fastqs -> 
				def meta = [:]
				meta.id = sample_protocol_id
				meta.sample_id = sample_id
				return [ meta, [fastqs].flatten() ]
			}

		assembly_input_ch.dump(pretty: true, tag: "assembly_input_ch")

		metaT_megahit(assembly_input_ch, "stage1")

		metaT_trinity(assembly_input_ch, "stage1")

		quast(metaT_megahit.out.contigs.mix(metaT_trinity.out.contigs))

		bwa_index(
			metaT_megahit.out.contigs
				.map { sample, contigs -> 
					def meta = sample.clone()
					sample.assembler = "megahit"
					return [ sample, contigs ]
				}
				.mix(
					metaT_trinity.out.contigs
						.map { sample, contigs -> 
							def meta = sample.clone()
							meta.assembler = "trinity"
							return [ sample, contigs ]
						}
				)		
		)

		bwa2assembly(
			downstream_fq_ch
				.map { sample, fastqs -> [ sample.id.replaceAll(/\.singles$/, ""), sample, fastqs ] }
				.combine(bwa_index.out.index, by: 0)
				.map { sample_id, sample, fastqs, index -> 
					def meta = sample.clone()
					meta.index_id = sample_id
					return [ meta, fastqs, index ]
				}
		)

		if (do_preprocessing && params.run_qa) {
			collate_stats(counts_ch.collect())		
		}

	}

}
