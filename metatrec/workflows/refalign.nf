include { bowtie2_build; bowtie2_align } from "../../nevermore/modules/align/bowtie2"
include { hisat2_build; hisat2_align } from "../../nevermore/modules/align/hisat2"
include { merge_and_sort } from "../../nevermore/modules/align/helpers"
include { prepare_fastqs } from "../../nevermore/workflows/input"


workflow align_to_reference {

	take:
		samples_ch

	main:

		samples_ch.branch { sample ->
			euk: sample[1] == "eukaryote"
			prok: sample[1] == "prokaryote"
		}
		.set { samples_by_domain_ch }

		// Eukaryotes
		hisat2_build(samples_by_domain_ch.euk.map { meta, source, reads, contigs, genes -> [ meta, contigs ] })
		
		hisat2_build.out.index.dump(pretty: true, tag: "hisat2_build_ch")
		hisat2_input_ch = samples_by_domain_ch.euk
			.map { meta, source, reads, contigs, genes -> [ meta.id, meta, reads ] }
			.combine(
				hisat2_build.out.index.map { sample, index -> [ sample.id, index ] },
				by: 0
			)
			.map { sample_id, meta, reads, index -> [ meta, reads, index ]}
		hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")

		hisat2_align(hisat2_input_ch)

		// Prokaryotes
		bowtie2_build(samples_by_domain_ch.prok.map { meta, source, reads, contigs, genes -> [ meta, contigs ] })

		bowtie2_build.out.index.dump(pretty: true, tag: "bowtie2_build_ch")
		bowtie2_input_ch = samples_by_domain_ch.prok
			.map { meta, source, reads, contigs, genes -> [ meta.id, meta, reads ] }
			.combine(
				bowtie2_build.out.index.map { sample, index -> [ sample.id, index ] },
				by: 0
			)
			.map { sample_id, meta, reads, index -> [ meta, reads, index ]}
		bowtie2_input_ch.dump(pretty: true, tag: "bowtie2_input_ch")

		bowtie2_align(bowtie2_input_ch)

		
		downstream_fq_ch  //nevermore_main.out.fastqs
			// .map { sample, fastqs -> return tuple(sample.id.replaceAll(/\.singles$/, ""), sample, fastqs) }
			.map { sample, fastqs -> return tuple(sample.id, sample, fastqs) }
			.combine(
				hisat2_build.out.index
					.map { sample, index -> return tuple(sample.id, sample, index) },
				by: 0
			)
		hisat2_input_chx.dump(pretty: true, tag: "hisat2_input_chx")

		hisat2_input_ch = hisat2_input_chx
			.map { sample_id, sample_fq, fastqs, sample_ix, index ->
				def meta = sample_ix.clone()
				meta.id += "." + sample_fq.id.replaceAll(/SAMEA[0-9]+_METAT/, "")
				// if (sample_fq.id.endsWith(".singles")) {
				// 	meta.id += ".singles"
				// }
				meta.sample_id = sample_ix.sample_id
				return tuple(meta, fastqs, index, "hisat2")
			}
		hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")


		/***/

		// downstream_fq_ch = nevermore_main.out.fastqs
			
		// downstream_fq_ch.dump(pretty: true, tag: "nvm_main_out_ch")

		// hisat2_input_ch = Channel.empty()
		// hisat2_input_chx = nevermore_main.out.fastqs
		// 	.map { sample, fastqs -> return tuple(sample.id.replaceAll(/\.singles$/, ""), sample, fastqs) }
		// 	.combine(
		// 		hisat2_build.out.index
		// 			.map { sample, index -> return tuple(sample.id, sample, index) },
		// 		by: 0
		// 	)
		// hisat2_input_chx.dump(pretty: true, tag: "hisat2_input_chx")

		// hisat2_input_ch = hisat2_input_chx
		// 	.map { sample_id, sample_fq, fastqs, sample_ix, index  ->
		// 		def meta = sample_ix.clone()
		// 		// meta.id = sample_ix.id
		// 		if (sample_fq.id.endsWith(".singles")) {
		// 			meta.id += ".singles"
		// 		}
		// 		// meta.id = sample_fq.id
		// 		meta.sample_id = sample_ix.sample_id
		// 		return tuple(meta, fastqs, index, "hisat2")
		// 	}
		// hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")
		/* HISAT2 */
		// hisat2_build.out.index.dump(pretty: true, tag: "hisat2_build_ch")
		// hisat2_input_chx = downstream_fq_ch  //nevermore_main.out.fastqs
		// 	// .map { sample, fastqs -> return tuple(sample.id.replaceAll(/\.singles$/, ""), sample, fastqs) }
		// 	.map { sample, fastqs -> return tuple(sample.id, sample, fastqs) }
		// 	.combine(
		// 		hisat2_build.out.index
		// 			.map { sample, index -> return tuple(sample.id, sample, index) },
		// 		by: 0
		// 	)
		// hisat2_input_chx.dump(pretty: true, tag: "hisat2_input_chx")

		// hisat2_input_ch = hisat2_input_chx
		// 	.map { sample_id, sample_fq, fastqs, sample_ix, index ->
		// 		def meta = sample_ix.clone()
		// 		meta.id += "." + sample_fq.id.replaceAll(/SAMEA[0-9]+_METAT/, "")
		// 		// if (sample_fq.id.endsWith(".singles")) {
		// 		// 	meta.id += ".singles"
		// 		// }
		// 		meta.sample_id = sample_ix.sample_id
		// 		return tuple(meta, fastqs, index, "hisat2")
		// 	}
		// hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")

		// /* BOWTIE2 */
		// bowtie2_build.out.index.dump(pretty: true, tag: "bowtie2_build_ch")
		// bowtie2_input_chx = downstream_fq_ch  //nevermore_main.out.fastqs
		// 	// .map { sample, fastqs -> return tuple(sample.id.replaceAll(/\.singles$/, ""), sample, fastqs) }
		// 	.map { sample, fastqs -> return tuple(sample.id, sample, fastqs) }
		// 	.combine(
		// 		bowtie2_build.out.index
		// 			.map { sample, index -> return tuple(sample.id, sample, index) },
		// 		by: 0
		// 	)
		// bowtie2_input_chx.dump(pretty: true, tag: "bowtie2_input_chx")

		// bowtie2_input_ch = bowtie2_input_chx
		// 	.map { sample_id, sample_fq, fastqs, sample_ix, index ->
		// 		def meta = sample_ix.clone()
		// 		meta.id += "." + sample_fq.id.replaceAll(/SAMEA[0-9]+_METAT/, "")
		// 		meta.id += ".b"
		// 		// if (sample_fq.id.endsWith(".singles")) {
		// 		// 	meta.id += ".singles"
		// 		// }
		// 		meta.sample_id = sample_ix.sample_id + ".b"
		// 		return tuple(meta, fastqs, index, "bowtie2")
		// 	}
		// bowtie2_input_ch.dump(pretty: true, tag: "bowtie2_input_ch")
		/***/

		// fastq_ch.dump(pretty: true, tag: "alt_fastq_ch")

		// fastq_ch
		// 	.branch {
		// 		hisat2: it[3] == "hisat2"
		// 		bowtie2: it[3] == "bowtie2"
		// 	}
		// 	.set { align_input_ch }

		// align_input_ch.hisat2.dump(pretty: true, tag: "align_input_ch.hisat2")
		// align_input_ch.bowtie2.dump(pretty: true, tag: "align_input_ch.bowtie2")

		// hisat2_align(align_input_ch.hisat2.map { sample, fastqs, index, aligner -> return tuple(sample, fastqs, index) })
		// bowtie2_align(align_input_ch.bowtie2.map { sample, fastqs, index, aligner -> return tuple(sample, fastqs, index) })
		
		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = hisat2_align.out.bam.mix(bowtie2_align.out.bam)
			.map { sample, bam ->
				def meta = sample.clone()
				meta.id = meta.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(meta, bam)
			}
			.groupTuple(sort: true, size: 2, remainder: true)

		aligned_ch.dump(pretty: true, tag: "aligned_ch")

		merge_and_sort(aligned_ch, (params.do_name_sort != null && params.do_name_sort))

	emit:
		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}