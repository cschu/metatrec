include { bowtie2_build; bowtie2_align } from "../../nevermore/modules/align/bowtie2"
include { hisat2_build; hisat2_align } from "../../nevermore/modules/align/hisat2"
include { merge_and_sort } from "../../nevermore/modules/align/helpers"


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
		hisat2_build(
			samples_by_domain_ch.euk
				.map { meta, source, reads, contigs, genes -> 					
					return [ meta.id.replaceAll(/\.singles$/, ""), contigs ] 
				}
				.unique()
				.map { sample_id, contigs ->
					def meta = [:]
					meta.id = sample_id
					return [ meta, contigs ]
				}
		)
		
		hisat2_build.out.index.dump(pretty: true, tag: "hisat2_build_ch")
		hisat2_input_ch = samples_by_domain_ch.euk
			.map { meta, source, reads, contigs, genes -> [ meta.id.replaceAll(/\.singles$/, ""), meta, reads ] }
			.combine(
				hisat2_build.out.index.map { sample, index -> [ sample.id, index ] },
				by: 0
			)
			.map { sample_id, meta, reads, index -> [ meta, reads, index ]}
		hisat2_input_ch.dump(pretty: true, tag: "hisat2_input_ch")

		hisat2_align(hisat2_input_ch)

		// Prokaryotes
		bowtie2_build(
			samples_by_domain_ch.prok
				.map { meta, source, reads, contigs, genes -> 					
					return [ meta.id.replaceAll(/\.singles$/, ""), contigs ] 
				}
				.unique()
				.map { sample_id, contigs ->
					def meta = [:]
					meta.id = sample_id
					return [ meta, contigs ]
				}
		)

		bowtie2_build.out.index.dump(pretty: true, tag: "bowtie2_build_ch")
		bowtie2_input_ch = samples_by_domain_ch.prok
			.map { meta, source, reads, contigs, genes -> [ meta.id.replaceAll(/\.singles$/, ""), meta, reads ] }
			.combine(
				bowtie2_build.out.index.map { sample, index -> [ sample.id, index ] },
				by: 0
			)
			.map { sample_id, meta, reads, index -> [ meta, reads, index ]}
		bowtie2_input_ch.dump(pretty: true, tag: "bowtie2_input_ch")

		bowtie2_align(bowtie2_input_ch)
		
		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = hisat2_align.out.bam.mix(bowtie2_align.out.bam)
			.map { sample, bam ->
				def meta = sample.clone()
				meta.id = meta.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return [ meta.id, meta, bam ]
			}
			.groupTuple(by:0, size: 2, remainder: true)
			.map { sample_id, samples, bamfiles -> 
				def meta = [:]
				meta.id = sample_id
				meta.library_source = samples[0].library_source
				meta.library = samples[0].library
				return [ meta, bamfiles ]
			}

		aligned_ch.dump(pretty: true, tag: "aligned_ch")

		merge_and_sort(aligned_ch, (params.do_name_sort != null && params.do_name_sort))

	emit:
		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}