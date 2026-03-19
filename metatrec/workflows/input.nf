include { prepare_fastqs } from "../../nevermore/workflows/input"

workflow handle_input {		

	main:
		samples_ch = Channel.empty()

		if (params.samplesheet) {
			samples_ch = Channel
				.fromPath(params.samplesheet)
				.splitCsv(sep: '\t', header: [ "id", "source", "r1", "r2", "singles", "contigs", "genes" ])
				.map { row -> 
					def meta = [:]
					meta.id = row.id
					// meta.source = row.source
					meta.library_source = "metaT"
					def reads = null
					if (row.r1 != null && row.r2 != null) {
						meta.is_paired = true
						meta.library = "paired"
						// reads = [ row.r1.split(",").collect({it -> file(it)}), row.r2.split(",").collect({it -> file(it)}) ]
						// reads = [ row.r1, row.r2 ]
						reads = row.r1.split(",") + row.r2.split(",")
					} else {
						meta.is_paired = false
						meta.library = "single"	
						meta.id += ".singles"					
						if (row.r1 != null) {
							// reads = [ row.r1.split(",") ].flatten()
							// reads = [ row.r1 ]
							reads = row.r1.split(",")
						} else if (row.r2 != null) {
							// reads = [ row.r2.split(",") ].flatten()
							// reads = [ row.r2 ]
							reads = row.r2.split(",")
						} else if (row.singles != null) {
							// reads = [ row.singles.split(",") ].flatten()
							// reads = [ row.singles ]
							reads = row.singles.split(",")
						}
					}
					return [ meta, row.source, reads, row.contigs, row.genes ]
				}
				.filter { it[2] != null }
			
			samples_ch.dump(pretty: true, tag: "handle_input_samples_ch")

			prepare_fastqs(
				samples_ch.map { meta, source, reads, contigs, genes -> [ meta.id, reads, params.remote_input_dir, null ] }
			)

			prepare_fastqs.out.pairs.mix(prepare_fastqs.out.singles).dump(pretty: true, tag: "prepare_fastqs")

			samples_ch = samples_ch.map { sample -> [ sample[0].id, sample ] }
				.combine(prepare_fastqs.out.pairs.mix(prepare_fastqs.out.singles), by: 0)
				.map { sample_id, sample, reads -> [ sample[0], sample[1], reads, sample[3], sample[4] ] }
			samples_ch.dump(pretty: true, tag: "samples_ch")
	
		} else {

			// fastq_input(
			// 	Channel.fromPath(input_dir + "/*", type: "dir")
			// 		.filter { !params.ignore_dirs.split(",").contains(it.name) },
			// 	Channel.of(null)
			// )

			// fastq_ch = fastq_input.out.fastqs
			// 	.map { sample, files -> 
			// 		def meta = [:]
			// 		meta.id = sample.id
			// 		meta.sample_id = sample.id.replaceAll(/_.*/, "")
			// 		meta.library_source = "metaT"
			// 		meta.is_paired = sample.is_paired
			// 		meta.library = sample.library
					
			// 		return [ sample.id, meta, files ]
			// 	}
	
			// fastq_ch.dump(pretty: true, tag: "fastq_ch")

			// genes_ch = Channel.fromPath(params.annotation_input_dir + "/**.{fna,ffn}.gz")
			// 	.map { file -> [ file.name.replaceAll(/\.psa_megahit.prodigal.fna.gz$/, ""), file ] }
			// 	// .map { sample_id, file -> [ sample_id.replaceAll(/_.*/, ""), file ] }
			// 	// .map { sample_id, file -> 
			// 	// 	def meta = [:]
			// 	// 	meta.id = sample_id
			// 	// 	meta.sample_id = sample_id.replaceAll(/_.*/, "")
			// 	// 	return [ meta, file ]
			// 	// }

			// contigs_ch = Channel.fromPath(params.assembly_input_dir + "/**.fa.gz")
			// 	.map { file -> [ file.name.replaceAll(/-assembled.fa.gz$/, ""), file ] }
			// 	// .map { sample_id, file -> [ sample_id.replaceAll(/_.*/, ""), file ] }
			// 	// .map { sample_id, file ->
			// 	// 	def meta = [:]
			// 	// 	meta.id = sample_id
			// 	// 	meta.sample_id = sample_id.replaceAll(/_.*/, "")
			// 	// 	return [ meta, file ]
			// 	// }
			
			// samples_ch = fastq_ch
			// 	.join(contigs_ch, by: 0)
			// 	.join(genes_ch, by: 0)
			// 	.map { sample_id, meta, reads, contigs, genes -> [ meta, reads, contigs, genes ] }
		}

	emit:
		samples = samples_ch

}