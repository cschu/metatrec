#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_align } from "./nevermore/workflows/align"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/collate"
include { kallisto_index; kallisto_quant} from "./nevermore/modules/profilers/kallisto"


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
			return tuple(meta, file)
		}
	
	annotation_ch.dump(pretty: true, tag: "annotation_ch")

	fastq_ch = fastq_input.out.fastqs

	// fastq_ch.dump(pretty: true, tag: "fastq_ch")

	kallisto_index(
		annotation_ch
	)
	
	// nevermore_main(fastq_ch)

	// counts_ch = nevermore_main.out.readcounts

	// if (do_preprocessing && params.run_qa) {
	// 	collate_stats(counts_ch.collect())		
	// }

}
