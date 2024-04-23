process hisat2_build {
    container "docker://quay.io/biocontainers/hisat2:2.2.1--hdbdd923_6"

    input:
    tuple val(sample), path(genomeseq)

    output:
    tuple val(sample), path("${sample.id}/hisat2/${sample.id}*")

    script:
    """
    mkdir -p ${sample.id}/hisat2/

    hisat2-build -f ${genomeseq} ${sample.id}/hisat2/${sample.id}
    """



}