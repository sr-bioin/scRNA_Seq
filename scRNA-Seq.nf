#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.outdir = "./data"
params.script = "bin/ssRNA_Seurat.R"

process RUN_SEURAT {
    tag "$sample_id"
    publishDir "results/${sample_id}", mode: 'copy'

    input:
    tuple path(matrix_dir), val(sample_id)

    output:
    path "*.rds"
    path "*.pdf"

    script:
    """
	Rscript ${params.script} $matrix_dir $sample_id
    """
}

workflow {
    samples_ch = Channel
        .fromPath("${params.outdir}/*", type: 'dir')
        .filter { dir ->
            file("${dir}/matrix.mtx.gz").exists() &&
            file("${dir}/barcodes.tsv.gz").exists() &&
            file("${dir}/features.tsv.gz").exists()
        }
        .map { path_obj -> 
            def sample_id = path_obj.getName()
            tuple(path_obj, sample_id)
        }

    RUN_SEURAT(samples_ch)
}
