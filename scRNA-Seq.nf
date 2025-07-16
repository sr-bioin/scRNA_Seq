#!/usr/bin/env nextflow

/*
 * Single cell RNA-seq Analysis Pipeline
 * 
 */

nextflow.enable.dsl = 2

params.reads = "data/*_R{1,2}_001.fastq.gz"
params.outdir = "results"
params.genomeDir = "references/genomeDir"
params.max_cpus = 8
params.max_memory = '32.GB'
params.max_time = '12.h'

log.info """
Single-cell RNA-seq Pipeline
============================
Reads          : ${params.reads}
Reference      : ${params.genomeDir}
Outdir         : ${params.outdir}
"""

Channel
    .fromFilePairs(params.reads, flat: true)
    .set { read_pairs }

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads.join(' ')}
    """
}

process star_solo {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_Solo", emit: solo_out

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${params.genomeDir} \\
         --readFilesIn $r1 $r2 \\
         --readFilesCommand zcat \\
         --soloType CB_UMI_Simple \\
         --soloCBstart 1 --soloCBlen 16 \\
         --soloUMIstart 17 --soloUMIlen 10 \\
         --soloFeatures Gene \\
         --outFileNamePrefix ${sample_id}_Solo/
    """
}

process seurat_r_analysis {
    publishDir "${params.outdir}/seurat", mode: 'copy'

    input:
    path solo_output_dir

    output:
    path "seurat_output/*"

    script:
    """
    Rscript run_seurat.R ${solo_output_dir}
    """
}

workflow {
    fastqc(read_pairs)
    star_solo(read_pairs)
    seurat_r_analysis(star_solo.out.solo_out)
}
