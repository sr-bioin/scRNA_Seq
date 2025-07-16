#!/usr/bin/env nextflow

/*
 * Single cell RNA-seq Analysis Pipeline
 * 
 */
nextflow.enable.dsl = 2
params.reads = 
params.outdir = "results"
params.genomeDir = "references/genomeDir"
params.max_cpus = 8
params.max_memory = '32.GB'
params.max_time = '12.h'
