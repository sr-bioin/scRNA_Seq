library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(presto)
library(destiny)


# Load sample
# Load count matrix
counts <- readMM("data/DS1/matrix.mtx.gz")
barcodes <- read.table("data/DS1/barcodes.tsv.gz", stringsAsFactors=FALSE)[,1]
features <- read.csv("data/DS1/features.tsv.gz", stringsAsFactors=FALSE, sep="\t", header=FALSE)

# # Assign row and column names
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes
