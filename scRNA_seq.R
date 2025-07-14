library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(presto)
library(destiny)
library(ggplot2)

# Public scRNA-seq data of human cerebral organoids were used [paper](https://www.nature.com/articles/s41586-019-1654-9).
# Load sample
# Load count matrix
counts <- readMM("data/DS1/matrix.mtx.gz")
barcodes <- read.table("data/DS1/barcodes.tsv.gz", stringsAsFactors=FALSE)[,1]
features <- read.csv("data/DS1/features.tsv.gz", stringsAsFactors=FALSE, sep="\t", header=FALSE)

# # Assign row and column names
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes

HCO_dat <- CreateSeuratObject(counts, project="DS1")
HCO_dat <- CreateSeuratObject(counts, project = "DS1", min.cells = 3, min.features = 200)

# Quality control
HCO_dat[["percent.mt"]] <- PercentageFeatureSet(HCO_dat, pattern = "^MT[-\\.]")
ggsave(file="percentagefeature.png", width=800, height=600) # Save figures in working directory
Plot(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()  # Save figures in working directory
VlnPlot(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)


