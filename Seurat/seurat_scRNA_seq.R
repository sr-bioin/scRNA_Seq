#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
matrix_dir <- args[1]
sample_id <- args[2]

# Seurat sCRNA analysis
install.packages("Seurat")
install.packages("Matrix")
install.packages("presto")
install.packages("destiny")
install.packages("patchwork")
install.packages("dplyr")

library(Seurat)
library(Matrix)
library(patchwork)
library(dplyr)
library(presto)
library(destiny)
library(ggplot2)
library(voxhunt)

# Public scRNA-seq data of human cerebral organoids were used [paper](https://www.nature.com/articles/s41586-019-1654-9)
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
# Save figures in working directory
VlnPlot_QC(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", 
                             "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave(file="Quality_control.pdf", width=15, height=10)
dev.off()  # Save figures in working directory

# Scatter plots -- library(patchwork)
plot1_Scp <- FeatureScatter(HCO_dat, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  theme(legend.position="none")
plot2_Scp <- FeatureScatter(HCO_dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position="none")
plot1_Scp + plot2_Scp

# Filter low-quality cells
HCO_dat <- subset(HCO_dat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
plot_filter <- VlnPlot(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
plot_filter

# Normalization
HCO_dat <- NormalizeData(HCO_dat, normalization.method = "LogNormalize", scale.factor = 10000)
HCO_dat <- NormalizeData(HCO_dat)

# Plot top variable genes
HCO_dat <- FindVariableFeatures(HCO_dat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HCO_dat), 10)
plot1_T10 <- VariableFeaturePlot(HCO_dat)
plot2_T10 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, 
                     xnudge = 0, ynudge = 0)
plot1_T10 + plot2_T10

# Scaling and regression
HCO_dat <- ScaleData(HCO_dat)
HCO_dat <- ScaleData(HCO_dat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

# Linear dimension reduction and principal component analysis (PCA)
HCO_dat <- RunPCA(HCO_dat, npcs = 50)
plot_LinD <- ElbowPlot(HCO_dat, ndims = ncol(Embeddings(HCO_dat, "pca")))
plot_LinD

# Heatmap for PCs
plot_PCA <- Heatmap(HCO_dat, dims = 1:21, cells = 200, balanced = TRUE)
plot_PCA

#Clustering
HCO_dat <- FindNeighbors(HCO_dat, dims = 1:10)
HCO_dat <- FindClusters(HCO_dat, resolution = 0.5)

# Non-linear dimensional
# Run UMAP and tSNE
HCO_dat <- RunTSNE(HCO_dat, dims = 1:20)
HCO_dat <- RunUMAP(HCO_dat, dims = 1:20)

plot1_tsne <- TSNEPlot(HCO_dat, reduction = "tsne", label = TRUE) + ggtitle("t-SNE")
plot2_umap <- UMAPPlot(HCO_dat, reduction = "umap", label = TRUE) + ggtitle("umap")
plot1_tsne + plot2_umap

# Cell Type Marker Visualization
plot1_ctMarker <- FeaturePlot(HCO_dat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), 
                     reduction = "tsne", ncol = 3)
plot2_ctMarker <- FeaturePlot(HCO_dat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), 
                     reduction = "umap", ncol = 3)
plot1_ctMarker
plot2_ctMarker


# Annotation cell clusters
# Known marker genes
ct_markers <- c("MKI67","NES","DCX","FOXG1","DLX2","DLX5","ISL1","SIX3","EMX1","PAX6","GLI3","EOMES",#"NEUROD6","RSPO3","OTX2","LHX9","TFAP2A") 
plot1_knMarker <- DoHeatmap(HCO_dat, features = ct_markers) #+ NoLegend()
plot1_knMarker

# Cluster mark identification
# With Seurat
cl.markers <- FindAllMarkers(HCO_dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cl.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
plot_clmarker <- DoHeatmap(HCO_dat, features = top10$gene) + NoLegend()
plot_clmarker
                
# With presto (faster than seurat)
cl_markers_presto <- wilcoxauc(HCO_dat)
cl_markers_presto %>%
  filter(logFC > log(1.2) & pct_in > 0.2 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group = TRUE) %>%
  top_n(n = 2, wt = logFC)

# Top 10 markers heatmap
top10_cl_markers_presto <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
plot_clmarker_presto <- DoHeatmap(HCO_dat, features = top10_cl_markers_presto$gene) + NoLegend()

#Exploring Specific genes
plot1_spgene <- FeaturePlot(HCO_dat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2_spgene <- VlnPlot(HCO_dat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1_spgene + plot2_spgene # + plot_layout(widths = c(1, 2))

# Renaming cell cluster labels by the annotation
new_names <- setNames(c("Dorsal telen. NPC",
                        "Midbrain-hindbrain boundary neuron",
                        "Dorsal telen. neuron",
                        "Dien. and midbrain excitatory neuron",
                        "MGE-like neuron","G2M dorsal telen. NPC",
                        "Dorsal telen. IP","Dien. and midbrain NPC",
                        "Dien. and midbrain IP and excitatory early neuron",
                        "G2M Dien. and midbrain NPC",
                        "G2M dorsal telen. NPC",
                        "Dien. and midbrain inhibitory neuron",
                        "Dien. and midbrain IP and early inhibitory neuron",
                        "Ventral telen. neuron",
                        "Unknown 1",
                        "Unknown 2",
                        "Unknown 3"))
names(new_names) <- levels(HCO_dat)
HCO_dat <- RenameIdents(HCO_dat, new_names)
plot_newIdent <- DimPlot(HCO_dat, reduction = "umap", label = TRUE) + NoLegend()
plot_newIdent
                
# Pseudotemporal cell ordering and variable genes
HCO_dat_dorsal <- subset(HCO_dat, subset = RNA_snn_res.1 %in% c(0, 2, 5, 6, 10))
HCO_dat_dorsal <- FindVariableFeatures(HCO_dat_dorsal, nfeatures = 2000)

# Remove cell cycle genes; cc.genes list is automatically imported by Seurat
VariableFeatures(HCO_dat) <- setdiff(VariableFeatures(HCO_dat), unlist(cc.genes))

# PCA and UMAP before cell cycle correction
HCO_dat_dorsal <- RunPCA(HCO_dat_dorsal) %>% RunUMAP(dims = 1:20)
plot_PCAUMAP_bcc <- FeaturePlot(HCO_dat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)
plot_PCAUMAP_bcc

# Generate cell-cycle-related scores for every cell
HCO_dat_dorsal <- CellCycleScoring(HCO_dat_dorsal,
                                   g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes, 
                                  set.ident = TRUE)
# Regress out cell cycle
HCO_dat_dorsal <- ScaleData(HCO_dat_dorsal, vars.to.regress = c("S.Score", 
                                                                "G2M.Score"))
# Run PCA and UMAP after cell cycle correction
HCO_dat_dorsal <- RunPCA(HCO_dat_dorsal) %>% RunUMAP(dims = 1:20)
plot_PCAUMAP_acc <- FeaturePlot(HCO_dat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)
plot_PCAUMAP_acc

# Run diffusion map to order the cell
# Create diffusion map
diff_map <- DiffusionMap(Embeddings(HCO_dat_dorsal, "pca")[,1:20])

# Compute diffusion pseudotime (DPT)
dpt <- DPT(dm)

# Rank pseudotime values
HCO_dat_dorsal$dpt <- rank(dpt$dpt)

# Plot pseudotime and gene expression
plot_diff_map <- FeaturePlot(HCO_dat_dorsal, c("dpt","GLI3","EOMES","NEUROD6"), ncol=4)
plot_diff_map


# Gene Expression vs Pseudotime
if (is(HCO_dat_dorsal[['RNA']], 'Assay5')){
  expr <- LayerData(HCO_dat_dorsal, assay = "RNA", layer = "data")
} else{
  expr <- HCO_dat_dorsal[['RNA']]@data
}

# Plot expression trends across pseudotime
plot_GLI_Exp <- qplot(HCO_dat_dorsal$dpt, as.numeric(expr["GLI3",]),
               xlab="Dpt", ylab="Expression", main="GLI3") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot_EOM_Exp <- qplot(HCO_dat_dorsal$dpt, as.numeric(expr["EOMES",]),
               xlab="Dpt", ylab="Expression", main="EOMES") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot_NEUR_Exp <- qplot(HCO_dat_dorsal$dpt, as.numeric(expr["NEUROD6",]),
               xlab="Dpt", ylab="Expression", main="NEUROD6") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot_GLI_Exp + plot_EOM_Exp + plot_NEUR_Exp

# Save the results
saveRDS(HCO_dat, file="DS1/HCO_dat_obj_all.rds")
saveRDS(HCO_dat_dorsal, file="DS1/HCO_dat_obj_dorsal.rds")
HCO_dat <- readRDS("DS1/HCO_dat_obj_all.rds")
HCO_dat_dorsal <- readRDS("DS1/HCO_dat_obj_dorsal.rds")
