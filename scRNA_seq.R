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

HCO_dat <- CreateSeuratObject(counts, project="DS1")
HCO_dat <- CreateSeuratObject(counts, project = "DS1", min.cells = 3, min.features = 200)

# Quality control
HCO_dat[["percent.mt"]] <- PercentageFeatureSet(HCO_dat, pattern = "^MT[-\\.]")
png(file="percentagefeature.png", width=800, height=600) # Save figures in working directory
VlnPlot(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()  # Save figures in working directory
VlnPlot(HCO_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# Scatter plots -- library(patchwork)
plot1 <- FeatureScatter(HCO_dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HCO_dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter low-quality cells
HCO_dat <- subset(HCO_dat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

# Normalization
HCO_dat <- NormalizeData(HCO_dat)
HCO_dat <- FindVariableFeatures(HCO_dat, selection.method = "vst", nfeatures = 2000)

# Plot top variable genes
top10 <- head(VariableFeatures(HCO_dat), 10)
plot1 <- VariableFeaturePlot(HCO_dat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

# Scaling and regression
HCO_dat <- ScaleData(HCO_dat)
HCO_dat <- ScaleData(HCO_dat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

# Linear dimension reduction and principal component analysis (PCA)
HCO_dat <- RunPCA(HCO_dat, npcs = 50)
ElbowPlot(HCO_dat, ndims = ncol(Embeddings(HCO_dat, "pca")))

# Heatmap for PCs
PCHeatmap(HCO_dat, dims = 1:20, cells = 200, balanced = TRUE, ncol = 4)

# Run UMAP and tSNE
HCO_dat <- RunTSNE(HCO_dat, dims = 1:20)
HCO_dat <- RunUMAP(HCO_dat, dims = 1:20)

plot1 <- TSNEPlot(HCO_dat) + ggtitle("t-SNE")
plot2 <- UMAPPlot(HCO_dat) + ggtitle("UMAP")
plot1 + plot2

# Cell Type Marker Visualization
plot1 <- FeaturePlot(HCO_dat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), 
                     reduction = "tsne", ncol = 3)
plot2 <- FeaturePlot(HCO_dat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), 
                     reduction = "umap", ncol = 3)
plot1
plot2

#Clustering cells
HCO_dat <- FindNeighbors(HCO_dat, dims = 1:20)
HCO_dat <- FindClusters(HCO_dat, resolution = 1)

# Visualize clusters
DimPlot(HCO_dat, reduction = "tsne", label = TRUE) + 
  DimPlot(HCO_dat, reduction = "umap", label = TRUE)

# Annotation cell clusters
# Known marker genes
ct_markers <- c("MKI67","NES","DCX","FOXG1", # Proliferative, NPC, Neuron, Telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # Ventral telencephalon
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # Dorsal telencephalon
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # Non-telencephalon

DoHeatmap(HCO_dat, features = ct_markers) + NoLegend()

# Cluster mark identification
# With Seurat
cl_markers <- FindAllMarkers(HCO_dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(HCO_dat, features = top10_cl_markers$gene) + NoLegend()

# With presto (fast)
cl_markers_presto <- wilcoxauc(HCO_dat)
cl_markers_presto %>%
  filter(logFC > log(1.2) & pct_in > 0.2 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group = TRUE) %>%
  top_n(n = 2, wt = logFC)

# Top 10 markers heatmap
top10_cl_markers_presto <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(HCO_dat, features = top10_cl_markers_presto$gene) + NoLegend()

#Exploring Specific genes
plot1 <- FeaturePlot(HCO_dat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(HCO_dat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))
