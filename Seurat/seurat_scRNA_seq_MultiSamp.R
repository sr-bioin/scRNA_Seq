#===========================================================================================================
#Multiple dataset scRNA-seq analysis

library(Seurat)
library(dplyr)
library(patchwork)
library(harmony)
library(rliger)
library(CSS)
library(SeuratWrappers)
library(gplots)
library(qlcMatrix)

# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

# Merge two data sets
ds1_ds2_merege <- merge(HCO_data_DS1, HCO_data_DS2)%>%
  FindVariableFeatures(nfeatures = 3000)%>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)
plot_ds12_merge1 <- DimPlot(ds1_ds2_merege,group.by="orig.ident")
plot_ds12_merge2 <- FeaturePlot(ds1_ds2_merege, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
plot_ds12_merge1 + plot_ds12_merge2 + plot_layout(widths = c(1.5, 2))

ds2_ds3_merege <- merge(HCO_data_DS2, HCO_data_DS3)%>%
  FindVariableFeatures(nfeatures = 3000)%>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)
plot_ds23_merge1 <- DimPlot(ds2_ds3_merege,group.by="orig.ident")
plot_ds23_merge2 <- FeaturePlot(ds1_ds2_merege, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
plot_ds23_merge1 + plot_ds23_merge2 + plot_layout(widths = c(1.5, 2))

#=======================================================================================
# 1). Seurat integration
HCO_data_DS1 <- NormalizeData(HCO_data_DS1) %>% FindVariableFeatures(nfeatures = 3000)
HCO_data_DS2 <- NormalizeData(HCO_data_DS2) %>% FindVariableFeatures(nfeatures = 3000)

# Merge objects into a list
HCO_seurat_LIST <- list(HCO_data_DS1, HCO_data_DS2)

# Normalize and find variable features
HCO_seurat_LIST <- lapply(seurat_LIST, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  return(x)
})

# Find integration anchors
HCO_anchors <- FindIntegrationAnchors(object.list = HCO_seurat_LIST, dims = 1:30)

# Integrate data, express level correction
HCO_seurat_INTEGRATE <- IntegrateData(anchorset = HCO_anchors, dims = 1:30)

# Downstream analysis
HCO_seurat <- ScaleData(HCO_seurat_INTEGRATE) %>% RunPCA() %>% RunUMAP(dims = 1:30)
HCO_seurat <- RunPCA(HCO_seurat, npcs = 50)
HCO_seurat <- RunUMAP(HCO_seurat, dims = 1:20)
HCO_seurat <- FindNeighbors(HCO_seurat, dims = 1:20) %>% FindClusters(resolution = 0.6)


DefaultAssay(HCO_seurat) <- "RNA"
plot1 <- UMAPPlot(HCO_seurat, group.by="orig.ident")
plot2 <- UMAPPlot(HCO_seurat, label = T)
plot3 <- FeaturePlot(HCO_seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
Plot_seurat_combine <- ((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
Plot_seurat_combine

#------------------------------------------------------------------------------------------
# 2). Harmony integration

Seurat_combine <- merge(HCO_data_DS1, HCO_data_DS2) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

harmony <- RunHarmony(Seurat_combine, group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50)
harmony <- RunUMAP(harmony, reduction = "harmony", dims = 1:20)
harmony <- FindNeighbors(harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.6)

# Save the object
saveRDS(harmony, file="integrated_harmony.rds")

# Visualize the integration results similar to before
plot1 <- UMAPPlot(harmony, group.by="orig.ident")
plot2 <- UMAPPlot(harmony, label = T)
plot3 <- FeaturePlot(harmony, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
plot_harmony_combine <- ((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))
plot_harmony_combine <- 

#===========================================================================================================
# Combine three datasets

library(patchwork)

pbmc.big <- merge(HCO_data_DS1, y = c(HCO_data_DS2, HCO_data_DS3), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
pbmc.big


#==============================================================================
# Analysis of a new data using annotated reference data set

# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

# Reference data set
seurat_ref <- readRDS("data/ref_seurat_obj.rds")
plot1 <- UMAPPlot(seurat_ref, group.by="branch")
plot2 <- UMAPPlot(seurat_ref, group.by="celltype")
plot3 <- FeaturePlot(seurat_ref, c("SOX2","DCX","FOXG1","EMX1","DLX2","LHX9"),
                     ncol=3, pt.size = 0.1)
plot_reference <- ((plot1 / plot2) | plot3) + plot_layout(width = c(1,3))
plot_reference

# Transcriptome similarity on cell cluster level
# Average transcriptome profiles for every annotated cell type in the reference 
# data set and every cell cluster in the query data set.
avg_expr_ref <- sapply(sort(unique(seurat_ref$celltype)), function(ct) 
  rowMeans(seurat_ref@assays$RNA@data[,which(seurat_ref$celltype == ct)] ))

avg_expr_ds1 <- sapply(levels(HCO_data_DS1@active.ident), function(ct) {
  cells <- which(HCO_data_DS1@active.ident == ct)
  rowMeans(LayerData(HCO_data_DS1, layer = "data")[, cells, drop = FALSE])})

# Calculate pairwise Spearman correlation across those genes
genes2cor <- intersect(VariableFeatures(seurat_ref), rownames(HCO_data_DS1))
corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="spearman")

library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=1, margins=c(14,15),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, 
          cexCol=0.8, col=colorRampPalette(rev(c("#b2182b","#d6604d",
           "#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))

#==============================================================================
# 1. Transcriptome similarity on cell level
ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)

library(presto)
ranked_expr_ds1 <- rank_matrix(LayerData(HCO_data_DS1, layer = "data")[genes2cor, ])$X_ranked

# Ranking algorithm by ourselves. 
rank_matrix <- function (mat) 
{
  if (is.matrix(mat) | is.data.frame(mat)) {
    ranked_mat <- apply(mat, 2, rank)
  }
  else {
    df_mat <- Matrix::summary(mat)
    dfs_mat <- split(df_mat, df_mat$j)
    df_mat_ranked <- do.call(rbind, lapply(dfs_mat, function(df) {
      num_zeros <- nrow(mat) - nrow(df)
      ranks_nonzero <- rank(df$x)
      df$x <- ranks_nonzero + num_zeros - (1 + num_zeros)/2
      return(df)
    }))
    ranked_mat <- sparseMatrix(i = df_mat_ranked$i, j = df_mat_ranked$j, 
                               x = df_mat_ranked$x, dims = dim(mat), dimnames = dimnames(mat))
  }
  return(ranked_mat)
}
ranked_expr_ds1 <- rank_matrix(HCO_data_DS1@assays$RNA@data[genes2cor,])
ranked_expr_ds1 <- rank_matrix(LayerData(HCO_data_DS1, layer = "data")[genes2cor, ])


# Calculate Pearson correlation between two sparse matrix or one sparse matrix and one dense matrix,
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
HCO_data_DS1$celltype_maxcor <- ct_maxcor


plot1 <- UMAPPlot(HCO_data_DS1, label=T)
plot2 <- UMAPPlot(HCO_data_DS1, group.by="celltype_maxcor", label=T) + 
  coord_cartesian(clip = "off")
plot1
plot2


# Summarize the cell-level similarities to the query cell clusters
corr2ref_scaled <- scale(t(corr2ref_cell))
corr2ref_sum2cl <- t(sapply(levels(HCO_data_DS1@active.ident), function(cl)
  rowMeans(corr2ref_scaled[,which(HCO_data_DS1@active.ident == cl)]) ))
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, 
          cexCol=0.8,col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582",
          "#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))


# 2. Seurat-based label transfer
# Update reference to v5 format
seurat_ref_v5 <- UpdateSeuratObject(seurat_ref) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

# Query data
HCO_data_DS1 <- NormalizeData(HCO_data_DS1) %>%
  FindVariableFeatures() %>%
  ScaleData()

# Find anchors
anchors <- FindTransferAnchors(
  reference = seurat_ref_v5,
  query = HCO_data_DS1,
  dims = 1:30,
  reference.reduction = "pca"
)

# Transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = seurat_ref_v5$celltype,
  dims = 1:30,
  reference = seurat_ref_v5,  # Explicitly provide reference object
  query = HCO_data_DS1    # Explicitly provide query object
)

# Add predictions to query
HCO_data_DS1$celltype_transfer <- predictions$predicted.id

# Visualize
plot1 <- UMAPPlot(HCO_data_DS1, label=T, label.size = 4, repel=TRUE)
plot2 <- UMAPPlot(HCO_data_DS1, group.by="celltype_transfer", 
                  label=T, label.size = 4, repel=TRUE)
Seurat-based_label_transfer <- plot1 | plot2 
ggsave(file="Seurat-based_label_transfer.pdf", width=18, height=5)
ggsave(file="Seurat-based_label_transfer.jpg", width=16, height=6)
dev.off()  # Save figures in working directory

# ================================================================================
# NOt working need to check
pred_matrix <- as.matrix(GetAssayData(predictions, slot = "data"))
pred_matrix <- as.matrix(GetAssayData(predictions, slot = "data"))

# Check and align identities
ident_vec <- HCO_data_DS1@active.ident
if (!is.factor(ident_vec)) ident_vec <- as.factor(ident_vec)

pred_scores_sum2cl <- do.call(rbind, lapply(levels(ident_vec), function(cl) {
  cells <- which(ident_vec == cl)
  colMeans(pred_matrix[cells, , drop = FALSE])
}))

# Add row and column names for clarity
rownames(pred_scores_sum2cl) <- levels(ident_vec)
colnames(pred_scores_sum2cl) <- colnames(pred_matrix)

library(gplots)

heatmap.2(pred_scores_sum2cl,
          scale = "none",
          trace = "none",
          key = FALSE,
          keysize = 0.5,
          margins = c(15, 17),
          labRow = colnames(avg_expr_ds1),    # should match 11 rows
          labCol = unique(seurat_ref$celltype), # should match 14 columns
          cexRow = 0.8,
          cexCol = 0.8,
          col = colorRampPalette(rev(c(
            "#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7",
            "#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))



# Original one
pred_scores_sum2cl <- t(sapply(levels(HCO_data_DS1@active.ident), function(cl)
  colMeans(predictions[which(HCO_data_DS1@active.ident == cl),-c(1,ncol(predictions))]) ))

heatmap.2(pred_scores_sum2cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = unique(seurat_ref$celltype), 
          cexRow=0.8, cexCol=0.8, col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582",
          "#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))

#===================================================================================
# Cluster connectivity analysis with PAGA pipeline
library(Seurat)
library(Matrix)
library(anndata)
library(loomR)

# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

# Extract cell attributes from DS1
cell_attrs <- list(
  pca = Embeddings(HCO_data_DS1, "pca")[, 1:20],
  umap = Embeddings(HCO_data_DS1, "umap"),
  celltype = HCO_data_DS1@active.ident
)

# Extract data matrices
data_mat <- GetAssayData(HCO_data_DS1, assay = "RNA", layer = "data")
counts_mat <- GetAssayData(HCO_data_DS1, assay = "RNA", layer = "counts")

# Create loom file
loom <- create("data/DS1/HCO_dat_DS1_loom_obj.loom",
               data = data_mat,
               layers = list(counts = counts_mat),
               cell.attrs = cell_attrs)

loom$close_all()

# Prepare AnnData object
library(reticulate)
expr_data <- LayerData(HCO_data_DS1[["RNA"]], layer = "data")
counts_data <- LayerData(HCO_data_DS1[["RNA"]], layer = "counts")

adata <- AnnData(
  X = t(expr_data),
  obs = data.frame(celltype = HCO_data_DS1@active.ident,
                   row.names = colnames(HCO_data_DS1)),
  var = data.frame(row.names = rownames(HCO_data_DS1)),
  layers = list(counts = t(counts_data)),
  obsm = list(
    pca = Embeddings(HCO_data_DS1, "pca")[, 1:20],
    umap = Embeddings(HCO_data_DS1, "umap")
  )
)

adata$write_h5ad("data/DS1/HCO_dat_DS1_anndata_obj.h5ad")

# Activate Python environment and load Scanpy
Sys.setenv(RETICULATE_CONDA = "C:/Users/thapa/miniconda3/Scripts/conda.exe")
use_condaenv("scanpy-env", required = TRUE)
scanpy <- import("scanpy")
py_config()

# Load AnnData object
adata_DS1 <- scanpy$read("data/DS1/HCO_dat_DS1_anndata_obj.h5ad")

# Run PAGA analysis
scanpy$pp$neighbors(adata_DS1, n_neighbors = 20L, use_rep = 'pca')
scanpy$tl$paga(adata_DS1, groups = 'celltype')
adata_DS1$write_h5ad("data/DS1/HCO_dat_DS1_anndata_obj.h5ad")

# Plot PAGA graph
plt <- import("matplotlib")
plt$use("Agg", force = TRUE) # need matplotlib, doesnot work if .pyplot
scanpy$pl$paga(adata_DS1,
               color = 'celltype',
               fontsize = 6,
               frameon = FALSE,
               save = "DS1_paga_NEW.png")

#==================================================================================================
# Pseudotime reconstruction 

library(Seurat)
library(destiny)

# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

dm <- DiffusionMap(Embeddings(HCO_data_DS1, "pca")[,1:20])
suppressWarnings({
  dm <- DiffusionMap(Embeddings(HCO_data_DS1, "pca")[, 1:20])
})
dpt <- DPT(dm)
HCO_data_DS1$dpt <- rank(dpt$dpt)

FeaturePlot(HCO_data_DS1, c("dpt","SOX2","NHLH1","DCX"), ncol=4)

#------
HCO_data_DS1$dpt <- max(HCO_data_DS1$dpt) - HCO_data_DS1$dpt
FeaturePlot(HCO_data_DS1, c("dpt","SOX2","NHLH1","DCX"), ncol=4)

set.seed(12345)
idx <- sample(which(HCO_data_DS1@active.ident %in% c('Dorsal telen. NPC',
                                               'G2M dorsal telen. NPC',
                                               'Dien. and midbrain NPC',
                                               'G2M Dien. and midbrain NPC')),3)
dpt2 <- DPT(dm, tips=idx)
HCO_data_DS1$dpt2 <- rank(dpt2$dpt)
FeaturePlot(HCO_data_DS1, c("dpt","dpt2"), ncol=2)

tips_cand <- sapply(1:100, function(i){ random_root(dm) })
idx_NPC <- which(HCO_data_DS1@active.ident %in% c('Dorsal telen. NPC',
                                            'G2M dorsal telen. NPC',
                                            'Dien. and midbrain NPC',
                                            'G2M Dien. and midbrain NPC'))
tips_cand <- as.numeric(names(which.max(table(tips_cand[tips_cand %in% idx_NPC]))))
dpt3 <- DPT(dm, tips=tips_cand)
HCO_data_DS1$dpt3 <- rank(dpt3$dpt)

FeaturePlot(HCO_data_DS1, c("dpt","dpt2", "dpt3"), ncol=3)

#======================================================================================
# RNA velocity analysis

library(Seurat)
library(Matrix)
# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

mats <- readRDS("data/DS1/mats_dropest.rds")
# Check what barcodes are present in the Seurat object and dropEst matrices
seurat_barcodes <- colnames(HCO_data_DS1)

# Check structure of mats
str(mats)

# Identify intersection of barcodes to prevent subscript errors
intersected_barcodes <- Reduce(intersect, list(
  seurat_barcodes,
  colnames(mats$exon),
  colnames(mats$intron),
  colnames(mats$spanning)
))

# Subset dropEst matrices to only those barcodes in Seurat object
mats <- lapply(mats, function(mat) mat[, intersected_barcodes, drop = FALSE])

mats$exon <- as(as.matrix(mats$exon), "dgCMatrix")
mats$intron <- as(as.matrix(mats$intron), "dgCMatrix")

# Create the loom file:
library(loomR)
cell_attrs <- list(pca = Embeddings(HCO_data_DS1,"pca")[,1:20],
                   umap = Embeddings(HCO_data_DS1,"umap"),
                   celltype = HCO_data_DS1@active.ident)
shared_genes <- intersect(rownames(mats$exon), rownames(mats$intron))
loom <- loomR::create("data/DS1/loom_obj_scvelo.loom",
                      data = mats$exon[shared_genes, , drop = FALSE],
                      layers = list(spliced = mats$exon[shared_genes, , drop = FALSE],
                                    unspliced = mats$intron[shared_genes, , drop = FALSE]),
                      cell.attrs = as.list(HCO_data_DS1@meta.data[intersected_barcodes, , drop = FALSE]))
loom$close_all()

# Create the h5ad file:
HCO_data_DS1 <- subset(HCO_data_DS1, cells = intersected_barcodes) # to match the cells

library(anndata)

shared_genes <- intersect(rownames(mats$exon),rownames(mats$intron))
adata <- AnnData(X = t(mats$exon[shared_genes,]),
                 obs = data.frame(HCO_data_DS1@meta.data, celltype=HCO_data_DS1@active.ident),
                 var = NULL,
                 layers = list(spliced = t(mats$exon[shared_genes,]),
                               unspliced = t(mats$intron[shared_genes,])),
                 obsm = list(pca = Embeddings(HCO_data_DS1,"pca")[,1:20],
                             umap = Embeddings(HCO_data_DS1,"umap"))
)
adata$write_h5ad("data/DS1/anndata_obj_scvelo.h5ad")


library(reticulate)
py_install("scvelo", pip=T)
scvelo <- import("scvelo")
scanpy <- import("scanpy")

# Run the RNA velocity.

adata_DS1 <- scanpy$read_loom("data/DS1/loom_obj_scvelo.loom") # option1: load the loom file
adata_DS1 <- scanpy$read_h5ad("data/DS1/anndata_obj_scvelo.h5ad")
print(adata_DS1$layers$keys()) # to check layer exists

scvelo$pp$filter_and_normalize(adata_DS1,
                               min_shared_counts=as.integer(10),
                               n_top_genes=as.integer(3000))
scvelo$pp$moments(adata_DS1,
                  n_neighbors = as.integer(30),
                  use_rep = "pca")
scvelo$tl$velocity(adata_DS1)
scvelo$tl$velocity_graph(adata_DS1)

# Visualize the velocity estimates.
library(reticulate)
scvelo <- import("scvelo")
plt <- import("matplotlib.pyplot")

# Set display parameters
plt$rcParams["figure.dpi"] <- 200L
plt$rcParams$figure.figsize <- c(12, 8)  # Wider figure
plt$rcParams["font.size"] <- 15L
plt$rcParams["legend.fontsize"] <- 10L
plt$rcParams["legend.title_fontsize"] <- 20L

# Critical legend spacing adjustments
plt$rcParams["legend.labelspacing"] <- 0.5  # More vertical space
plt$rcParams["legend.handletextpad"] <- 0.5  # Symbol-to-text spacing
plt$rcParams["legend.borderpad"] <- 1     # Legend box padding

# Generate plot with proper legend control
scvelo$pl$velocity_embedding_stream(
  adata_DS1,
  basis = "umap",
  color = "celltype",
  linewidth = 0.8,
  density = 3,  # line density
  arrow_size = 2,
  size = 200,
  legend_loc = "on data",  # Moved to right side for more space
  legend_fontsize = 15L,
  fontsize = 30L,
  show = FALSE
)

# Manually adjust legend after creation
current_legend <- plt$gca()$get_legend()
#current_legend$set_bbox_to_anchor(c(1.25, 0.5))  # Move legend further right
#current_legend$set_ncols(1)  # Force single column

# Adjust plot margins
plt$subplots_adjust(
  left = 0.1,
  right = 0.7,  # More space for legend
  top = 0.95,
  bottom = 0.1
)
plt$grid(FALSE)
# Save with high quality
plt$savefig(
  "velocity_plot.png",
  dpi = 200L,
  bbox_inches = "tight",
  pad_inches = 0.4
)
#-------------------------------------------------------------------------------

# velocity pseudotime

# Load individual Seurat objects
HCO_data_DS1 <- readRDS("data/DS1/HCO_dat_obj_all.rds")
HCO_data_DS2 <- readRDS("data/DS2/HCO_dat_obj_all.rds")
HCO_data_DS3 <- readRDS("data/DS3/HCO_dat_obj_all.rds")

scvelo$tl$velocity_pseudotime(adata_DS1)
scvelo$tl$recover_dynamics(adata_DS1)
scvelo$tl$recover_latent_time(adata_DS1)

# UMAP as feature plots for visualization and comparison.
HCO_data_DS1$velocity_pseudotime <- adata_DS1[['obs']]$velocity_pseudotime
HCO_data_DS1$latent_time <- adata_DS1[['obs']]$latent_time

FeaturePlot(HCO_data_DS1,
            c("velocity_pseudotime", "latent_time")) & NoAxes()


#  AnnData and the Seurat object with the velocity-based pseudotime information
adata_DS1$write_h5ad('data/DS1/anndata_DS1_obj_scvelo_PsT.h5ad')
saveRDS(HCO_data_DS1, file='data/DS1/HCO_data_DS1_obj_all_PsT.rds')
