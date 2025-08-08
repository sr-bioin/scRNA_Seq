<h2> Nextflow pipeline for Single Cell sequence analysis </h2>
## <h3>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h3> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

## <h3>2). Seurat: R toolkit for single cell genomics. </h3> Detailed information can be found here (https://satijalab.org/seurat)

## Seurat pipeline</h3>
### <ins>Single dataset scRNA-seq analysis</ins>
1. Data Loading & Preprocessing</br>
&emsp;a) Create a Seurat object</br>
&emsp;b).Quality control</br>
2. Normalization
3. Feature selection
4. Data Scaling & Regressing Out Confounders
5. Dimensionality reduction using principal component analysis (PCA)
6. Non-linear dimension reduction (UMAP/t-SNE)
7. Cluster the cells
8. Annotate cell clusters
9. Pseudotemporal cell ordering


### <ins>Multiple dataset scRNA-seq analysis</ins>
<h3>1). Merge two data sets</h3>
It involves combining the gene expression data of two or more single-cell RNA sequencing (scRNA-seq) datasetsand associated metadata from separate experiments into a unified dataset for downstream analysis. This process is often a prerequisite for data integration, which aims to correct for technical variations and enable the comparison of biological features across datasets. There are different ways to integrate datasets, eg. Seurat, Harmony, Liger, MNN,CSS etc.Seurat and Harmony integration are described below.<br>
   
#### i) _Seurat integration_:
Seurat v5 offers a streamlined approach to data integration, primarily through the IntegrateLayers function, which performs integration 
    in a low-dimensional space.Before integration, the data first has to be split into individual samples (i.e. a separate count matrix for each sample). Next each sample is normalized using SCTransform, and then PCA is performed to reduce the dimensionality of the 
    expression data.<br>
      <img width="700" height="400" alt="Surat_mapped" src="https://github.com/user-attachments/assets/41d005a3-be37-43e1-a33a-cec61bb2b6de" /></br>
#### ii). _Harmony integration_
It is an algorithm for robust, scalable, and flexible multi-dataset integration to meet four key challenges of unsupervised scRNAseq joint embedding: scaling to large datasets, identification of both broad populations and fine-grained subpopulations, flexibility to accommodate complex experimental design, and the power to integrate across modalities. </br>
      <img width="700" height="400" alt="harmony" src="https://github.com/user-attachments/assets/2f169423-b5a9-4ac1-b07f-65987cb170cd" />
<h3> 2). Annotate query datasets using reference data</h3>
The transfer of cell type labels from pre-annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of    publicly available annotated datasets which can be used as reference, as well as the number of computational methods for cell type label transfer are    constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset are needed. <br>   

#### i). _Transcriptome similarity on cell cluster level_
#### ii). _Transcriptome similarity on cell level_


&emsp; _2). Seurat-based label transfer_</br>
<img width="600" height="900" alt="Seurat-based_label_transfer_SCT]" src="https://github.com/user-attachments/assets/a77c36c2-59bc-4b92-87dc-feb5e3354a9d"/></br>


<h3> 3) Advanced analysis for scRNA-seq data</h3>

#### i). _Cluster connectivity analysis with PAGA_
Partition-based graph abstraction (PAGA) provides an interpretable graph-like map of the arising data manifold, based on estimating connectivity of manifold partitions (https://github.com/theislab/paga). PAGA maps preserve the global topology of data, allow analyzing data at different resolutions, and result in much higher computational efficiency of the typical exploratory data analysis workflow. <br>
<img width="600" height="600" alt="DS1_paga_SCT" src="https://github.com/user-attachments/assets/50485134-f44f-4fe2-831c-3bf843bcc5e8" />

#### ii). _Pseudotime reconstruction without subseting into an unbranched trajectory_</br>
#### iii). _RNA velocity analysis_
RNA velocity is a high-dimensional vector that predicts the future state of individual cells on a timescale of hours. It is the time derivative of gene expression state (ds/dt with s representing the high-dimensional expression state and t time) and is widely used to infer temporal dynamics in single-cell gene expression data.
<img width="600" height="450" alt="velocity_plot_SCT" src="https://github.com/user-attachments/assets/75690caa-b327-40ac-942c-c2a745da356f" />

#### iv). _Trajectory analysis with CellRank_
Trajectory analysis with CellRank in single-cell transcriptomics is a powerful method for studying cell fate decisions and cellular dynamics. It combines single-cell gene expression data with RNA velocity information to reconstruct directed cell state trajectories, revealing how cells transition between different states and ultimately differentiate or reprogram.
<img width="600" height="550" alt="terminal_states_highres" src="https://github.com/user-attachments/assets/821494ac-122d-446b-b6b4-63d4ddfc0712" /></br>
#### v). _Cell communication analysis_

### <ins>Multiple dataset scRNA-seq analysis</ins> 
<img width="600" height="333" alt="DS1DS2DS3_integrationN" src="https://github.com/user-attachments/assets/8ae5fda2-ca09-4a8d-acb5-7303bc83f8d2" />


<h3>Resources</h3>
Seurat https://github.com/satijalab/seurat</br>
SeuratData https://github.com/satijalab/seurat-data</br>
scVelo Estimating RNA Velocity using Seurat and scVelo https://scvelo.readthedocs.io</br>
Velocity Estimating RNA Velocity using Seurat https://velocyto.org</br>
Harmoney Single cell genomics datasets integration https://github.com/immunogenomics/harmony</br>\

