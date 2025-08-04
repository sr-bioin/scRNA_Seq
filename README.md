<h2> Nextflow pipeline for Single Cell sequence analysis </h2>
## <h3>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h3> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

## <h3>2). Seurat: R toolkit for single cell genomics. </h3> Detailed information can be found here (https://satijalab.org/seurat)

## Seurat pipeline</h3>
### Single dataset scRNA-seq analysis
 

### Multiple dataset scRNA-seq analysis
#### 1). Merge two data sets
<h4>1). Merge two data sets</h3>
   It involves combining the gene expression data of two or more single-cell RNA sequencing (scRNA-seq) datasetsand associated metadata from separate experiments into a unified dataset for downstream analysis. This process is often a prerequisite for data integration, which aims to correct for technical variations and enable the comparison of biological features across datasets. There are different ways to integrate datasets, eg. Seurat, Harmony, Liger, MNN,CSS etc.Seurat and Harmony integration are described below.<br>
##### 1). Seurat integration
      <br>&emsp; 1). Seurat integration<br>
        Seurat v5 offers a streamlined approach to data integration, primarily through the IntegrateLayers function, which performs integration in a
    low-dimensional space.Before integration, the data first has to be split into individual samples (i.e. a separate count matrix for each sample).
    Next each sample is normalized using SCTransform, and then PCA is performed to reduce the dimensionality of the expression data.<br>
      <img width="700" height="400" alt="Surat_mapped" src="https://github.com/user-attachments/assets/41d005a3-be37-43e1-a33a-cec61bb2b6de" /></br>
     <br> &emsp; 2). Harmony integration</br>
    <br> Harmony, an algorithm for robust, scalable, and flexible multi-dataset integration to meet four key challenges of unsupervised scRNAseq joint embedding: scaling to large datasets, identification of both broad populations and fine-grained subpopulations, flexibility to accommodate complex experimental design, and the power to integrate across modalities. </br>
      <img width="700" height="400" alt="harmony" src="https://github.com/user-attachments/assets/2f169423-b5a9-4ac1-b07f-65987cb170cd" />
    <h3><br>2). Annotate query datasets using reference data </h3>
The transfer of cell type labels from pre-annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of    publicly available annotated datasets which can be used as reference, as well as the number of computational methods for cell type label transfer are    constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset are needed.    
          <br>&emsp; 1a). Transcriptome similarity on cell cluster level</br>
          &emsp; 1b). Transcriptome similarity on cell level</br>
      <img width="700" height="600" alt="Transcriptome similarity on cell level2" src="https://github.com/user-attachments/assets/75e2b7bc-1a68-477d-bcf5-5dfe845a7bd0" /></br>
          &emsp; 2). Seurat-based label transfer</br>
      <img width="700" height="400" alt="Seurat" src="https://github.com/user-attachments/assets/0d5dcd00-8bdc-4dbb-a6bd-18a3f36136da">

 



3) Advanced analysis for scRNA-seq data</br>
    &ensp; 1). Cluster connectivity analysis with PAGA</br>
    &ensp; 2). Pseudotime reconstruction without subseting into an unbranched trajectory</br>
   &ensp;  3). RNA velocity analysis</br>
    <img width="600" height="500" alt="Velocity-Plot" src="https://github.com/user-attachments/assets/44d90282-b064-4ac8-9b33-83e2e664475a" /></br>
    &ensp; 4). Trajectory analysis with CellRank

Trajectory analysis with CellRank in single-cell transcriptomics is a powerful method for studying cell fate decisions and cellular dynamics. It combines single-cell gene expression data with RNA velocity information to reconstruct directed cell state trajectories, revealing how cells transition between different states and ultimately differentiate or reprogram.
    <img width="600" height="550" alt="terminal_states_highres" src="https://github.com/user-attachments/assets/821494ac-122d-446b-b6b4-63d4ddfc0712" /></br>

   &ensp;  5). Cell communication analysis</br>

   
<h3>Resources</h3>
Seurat https://github.com/satijalab/seurat</br>
SeuratData https://github.com/satijalab/seurat-data</br>
scVelo Estimating RNA Velocity using Seurat and scVelo https://scvelo.readthedocs.io</br>
Velocity Estimating RNA Velocity using Seurat https://velocyto.org</br>
Harmoney Single cell genomics datasets integration https://github.com/immunogenomics/harmony</br>\

1. Hello\
    This sentence is indented to match the heading.

2. Heading\
    This sentence is indented to match the heading.
    This sentence is indented to match the heading.





# Header1
## Header 2
### Multiple dataset scRNA-seq analysis

#### 1). Merge two data sets

It involves combining the gene expression data of two or more single-separate experiments into a unified dataset for downstream analysis. Correct for technical variations and enable the comparison of biological conditions. Seurat, Harmony, LIGER, MNN, CSS, etc. Seurat and Harmony integration approaches  

**a). Seurat integration**\
    Seurat v5 offers a streamlined approach to data integration, primarily in low-dimensional space.  
    **Seurat v5** is now properly aligned under the "Seurat integration" heading.
