<h2> Nextflow pipeline for Single Cell sequence analysis </h2>
## <h3>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h3> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

## <h3>2). Seurat: R toolkit for single cell genomics. </h3> Detailed information can be found here (https://satijalab.org/seurat)

## Seurat pipeline</h3>
### <ins>Single dataset scRNA-seq analysis</ins>
 

### <ins>Multiple dataset scRNA-seq analysis</ins>
<h3>1). Merge two data sets</h3>
It involves combining the gene expression data of two or more single-cell RNA sequencing (scRNA-seq) datasetsand associated metadata from separate experiments into a unified dataset for downstream analysis. This process is often a prerequisite for data integration, which aims to correct for technical variations and enable the comparison of biological features across datasets. There are different ways to integrate datasets, eg. Seurat, Harmony, Liger, MNN,CSS etc.Seurat and Harmony integration are described below.<br>
   
#### i) _Seurat integration_:
Seurat v5 offers a streamlined approach to data integration, primarily through the IntegrateLayers function, which performs integration 
    in a low-dimensional space.Before integration, the data first has to be split into individual samples (i.e. a separate count matrix for 
    each sample). Next each sample is normalized using SCTransform, and then PCA is performed to reduce the dimensionality of the 
    expression data.<br>
      <img width="700" height="400" alt="Surat_mapped" src="https://github.com/user-attachments/assets/41d005a3-be37-43e1-a33a-cec61bb2b6de" /></br>
#### ii). _Harmony integration_
Harmony, an algorithm for robust, scalable, and flexible multi-dataset integration to meet four key challenges of unsupervised scRNAseq joint embedding: scaling to large datasets, identification of both broad populations and fine-grained subpopulations, flexibility to accommodate complex experimental design, and the power to integrate across modalities. </br>
      <img width="700" height="400" alt="harmony" src="https://github.com/user-attachments/assets/2f169423-b5a9-4ac1-b07f-65987cb170cd" />
<h3> 2). Annotate query datasets using reference data</h3>
The transfer of cell type labels from pre-annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of    publicly available annotated datasets which can be used as reference, as well as the number of computational methods for cell type label transfer are    constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset are needed. <br>   

#### 1a). _Transcriptome similarity on cell cluster level_
####   b). _Transcriptome similarity on cell level_

      <img width="700" height="600" alt="Transcriptome similarity on cell level2" src="https://github.com/user-attachments/assets/75e2b7bc-1a68-477d-bcf5-5dfe845a7bd0" /></br>
&emsp; _2). Seurat-based label transfer_</br>
      <img width="700" height="400" alt="Seurat" src="https://github.com/user-attachments/assets/0d5dcd00-8bdc-4dbb-a6bd-18a3f36136da">
<h3> 3) Advanced analysis for scRNA-seq data</h3>

#### 1). _Cluster connectivity analysis with PAGA_</br>
#### 2). _Pseudotime reconstruction without subseting into an unbranched trajectory_</br>
#### 3). _RNA velocity analysis_
    <img width="600" height="500" alt="Velocity-Plot" src="https://github.com/user-attachments/assets/44d90282-b064-4ac8-9b33-83e2e664475a" /></br>
#### 4). _Trajectory analysis with CellRank_
Trajectory analysis with CellRank in single-cell transcriptomics is a powerful method for studying cell fate decisions and cellular dynamics. It combines single-cell gene expression data with RNA velocity information to reconstruct directed cell state trajectories, revealing how cells transition between different states and ultimately differentiate or reprogram.
    <img width="600" height="550" alt="terminal_states_highres" src="https://github.com/user-attachments/assets/821494ac-122d-446b-b6b4-63d4ddfc0712" /></br>

&ensp; 5). Cell communication analysis</br>

   
<h3>Resources</h3>
Seurat https://github.com/satijalab/seurat</br>
SeuratData https://github.com/satijalab/seurat-data</br>
scVelo Estimating RNA Velocity using Seurat and scVelo https://scvelo.readthedocs.io</br>
Velocity Estimating RNA Velocity using Seurat https://velocyto.org</br>
Harmoney Single cell genomics datasets integration https://github.com/immunogenomics/harmony</br>\

