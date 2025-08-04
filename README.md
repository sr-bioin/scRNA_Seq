<h2> Nextflow pipeline for Single Cell sequence analysis </h2>

<h3>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h3> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

<h3>2). Seurat: R toolkit for single cell genomics. </h3> More information can be found here (https://satijalab.org/seurat)

## Seurat pipeline</h3>
### Single dataset scRNA-seq analysis


#### Multiple dataset scRNA-seq analysis


   &ensp; 1). Merge two data sets</br>
      &emsp; 1). Seurat integration</br>
      &emsp; 2). Harmony integration</br>
      <img width="600" height="300" alt="harmony" src="https://github.com/user-attachments/assets/2f169423-b5a9-4ac1-b07f-65987cb170cd" />

      
   &ensp; 2).Data annotation with reference data</br>

   The transfer of cell type labels from pre-annotated (reference) to newly collected data is an important task in single-cell data analysis. As the number of publicly available annotated datasets which can be used as reference, as well as the number of computational methods for cell type label transfer are constantly growing, rationals to understand and decide which reference design and which method to use for a particular query dataset are needed.
   
          &emsp; 1a). Transcriptome similarity on cell cluster level</br>
          &emsp; 1b). Transcriptome similarity on cell level</br>
      <img width="700" height="600" alt="Transcriptome similarity on cell level2" src="https://github.com/user-attachments/assets/75e2b7bc-1a68-477d-bcf5-5dfe845a7bd0" />
</br>
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


