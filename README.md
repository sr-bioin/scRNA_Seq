<h2> Nextflow pipeline for Single Cell sequence analysis </h2>

<h4>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h4> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

<h3>2). Seurat: R toolkit for single cell genomics. </h3> More information can be found here (https://satijalab.org/seurat)

<h3>Seurat pipeline</h3>
<h4>Single dataset scRNA-seq analysis</h4>


<h4>Multiple dataset scRNA-seq analysis</h4>


   &ensp; 1). Merge two data sets</br>
      &emsp; 1). Seurat integration</br>
      &emsp; 2). Harmony integration</br>
      
   &ensp; 2).Data annotation with reference data</br>
          &emsp; 1a). Transcriptome similarity on cell cluster level</br>
          &emsp; 1b). Transcriptome similarity on cell level</br>
      <img width="600" height="600" alt="umap_transfer_anchor" src="https://github.com/user-attachments/assets/b7eb68f3-f09e-48a8-bab0-6ee8546f24ca" /></br>
          &emsp; 2). Seurat-based label transfer</br>


3) Advanced analysis for scRNA-seq data</br>
    &ensp; 1). Cluster connectivity analysis with PAGA</br>
    &ensp; 2). Pseudotime reconstruction without subseting into an unbranched trajectory</br>
   &ensp;  3). RNA velocity analysis</br>
    <img width="600" height="500" alt="Velocity-Plot" src="https://github.com/user-attachments/assets/44d90282-b064-4ac8-9b33-83e2e664475a" /></br>
    &ensp; 4). Trajectory analysis with CellRank</br>
        - Trajectory analysis with CellRank in single-cell transcriptomics is a powerful method for studying cell fate decisions and cellular dynamics. It               - combines single-cell gene expression data with RNA velocity information to reconstruct directed cell state trajectories, revealing how cells                   - transition between different states and ultimately differentiate or reprogram."
    <img width="600" height="550" alt="terminal_states_highres" src="https://github.com/user-attachments/assets/821494ac-122d-446b-b6b4-63d4ddfc0712" /></br>

   &ensp;  5). Cell communication analysis</br>

   
<h3>Resources</h3>
Seurat https://github.com/satijalab/seurat</br>
SeuratData https://github.com/satijalab/seurat-data</br>
scVelo Estimating RNA Velocity using Seurat and scVelo https://scvelo.readthedocs.io</br>
Velocity Estimating RNA Velocity using Seurat https://velocyto.org</br>
