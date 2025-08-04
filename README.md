<h2> Nextflow pipeline for Single Cell sequence analysis </h2>

<h4>1). Cell Ranger pipeline: Sequencing data preprocessing, including base calling, mapping and read counting.</h4> More information can be found here (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

<h4>2). Seurat: R toolkit for single cell genomics. </h4> More information can be found here (https://satijalab.org/seurat)

a) Single dataset scRNA-seq analysis

b) Multiple dataset scRNA-seq analysis</br>
  Merge two data sets</br>
    1). Seurat integration</br>
    2). Harmony integration</br>
c) Data annotation with reference data</br>
    1a). Transcriptome similarity on cell cluster level</br>
    1b). Transcriptome similarity on cell level</br>
      <img width="600" height="600" alt="umap_transfer_anchor" src="https://github.com/user-attachments/assets/b7eb68f3-f09e-48a8-bab0-6ee8546f24ca" /></br>
    2. Seurat-based label transfer</br>
d) Advanced analysis for scRNA-seq data</br>
    1). Cluster connectivity analysis with PAGA</br>
    2). Pseudotime reconstruction without subseting into an unbranched trajectory</br>
    3). RNA velocity analysis</br>
    <img width="600" height="500" alt="Velocity-Plot" src="https://github.com/user-attachments/assets/44d90282-b064-4ac8-9b33-83e2e664475a" /></br>
    4). Trajectory analysis with CellRank</br>
    5). Cell communication analysis</br>

   
<h3>Resources</h3>
Seurat https://github.com/satijalab/seurat</br>
SeuratData https://github.com/satijalab/seurat-data</br>
scVelo Estimating RNA Velocity using Seurat and scVelo https://scvelo.readthedocs.io</br>
Velocity Estimating RNA Velocity using Seurat https://velocyto.org</br>
