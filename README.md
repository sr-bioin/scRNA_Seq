<h2> Nextflow pipeline for Single Cell sequence analysis </h2>

<h3>Quality control fo the raw sequencing data generated from the sequencing instrument to remove low-quality reads and adapter sequences.</h3>

1). FastQC is used to assess the quality distribution and other quality metrics of sequencing reads. It helps detect potential issues during the sequencing process, such as decreasing sequencing quality and adapter contamination. For more information, please visit the FastQC website.

2). MultiQC is a tool to create a single report with interactive plots for multiple bioinformatics analyses across many samples. MultiQC reports can describe multiple analysis steps and large numbers of samples within a single plot, and multiple analysis tools making it ideal for routine fast quality control.For more information, please visit the MultiQC website.

3). Trimmomatic: Used to remove adapter sequences, low-quality bases, and low-quality reads from sequencing reads. You can find more information and usage instructions on the Trimmomatic GitHub.


