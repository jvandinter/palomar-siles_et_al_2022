# README #

This Repository is a code deposition of the performed Ribo-seq analyses and raw figures in the FUr readthrough manuscript (link TBA). Code to generate the figures is found under **figures**. *SessionInfo.RMD* holds all required packages for the R code. Code to process Ribo-seq libraries to trimmed reads is located in the **pipeline** folder. Please note that our Ribo-seq pipeline is still in active development and that this is merely a historical snapshot to replicate our analyses. If your are interested in our current pipeline, please contact us.

### Subfolders ###

- Data
- Figures
- Pipeline

### Overview of the Ribo-seq pipeline ###

Required files are in the Ribo-seq config located in the **data** folder. Our local computer cluster uses the SLURM workload manager. In addition, we use the LMOD system to load specific versions of packages.

1. Quality control and 3' adapter trimming with **TrimGalore**
2. Remove unwanted RNAs with **bowtie2**
3. Map ribosome fragments with **STAR**
4. Run **RiboseQC** for Ribo-seq QC
5. Use **FeatureCount** from the subread package to quantify CDS counts

Before running the pipeline, make sure you have the following files:

- Filled in config file (**pipeline/mrp.config**)
- BOWTIE2 index of contaminant sequences that should be removed (rRNA, tRNA, snRNA, snoRNA, mtDNA)
- STAR index of your reference genome of choice (we used ensembl v102)
- Genome annotation file for your reference genome (we used ensembl v102)
- Create annotation files for RiboseQC (see https://github.com/ohlerlab/RiboseQC)

### Contributors ###

Contact us for more information regarding our analyses:

* j.t.vandinter-3@prinsesmaximacentrum.nl
* d.a.hofman-3@prinsesmaximacentrum.nl
* s.vanheesch@prinsesmaximacentrum.nl
