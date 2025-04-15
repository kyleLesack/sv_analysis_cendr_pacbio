# Software Requirements

Note: The [yaml](./yaml/) directory includes conda yaml environment files. Most tools are available as Conda packages, which I recommend highly for replicating the pipeline. The required tools and the versions we used are listed below:

## Sequence Alignment and Structural Variant Calling

The analysis of *Caenorhabditis elegans* structural variants (SVs) was performed using the following tools:

* NGMLR (v0.2.7)
* Picard (v2.27.5)
* Sniffles (v.2.0.7)
* Jasmine (v1.1.5) 
* Iris (v1.0.4)
* Python (v3.10)

## Structural Variant Annotation

SV annotation was performed using the following tools:

* Variant Effect Predictor (v.109.3)
	* To narrow down the SVs to those that were most likely to have significant phenotypic consequences, the VEP filter parameter was set to "IMPACT is HIGH".
* easyGSEA tool from eVITTA (v1.3.1) 
	* Overrepresented genes were annotated using the annotated gene sets from WormCat (V2)


## Snakemake pipeline

A Snakemake pipeline was created to manage the SV calling part of the analysis. Most tools are available as Conda packages (specified in the yaml directory). Python, and Linux shell scripts were used to process the analysis results. This produces a finalized merged dataset with high confidence SV calls. These VCF files should be available in the 3_jasmine/ngmlr/sniffles/8_final/sv_types/ directory if the pipeline ran correctly.


# Input data

*C. elegans PacBio* sequencing data from the [*Caenorhabditis* Natural Diversity Resource](https://caendr.org/) project were used and are available from the Sequence Read Archive: [PacBio sequencing data](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=12908562) 

The Snakemake pipeline expects the FASTQ files to be in the following locations:

0_input/fastq/XZ1516_all_reads.fastq

0_input/fastq/DL238_all_reads.fastq

0_input/fastq/JU1400_all_reads.fastq

0_input/fastq/JU2600_all_reads.fastq

0_input/fastq/MY2147_all_reads.fastq

0_input/fastq/ECA396_all_reads.fastq

0_input/fastq/EG4725_all_reads.fastq

0_input/fastq/JU2526_all_reads.fastq

0_input/fastq/JU310_all_reads.fastq

0_input/fastq/NIC526_all_reads.fastq

0_input/fastq/QX1794_all_reads.fastq

0_input/fastq/MY2693_all_reads.fastq

0_input/fastq/NIC2_all_reads.fastq

0_input/fastq/ECA36_all_reads.fastq


## The Snakemake pipeline expects the *C. elegans* reference genome to be in the following location:

0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa

