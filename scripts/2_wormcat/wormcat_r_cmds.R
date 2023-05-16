library(wormcat)
setwd("/bulk/worm_lab/mrkyle/sv_analysis_cendr_pacbio")
worm_cat_fun("VEP_GENES_FILE" , title = "VARIANT_TYPE Gene Set", output_dir = "OUTDIR", rm_dir = FALSE, annotation_file = "whole_genome_v2_nov-11-2021.csv",	input_type = "Wormbase.ID")
