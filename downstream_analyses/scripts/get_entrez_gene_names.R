library(tidyverse)
library(org.Ce.eg.db)
library(vcfR)

svs <- read.table(snakemake@input[[1]], stringsAsFactors = F, header = T)
svs <- svs %>% dplyr::rename(WORMBASE = Genes) %>% separate_longer_delim(WORMBASE, delim = ",") 

svs <- svs %>% mutate(genesymbol = mapIds(org.Ce.eg.db,
                                         keys = svs$WORMBASE,
                                         column = "SYMBOL",
                                         keytype = "WORMBASE"))

svs <- svs %>% group_by(SVName, Strains) %>%
  dplyr::summarise(Genes = paste(genesymbol, collapse = ","),
                   Wormbase = paste(WORMBASE, collapse = ","))

vcf <- read.vcfR(snakemake@input[[2]], verbose = FALSE )
vcf <- vcfR2tidy(vcf, info_only = T, info_fields = "SVLEN")
vcf_fields <- c("ID", "SVLEN")
svs <- vcf$fix %>% dplyr::select(all_of(vcf_fields)) %>% 
  dplyr::rename(SVName = ID) %>% 
  mutate(SVLEN = abs(SVLEN)) %>%
  right_join(svs, by = "SVName")

write_delim(svs, snakemake@output[[1]], delim = "\t", col_names = T)
