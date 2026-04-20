library(tidyverse)
HEADER <- c("Chr", "START", "STOP", "WBGENE", "SV_Chr", "SV_START", "SV_STOP", "SVName", "STRAINS", "OVERLAP")
overlapped_svs <- snakemake@input[[1]]
overlapped_svs <- read.table(overlapped_svs, stringsAsFactors = F, header = F)
colnames(overlapped_svs) <- HEADER

genes_of_interest <- snakemake@input[[2]]
genes_of_interest <- read.table(genes_of_interest, stringsAsFactors = F, header = T)
genes_of_interest<- genes_of_interest %>% select(all_of(c("SVName", "Genes", "Wormbase"))) %>%
  separate_longer_delim(c("Genes", "Wormbase"), delim = ",")

colnames(genes_of_interest) <- c("SVName", "Symbol", "WBGENE")

DESIRED_COLS <- c("SVName", "Symbol", "WBGENE", "gene_size", "sv_size", "OVERLAP", "overlap_prop", "STRAINS")
overlapped_svs <- overlapped_svs %>% right_join(genes_of_interest, join_by(WBGENE,SVName)) %>%
  mutate(STOP = STOP -1 ,SV_STOP = SV_STOP -1) %>%
  mutate(gene_size = STOP-START+1, sv_size = SV_STOP-SV_START) %>%
  mutate(OVERLAP = if_else(OVERLAP < gene_size, OVERLAP - 1, OVERLAP)) %>%
  mutate(overlap_prop = OVERLAP/gene_size) %>%
  dplyr::select(all_of(DESIRED_COLS)) %>%
  filter(!is.na(SVName)) %>%
  mutate(overlap_prop = round(overlap_prop, digits = 2)) %>%
  arrange(Symbol)

write_delim(overlapped_svs, snakemake@output[[1]], delim = "\t", col_names = T)
