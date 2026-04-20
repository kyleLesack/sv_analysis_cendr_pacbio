library(tidyverse)

columns <- c("Chr", "Start", "Stop", "Region", "Count")

dels <- read_delim(snakemake@input[["deletion_distribution"]], col_names = columns)
dels$SVType <- "DEL"
dels <- dels %>% select(all_of(c("Chr", "SVType", "Region", "Count"))) %>%
  pivot_wider(names_from = Region, values_from = Count) %>%
  mutate(Total = LARM + CENTER + RARM)

dups <- read_delim(snakemake@input[["duplication_distribution"]], col_names = columns)
dups$SVType <- "DUP"
dups <- dups %>% select(all_of(c("Chr", "SVType", "Region", "Count"))) %>%
  pivot_wider(names_from = Region, values_from = Count) %>%
  mutate(Total = LARM + CENTER + RARM)

invs <- read_delim(snakemake@input[["inversion_distribution"]], col_names = columns)
invs$SVType <- "INV"
invs <- invs %>% select(all_of(c("Chr", "SVType", "Region", "Count"))) %>%
  pivot_wider(names_from = Region, values_from = Count) %>%
  mutate(Total = LARM + CENTER + RARM)

all_svs <- rbind(dels, dups, invs)
all_svs <- all_svs %>% arrange(Chr, SVType) %>% mutate(PropArms = (LARM + RARM)/Total) %>%
  mutate(PropArms = round(PropArms, digits = 2))
write_delim(all_svs, snakemake@output[[1]], delim = ",", col_names = T)

                    