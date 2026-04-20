library(tidyverse)
header <- c("SVName", "Size")

# Deletions

del_size_hyper <- read_tsv(snakemake@input[["deletion_sizes_hyper"]],header)
del_size_hyper$Region <- "Hyperdivergent"
del_size_non_hyper <- read_tsv(snakemake@input[["deletion_sizes_non_hyper"]],header)
del_size_non_hyper$Region <- "Non-hyperdivergent"

del_mean_hyper <- round(mean(del_size_hyper$Size),0)
del_sd_hyper <- round(sd(del_size_hyper$Size),0)
del_mean_sd_hyper <- paste(del_mean_hyper, del_sd_hyper, sep = "±")
del_mean_non_hyper <- round(mean(del_size_non_hyper$Size),0)
del_sd_non_hyper <- round(sd(del_size_non_hyper$Size),0)
del_mean_sd_non_hyper <- paste(del_mean_non_hyper, del_sd_non_hyper, sep = "±")
del_row <- c("DEL", del_mean_sd_hyper, del_mean_sd_non_hyper)

# Duplications

dup_size_hyper <- read_tsv(snakemake@input[["duplication_sizes_hyper"]],header)
dup_size_hyper$Region <- "Hyperdivergent"
dup_size_non_hyper <- read_tsv(snakemake@input[["duplication_sizes_non_hyper"]],header)
dup_size_non_hyper$Region <- "Non-hyperdivergent"

dup_mean_hyper <- round(mean(dup_size_hyper$Size),0)
dup_sd_hyper <- round(sd(dup_size_hyper$Size),0)
dup_mean_sd_hyper <- paste(dup_mean_hyper, dup_sd_hyper, sep = "±")
dup_mean_non_hyper <- round(mean(dup_size_non_hyper$Size),0)
dup_sd_non_hyper <- round(sd(dup_size_non_hyper$Size),0)
dup_mean_sd_non_hyper <- paste(dup_mean_non_hyper, dup_sd_non_hyper, sep = "±")
dup_row <- c("DUP", dup_mean_sd_hyper, dup_mean_sd_non_hyper)

# Inversions

inv_size_hyper <- read_tsv(snakemake@input[["inversion_sizes_hyper"]],header)
inv_size_hyper$Region <- "Hyperdivergent"
inv_size_non_hyper <- read_tsv(snakemake@input[["inversion_sizes_non_hyper"]],header)
inv_size_non_hyper$Region <- "Non-hyperdivergent"

inv_mean_hyper <- round(mean(inv_size_hyper$Size),0)
inv_sd_hyper <- round(sd(inv_size_hyper$Size),0)
inv_mean_sd_hyper <- paste(inv_mean_hyper, inv_sd_hyper, sep = "±")
inv_mean_non_hyper <- round(mean(inv_size_non_hyper$Size),0)
inv_sd_non_hyper <- round(sd(inv_size_non_hyper$Size),0)
inv_mean_sd_non_hyper <- paste(inv_mean_non_hyper, inv_sd_non_hyper, sep = "±")
inv_row <- c("INV", inv_mean_sd_hyper, inv_mean_sd_non_hyper)

# Combined

all_svs <- data.frame(rbind(del_row,dup_row, inv_row))
colnames(all_svs) <- c("SVType", "Hyperdivergent (Mean ± SD)", "Non-Hyperdivergent (Mean ± SD)")
write_delim(all_svs, snakemake@output[[1]], delim = ",", col_names = T)
