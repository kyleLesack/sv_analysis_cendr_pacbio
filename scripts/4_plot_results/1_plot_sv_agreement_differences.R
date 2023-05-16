library(ggplot2)
library(tidyverse)
library(patchwork)

## VCF Level Agreement and Differences

# Open input data
pbsv_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

## Convert to data frames

### pbsv
pbsv_pbmm2.df <- data.frame(pbsv_pbmm2)
pbsv_pbmm2.df['Caller'] <- c("pbsv")
pbsv_pbmm2.df <-pbsv_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller, Intersection, Non.intersecting))
pbsv_pbmm2.df <- pbsv_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
pbsv_pbmm2.df['Aligner'] <- c("pbmm2")

### Sniffles
sniffles_pbmm2.df <- data.frame(sniffles_pbmm2)
sniffles_pbmm2.df['Caller'] <- c("Sniffles")
sniffles_pbmm2.df <-sniffles_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_pbmm2.df <- sniffles_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_pbmm2.df['Aligner'] <- c("pbmm2")

sniffles_minimap2.df <- data.frame(sniffles_minimap2)
sniffles_minimap2.df['Caller'] <- c("Sniffles")
sniffles_minimap2.df <-sniffles_minimap2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_minimap2.df <- sniffles_minimap2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_minimap2.df['Aligner'] <- c("Minimap2")

sniffles_ngmlr.df <- data.frame(sniffles_ngmlr)
sniffles_ngmlr.df['Caller'] <- c("Sniffles")
sniffles_ngmlr.df <-sniffles_ngmlr.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_ngmlr.df <- sniffles_ngmlr.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_ngmlr.df['Aligner'] <- c("NGMLR")

## SVIM
svim_pbmm2.df <- data.frame(svim_pbmm2)
svim_pbmm2.df['Caller'] <- c("SVIM")
svim_pbmm2.df <-svim_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_pbmm2.df <- svim_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_pbmm2.df['Aligner'] <- c("pbmm2")

svim_minimap2.df <- data.frame(svim_minimap2)
svim_minimap2.df['Caller'] <- c("SVIM")
svim_minimap2.df <-svim_minimap2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_minimap2.df <- svim_minimap2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_minimap2.df['Aligner'] <- c("Minimap2")

svim_ngmlr.df <- data.frame(svim_ngmlr)
svim_ngmlr.df['Caller'] <- c("SVIM")
svim_ngmlr.df <-svim_ngmlr.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_ngmlr.df <- svim_ngmlr.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_ngmlr.df['Aligner'] <- c("NGMLR")

## Reorder and rename labels
intersection.df <- rbind(pbsv_pbmm2.df, sniffles_minimap2.df, sniffles_ngmlr.df, sniffles_pbmm2.df,svim_minimap2.df, svim_ngmlr.df, svim_pbmm2.df)
col_order <- c("Caller", "Aligner", "Category", "Count")
intersection_reordered.df <- intersection.df[, col_order]
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Intersection", "Intersecting calls")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Non.intersecting", "Different calls")))
intersection_reordered.df$Count = as.numeric(intersection_reordered.df$Count)
intersection_reordered.df$Count <- round(intersection_reordered.df$Count, 0)

p_grouped_by_caller <- ggplot(data=intersection_reordered.df, aes(x=Aligner, y=Count, fill=Category)) + 
	geom_col(position=position_dodge()) + 
	facet_grid(~ Caller, scales = "free", space='free') 

p_grouped_by_caller <- p_grouped_by_caller + xlab("SV Calling and Alignment Method") + 
	scale_y_continuous(labels = scales::comma)  

## Add counts to columns
p_grouped_by_caller_counts <- p_grouped_by_caller + 
	geom_text(aes(x=Aligner, y = Count, label = Count), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

dir.create(file.path("5_plots","full_depth"), recursive = TRUE)
ggsave("1_sv_agreement_grouped_by_caller_counts_top_labels.png", plot = p_grouped_by_caller_counts, device = "png", path = "5_plots/full_depth", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Coordinate Level Agreement and Differences

## Open input data
pbsv_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
sniffles_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_minimap2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_ngmlr <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)
svim_pbmm2 <- read.delim(file = '4_results/full_depth/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_coords_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

## Convert to data frames

### pbsv
pbsv_pbmm2.df <- data.frame(pbsv_pbmm2)
pbsv_pbmm2.df['Caller'] <- c("pbsv")
pbsv_pbmm2.df <-pbsv_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller, Intersection, Non.intersecting))
pbsv_pbmm2.df <- pbsv_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
pbsv_pbmm2.df['Aligner'] <- c("pbmm2")

### sniffles
sniffles_pbmm2.df <- data.frame(sniffles_pbmm2)
sniffles_pbmm2.df['Caller'] <- c("Sniffles")
sniffles_pbmm2.df <-sniffles_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_pbmm2.df <- sniffles_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_pbmm2.df['Aligner'] <- c("pbmm2")

sniffles_minimap2.df <- data.frame(sniffles_minimap2)
sniffles_minimap2.df['Caller'] <- c("Sniffles")
sniffles_minimap2.df <-sniffles_minimap2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_minimap2.df <- sniffles_minimap2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_minimap2.df['Aligner'] <- c("Minimap2")

sniffles_ngmlr.df <- data.frame(sniffles_ngmlr)
sniffles_ngmlr.df['Caller'] <- c("Sniffles")
sniffles_ngmlr.df <-sniffles_ngmlr.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
sniffles_ngmlr.df <- sniffles_ngmlr.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
sniffles_ngmlr.df['Aligner'] <- c("NGMLR")

## Plot results
svim_pbmm2.df <- data.frame(svim_pbmm2)
svim_pbmm2.df['Caller'] <- c("SVIM")
svim_pbmm2.df <-svim_pbmm2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_pbmm2.df <- svim_pbmm2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_pbmm2.df['Aligner'] <- c("pbmm2")

svim_minimap2.df <- data.frame(svim_minimap2)
svim_minimap2.df['Caller'] <- c("SVIM")
svim_minimap2.df <-svim_minimap2.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_minimap2.df <- svim_minimap2.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_minimap2.df['Aligner'] <- c("Minimap2")

svim_ngmlr.df <- data.frame(svim_ngmlr)
svim_ngmlr.df['Caller'] <- c("SVIM")
svim_ngmlr.df <-svim_ngmlr.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Caller,Intersection, Non.intersecting))
svim_ngmlr.df <- svim_ngmlr.df %>% pivot_longer(!Caller, names_to = "Category", values_to = "Count")
svim_ngmlr.df['Aligner'] <- c("NGMLR")

## Reorder and rename labels
intersection.df <- rbind(pbsv_pbmm2.df, sniffles_minimap2.df, sniffles_ngmlr.df, sniffles_pbmm2.df,svim_minimap2.df, svim_ngmlr.df, svim_pbmm2.df)
col_order <- c("Caller", "Aligner", "Category", "Count")
intersection_reordered.df <- intersection.df[, col_order]
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Intersection", "Intersecting calls")))
intersection_reordered.df <- intersection_reordered.df %>% mutate_all(funs(str_replace(., "Non.intersecting", "Different calls")))
intersection_reordered.df$Count = as.numeric(intersection_reordered.df$Count)
intersection_reordered.df$Count <- round(intersection_reordered.df$Count, 0)

p_grouped_by_caller <- ggplot(data=intersection_reordered.df, aes(x=Aligner, y=Count, fill=Category)) + 
	geom_col(position=position_dodge()) + 
	facet_grid(~ Caller, scales = "free", space='free') 

p_grouped_by_caller <- p_grouped_by_caller + xlab("SV Calling and Alignment Method") + 
	scale_y_continuous(labels = scales::comma)

## Add counts to columns
p_grouped_by_caller_counts <- p_grouped_by_caller + 
	geom_text(aes(x=Aligner, y = Count, label = Count), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

dir.create(file.path("5_plots","full_depth","coords"), recursive = TRUE)
ggsave("1_sv_agreement_coords_grouped_by_caller_counts_top_labels.png", plot = p_grouped_by_caller_counts, device = "png", path = "5_plots/full_depth/coords", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

# Subsampled VCF Level Agreement and Disagreement

pbsv_pbmm2_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
pbsv_pbmm2_10X.df$Depth <- c("10X")
pbsv_pbmm2_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
pbsv_pbmm2_20X.df$Depth <- c("20X")
pbsv_pbmm2_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
pbsv_pbmm2_40X.df$Depth <- c("40X")
pbsv_pbmm2_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/pbsv/pbsv-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
pbsv_pbmm2_60X.df$Depth <- c("60X")

## pbsv

pbsv_pbmm2_subsampled.df <- rbind(pbsv_pbmm2_10X.df,pbsv_pbmm2_20X.df,pbsv_pbmm2_40X.df,pbsv_pbmm2_60X.df)
pbsv_pbmm2_subsampled.df['Aligner'] <- c("pbmm2")
pbsv_pbmm2_subsampled.df <-pbsv_pbmm2_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

p_pbsv <- ggplot(pbsv_pbmm2_subsampled.df, aes(x=Depth, y = Non.intersecting.proportion)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 1)
p_pbsv <- p_pbsv + xlab("pbsv") + ylab("Proportion Different SV Calls")

## Sniffles

sniffles_minimap2_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_minimap2_10X.df$Depth <- c("10X")
sniffles_minimap2_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_minimap2_20X.df$Depth <- c("20X")
sniffles_minimap2_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_minimap2_40X.df$Depth <- c("40X")
sniffles_minimap2_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/sniffles-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_minimap2_60X.df$Depth <- c("60X")
sniffles_minimap2_subsampled.df <- rbind(sniffles_minimap2_10X.df, sniffles_minimap2_20X.df, sniffles_minimap2_40X.df, sniffles_minimap2_60X.df)
sniffles_minimap2_subsampled.df$Aligner <- c("minimap2")
sniffles_minimap2_subsampled.df <- sniffles_minimap2_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

sniffles_ngmlr_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_ngmlr_10X.df$Depth <- c("10X")
sniffles_ngmlr_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_ngmlr_20X.df$Depth <- c("20X")
sniffles_ngmlr_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_ngmlr_40X.df$Depth <- c("40X")
sniffles_ngmlr_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/sniffles-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_ngmlr_60X.df$Depth <- c("60X")
sniffles_ngmlr_subsampled.df <- rbind(sniffles_ngmlr_10X.df, sniffles_ngmlr_20X.df, sniffles_ngmlr_40X.df, sniffles_ngmlr_60X.df)
sniffles_ngmlr_subsampled.df$Aligner <- c("ngmlr")
sniffles_ngmlr_subsampled.df <- sniffles_ngmlr_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

sniffles_pbmm2_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_pbmm2_10X.df$Depth <- c("10X")
sniffles_pbmm2_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_pbmm2_20X.df$Depth <- c("20X")
sniffles_pbmm2_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_pbmm2_40X.df$Depth <- c("40X")
sniffles_pbmm2_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/sniffles/sniffles-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
sniffles_pbmm2_60X.df$Depth <- c("60X")
sniffles_pbmm2_subsampled.df <- rbind(sniffles_pbmm2_10X.df, sniffles_pbmm2_20X.df, sniffles_pbmm2_40X.df, sniffles_pbmm2_60X.df)
sniffles_pbmm2_subsampled.df$Aligner <- c("pbmm2")
sniffles_pbmm2_subsampled.df <- sniffles_pbmm2_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

sniffles_all_aligners.df <- rbind(sniffles_minimap2_subsampled.df, sniffles_pbmm2_subsampled.df, sniffles_ngmlr_subsampled.df)
p_sniffles <- ggplot(sniffles_all_aligners.df, aes(x=Depth, y = Non.intersecting.proportion)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 1)
p_sniffles <- p_sniffles + xlab("Sniffles") + ylab("Proportion Non-Intersecting SV Calls")

## SVIM

svim_minimap2_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_minimap2_10X.df$Depth <- c("10X")
svim_minimap2_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_minimap2_20X.df$Depth <- c("20X")
svim_minimap2_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_minimap2_40X.df$Depth <- c("40X")
svim_minimap2_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_0/svim-minimap2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_minimap2_60X.df$Depth <- c("60X")
svim_minimap2_subsampled.df <- rbind(svim_minimap2_10X.df, svim_minimap2_20X.df, svim_minimap2_40X.df, svim_minimap2_60X.df)
svim_minimap2_subsampled.df$Aligner <- c("minimap2")
svim_minimap2_subsampled.df <- svim_minimap2_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

svim_ngmlr_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_ngmlr_10X.df$Depth <- c("10X")
svim_ngmlr_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_ngmlr_20X.df$Depth <- c("20X")
svim_ngmlr_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_ngmlr_40X.df$Depth <- c("40X")
svim_ngmlr_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_0/svim-ngmlr_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_ngmlr_60X.df$Depth <- c("60X")
svim_ngmlr_subsampled.df <- rbind(svim_ngmlr_10X.df, svim_ngmlr_20X.df, svim_ngmlr_40X.df, svim_ngmlr_60X.df)
svim_ngmlr_subsampled.df$Aligner <- c("ngmlr")
svim_ngmlr_subsampled.df <- svim_ngmlr_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

svim_pbmm2_10X.df <- data.frame(read.delim(file = '4_results/subsampled/10X/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_pbmm2_10X.df$Depth <- c("10X")
svim_pbmm2_20X.df <- data.frame(read.delim(file = '4_results/subsampled/20X/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_pbmm2_20X.df$Depth <- c("20X")
svim_pbmm2_40X.df <- data.frame(read.delim(file = '4_results/subsampled/40X/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_pbmm2_40X.df$Depth <- c("40X")
svim_pbmm2_60X.df <- data.frame(read.delim(file = '4_results/subsampled/60X/sv_intersection_agreement/svim/qual_0/svim-pbmm2_agreement_summary_total.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE))
svim_pbmm2_60X.df$Depth <- c("60X")
svim_pbmm2_subsampled.df <- rbind(svim_pbmm2_10X.df, svim_pbmm2_20X.df, svim_pbmm2_40X.df, svim_pbmm2_60X.df)
svim_pbmm2_subsampled.df$Aligner <- c("pbmm2")
svim_pbmm2_subsampled.df <- svim_pbmm2_subsampled.df %>% filter(., SV.TYPE == "ALL_SVS") %>% select(c(Aligner, Depth, Non.intersecting.proportion))

svim_all_aligners.df <- rbind(svim_minimap2_subsampled.df, svim_pbmm2_subsampled.df, svim_ngmlr_subsampled.df)
p_svim <- ggplot(svim_all_aligners.df, aes(x=Depth, y = Non.intersecting.proportion)) + geom_col( position=position_dodge()) + facet_grid(~ Aligner, scales = "free_x", space = "free_x", switch = "x") + ylim(0, 1) 
p_svim_no_y_title <- p_svim + xlab("svim")  + theme(axis.title.y=element_blank())
p_svim <- p_svim + xlab("svim") + ylab("Proportion Non-Intersecting SV Calls")

## Plot Sniffles and SVIM only
combined_plots_sniffles_svim <- p_sniffles + p_svim_no_y_title
combined_plots_sniffles_svim <- combined_plots_sniffles_svim 

## Plot pbsv, Sniffles, and SVIM

combined_plots_pbsv_sniffles_svim <- p_sniffles + p_svim_no_y_title + p_pbsv + plot_layout(nrow = 2, byrow = TRUE)
combined_plots_pbsv_sniffles_svim <- combined_plots_pbsv_sniffles_svim 

## Plot pbsv, Sniffles, and SVIM with counts

p_sniffles_counts <- p_sniffles + 
	geom_text(aes(x=Depth, y = Non.intersecting.proportion, label = Non.intersecting.proportion), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)

dir.create(file.path("5_plots","subsampled"), recursive = TRUE)
ggsave("2_sv_agreement_sniffles_counts_subsampled.png", plot = p_sniffles_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

p_svim_counts <- p_svim_no_y_title + geom_text(aes(x=Depth, y = Non.intersecting.proportion, label = Non.intersecting.proportion, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
ggsave("2_sv_agreement_svim_counts_subsampled.png", plot = p_sniffles_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

p_pbsv_counts <- p_pbsv + geom_text(aes(x=Depth, y = Non.intersecting.proportion, label = Non.intersecting.proportion, group = Depth), position = position_dodge(width = 1),vjust = -0.5, size = 2.5)
ggsave("2_sv_agreement_pbsv_counts.png", plot = p_pbsv_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")

combined_plots_pbsv_sniffles_svim_counts <- p_sniffles_counts + p_svim_counts + p_pbsv_counts + plot_layout(nrow = 2, byrow = TRUE)
combined_plots_pbsv_sniffles_svim_counts <- combined_plots_pbsv_sniffles_svim_counts + plot_annotation(title = 'Impact of Sequencing Depth on SV Agreement', theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("2_sv_agreement_pbsv_sniffles_svim_subsampled_counts.png", plot = combined_plots_pbsv_sniffles_svim_counts, device = "png", path = "5_plots/subsampled", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")
