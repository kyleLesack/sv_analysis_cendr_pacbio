library(ggplot2)
library(tidyverse)
library(patchwork)
setwd("/bulk/worm_lab/mrkyle/sv_analysis_cendr_pacbio/")

# Open input data
sv_call_summary <- read.delim(file = '8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/sv_vep_summary.csv', sep = ',',stringsAsFactors= TRUE, header = TRUE)

# Convert to data frames and change NAs to 0
sv_call_summary.df <- data.frame(sv_call_summary)
sv_call_summary.df <- sv_call_summary.df %>% dplyr::mutate(INV.VEP = replace_na(INV.VEP, 0))

# Select columns containing SV Counts (no VEP) for each strain
sv_call_summary.df <- sv_call_summary.df %>% select(c(Strain,DEL.VEP,DUP.VEP, INV.VEP))
sv_call_summary.df <- sv_call_summary.df %>% rename(DEL = DEL.VEP) %>% rename(DUP = DUP.VEP) %>% rename(INV = INV.VEP)
sv_call_summary_long.df <- pivot_longer(sv_call_summary.df, cols = 2:4, names_to ="SV Type", values_to = "Count")
sv_call_summary_long.df

p_grouped_by_strain <- ggplot(data=sv_call_summary_long.df, aes(x=Strain, y=Count, fill=`SV Type`)) + geom_col(position=position_dodge()) 
p_grouped_by_strain <- p_grouped_by_strain + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p_grouped_by_strain + theme(axis.text.x = element_text(angle = 22.5, hjust = 1,size = 12))
p_grouped_by_strain <- p_grouped_by_strain + xlab("Strain") + scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0.01,0.01)))  
p_grouped_by_strain

# Add counts to columns
p_grouped_by_strain_counts <- p_grouped_by_strain + geom_text(aes(x=Strain, y = Count, label = Count), position = position_dodge(width = 1), hjust = -0.25, size = 2.5, angle = 90) 
#p_grouped_by_strain <- p_grouped_by_strain + geom_text(aes(x=Strain, y = Count, label = Count), position = position_dodge(width = 1),vjust = -0.5, size = 2.5, angle = 90)
p_grouped_by_strain_counts <- p_grouped_by_strain_counts +  scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0.01,0.075)))  
p_grouped_by_strain_counts

# Save
dir.create(file.path("9_plots","sv_counts", "7_ins_to_dup"), recursive = TRUE)
ggsave("2_svs_vep_high_impact_grouped_strain.png", plot = p_grouped_by_strain, device = "png", path = "9_plots/sv_counts/7_ins_to_dup", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1715, units = "px")
ggsave("2_svs_vep_high_impact_grouped_strain_counts.png", plot = p_grouped_by_strain_counts, device = "png", path = "9_plots/sv_counts/7_ins_to_dup", scale = 1, dpi = 300, limitsize = TRUE, width = 2605, height = 1800, units = "px")
