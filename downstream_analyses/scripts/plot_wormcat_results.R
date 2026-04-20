library(tidyverse)
library(readr)
library(clipr)

cat1_file <- file.path(snakemake@input[["cat1"]])
cat1 <- read_csv(cat1_file)

cat1 <- cat1 %>% select(!overlapGenes) %>% mutate(pvalue_adjusted = p.adjust(pval, method='fdr')) %>%
  filter(pvalue_adjusted <= 0.05) %>% mutate(level = "C1") %>%
  mutate(pathway = str_replace(pathway, pattern = "GMT1_", replacement = "")) %>%
  arrange(desc(overlap)) %>% 
  mutate_if(is.character,as.factor) %>%
  mutate(Category = fct_reorder(pathway, overlap))

cat2_file <- file.path(snakemake@input[["cat2"]])
cat2 <- read_csv(cat2_file)

cat2 <- cat2 %>% select(!overlapGenes) %>% mutate(pvalue_adjusted = p.adjust(pval, method='fdr')) %>%
  filter(pvalue_adjusted <= 0.05) %>% mutate(level = "C2") %>%
  mutate(pathway = str_replace(pathway, pattern = "GMT1_", replacement = "")) %>%
  arrange(desc(overlap)) %>% 
  mutate_if(is.character,as.factor) %>%
  mutate(Category = fct_reorder(pathway, overlap))

cat3_file <- file.path(snakemake@input[["cat3"]])
cat3 <- read_csv(cat3_file)
cat3 <- cat3 %>% select(!overlapGenes) %>% mutate(pvalue_adjusted = p.adjust(pval, method='fdr')) %>%
  filter(pvalue_adjusted <= 0.05) %>% mutate(level = "C3") %>%
  mutate(pathway = str_replace(pathway, pattern = "GMT1_", replacement = "")) %>%
  arrange(desc(overlap)) %>% 
  mutate_if(is.character,as.factor) %>%
  mutate(Category = fct_reorder(pathway, overlap))

wormcat_all_levels <- rbind(cat1, cat2, cat3) %>%
  mutate(Category = fct_reorder(Category, overlap))

p <- ggplot(wormcat_all_levels, aes(x = overlap, y = Category, fill = level)) + geom_col() +
  facet_wrap(~ level, ncol = 1, scales = "free") +
  xlab("Overrepresented Genes") + geom_text(aes(label = formatC(padj, format = "e", digits = 2)),
                                            position=position_dodge(width=0.9), hjust=-.1, size = 8/.pt) +
  theme(legend.position="none") + ylab(NULL)
custom_palette <- c("#E69F00", "#56B4E9", "#999999")
p <- p + theme(panel.background = element_rect(fill = "white"), 
          panel.grid.major.x = element_line(color = "gray"), 
          panel.grid.minor.x = element_line(color = "gray"),
          panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values=custom_palette)
(p <- p + expand_limits(x = c(0,650)))
outfile <- file.path(snakemake@output[[1]])
ggsave(outfile, scale = 2)

outfile <- file.path(snakemake@output[[2]])
ggsave(outfile, scale = 2)