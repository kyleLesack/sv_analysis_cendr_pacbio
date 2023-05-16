#BiocManager::install("RIdeogram")
#library(RIdeogram)
# Read optimal window size https://www.biostars.org/p/17899/ and C. elegans gene density https://www.ncbi.nlm.nih.gov/books/NBK20219/
#require(RIdeogram)
setwd("/bulk/worm_lab/mrkyle/sv_analysis_cendr_pacbio")
#c_elegans_karyotype <- read.table("9_plots/sv_hotspots/c_elegans_karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
#gene_density <- GFFex(input = "9_plots/sv_hotspots/Caenorhabditis_elegans.WBcel235.109.gff3.gz", karyotype = "9_plots/sv_hotspots/c_elegans_karyotype.txt", feature = "gene", window = 100000)
#head(gene_density)
#del_density <- GFFex(input = "9_plots/sv_hotspots/ngmlr.subsampled_40X/sniffles/7_ins_to_dup/DEL.gff", karyotype = "9_plots/sv_hotspots/c_elegans_karyotype.txt", feature = "DEL", window = 100000)
#head(del_density)
#ideogram(karyotype = c_elegans_karyotype, overlaid = gene_density, label = LTR_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
#ideogram(karyotype = c_elegans_karyotype)
#convertSVG("chromosome.svg", device = "png")
#ideogram(karyotype = c_elegans_karyotype, overlaid = gene_density, width = 100)
#convertSVG("chromosome.svg", device = "png")
#ideogram(karyotype = c_elegans_karyotype, overlaid = del_density, width = 100)
#convertSVG("chromosome.svg", device = "png")
#ideogram(karyotype = c_elegans_karyotype, overlaid = gene_density, label = del_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
#ideogram(karyotype = c_elegans_karyotype, overlaid = gene_density, label = del_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
#data(human_karyotype, package="RIdeogram") #reload the karyotype data
#ideogram(karyotype = human_karyotype, overlaid = gene_density, label = LTR_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
#convertSVG("chromosome.svg", device = "png")

library(karyoploteR)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- genes(txdb)
head(all.genes)

if (interactive()) {
	if (!require("BiocManager"))
		install.packages("BiocManager")
	BiocManager::install("BSgenome.Celegans.UCSC.ce11")
}

library(BSgenome.Celegans.UCSC.ce11)
library(regioneR)
dd <- toGRanges("data.bed")
dd <- split(dd, f = dd$name)

genome <- getBSgenome("BSgenome.Celegans.UCSC.ce11")
kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11")
genes(genome)


