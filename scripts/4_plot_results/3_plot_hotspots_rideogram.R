library(RIdeogram)
# Read optimal window size https://www.biostars.org/p/17899/ and C. elegans gene density https://www.ncbi.nlm.nih.gov/books/NBK20219/

library(karyoploteR)
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


