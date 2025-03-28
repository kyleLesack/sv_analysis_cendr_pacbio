library(karyoploteR)
library(BSgenome.Celegans.UCSC.ce11)
library(regioneR)
genome <- getBSgenome("BSgenome.Celegans.UCSC.ce11")
dels_coords <- toGRanges("9_plots/sv_hotspots/ngmlr/sniffles/7_ins_to_dup/DEL.bed", genome="BSgenome.Celegans.UCSC.ce11" )
dups_coords <- toGRanges("9_plots/sv_hotspots/ngmlr/sniffles/7_ins_to_dup/DUP.bed", genome="BSgenome.Celegans.UCSC.ce11" )
dels_coords_40x <- toGRanges("9_plots/sv_hotspots/ngmlr.subsampled_40X/sniffles/7_ins_to_dup/DEL.bed", genome="BSgenome.Celegans.UCSC.ce11" )
dups_coords_40x <- toGRanges("9_plots/sv_hotspots/ngmlr.subsampled_40X/sniffles/7_ins_to_dup/DUP.bed", genome="BSgenome.Celegans.UCSC.ce11" )
head(dels_coords)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2)
kp <- kpPlotDensity(kp, dels_coords, window.size = 500000, col="#3388FF")
kp <- kpPlotDensity(kp, dups_coords, window.size = 500000, col="#ff4433",data.panel = 2)


kp_smaller_window <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2)
kp_smaller_window <- kpPlotDensity(kp_smaller_window, dels_coords, window.size = 100000, col="#3388FF")
kp_smaller_window <- kpPlotDensity(kp_smaller_window, dups_coords, window.size = 100000, col="#ff4433",data.panel = 2)

