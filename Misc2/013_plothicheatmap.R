## This script plot heatmap from .hic file
## The .hic file for MB288 was downloaded from GSE240410
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse240410
suppressPackageStartupMessages({
  library(plotgardener)
  library(RColorBrewer)
})

DIR="HIC/"
DIRO="HIC/hicT/"
sam <- "MB288"
setwd(DIR)

pdf(paste0(DIRO,sam,"_chr11_hicT10kb.pdf"),width = 4,height = 3, pointsize = 30)
pageCreate(width = 4, height = 3, default.units = "inches",showGuides = FALSE)

hicPlot <- plotHicTriangle(
  data = paste0("GSE240410_data/",sam,".allValidPairs.hic"),
  resolution = 10000,
  zrange = c(1, 100),
  chrom = "11", 
  chromstart = 31000000, chromend = 32000000,
  matrix="log2oe",
  palette = colorRampPalette(brewer.pal(n = 9, "YlOrRd")),
  assembly = "hg38",
  x = 2, y = 0.5, width = 3, height = 2,
  just = "top", default.units = "inches"
)
annoHeatmapLegend(
  plot = hicPlot, x = 3.5, y = 0.5,
  width = 0.12, height = 1.2,
  just = c("left", "top"), default.units = "inches"
)

annoGenomeLabel(
  plot = hicPlot, scale = "Mb", axis = "x",
  x = 0.5, y = 2.53, just = c("left", "top")
)
pageGuideHide()
dev.off()
q()
