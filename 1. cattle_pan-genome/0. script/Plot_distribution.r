library("karyoploteR")
library(rtracklayer)
setwd("E:\\OR_project\\GCF_002263795.3_ARS-UCD2.0_genomic.gff")
karyotype=read.table("kayotype.txt",head=F)
gff.file="GCF_002263795.3_ARS-UCD2.0_genomic.OR.gff"
features <- import(gff.file)
BS_genome <- toGRanges(karyotype)
kp <- plotKaryotype(genome=BS_genome)
table(features$type)
genes <- features[features$type=="gene" ]
pseudogene <- features[features$type=="pseudogene"]
pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
png( 
    filename = "OR_distribution.png", # 文件名称
    width = 30,            # 宽
    height = 10,           # 高
    units = "in",          # 单位
    bg = "white",          # 背景颜色
    res = 300)             # 分辨率
kp <- plotKaryotype(genome=BS_genome, ideogram.plotter = NULL, plot.type=2, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "Olfactory receptor genes or pseudogenes", cex=2)
kpPlotRegions(kp, data=genes, avoid.overlapping = FALSE, col="deepskyblue")
kpPlotRegions(kp, data=pseudogene, avoid.overlapping = FALSE, col="gold", data.panel=2)
dev.off()
