#,===============================================================================
#' Description:
#' Reproduce the Figures 2 b-e
#'
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Author: 
#' - Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Uses CopyNumberPlots and a set of matrices to plot Figure 2b-e
#****************************************************************************
################################################################################

#===============================================================================
#'  Creates an GRanges object to plot the AR enhancer 
#===============================================================================

library(CopyNumberPlots)
library(karyoploteR)
library(BSgenome)

genome_cn <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
# Zoom-in to only plot AR gene on the ideogram
zoom.region <- toGRanges(data.frame("chrX", 67.5e6, 67.8e6), genome = genome_cn)
kp <- plotKaryotype(chromosomes="chrX", zoom=zoom.region, genome = genome_cn)
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

# To zoom in a lager proportion of ChrX
zoom.region <- toGRanges(data.frame("chrX", 64.5e6, 69.5e6))
markers_chrX <- data.frame(chr=rep("chrX", 1), pos=c(66920021), labels=c("AR enhancer"))

# Plotting parameters
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 40
pp$data1height <- 500
pp$data1inmargin <- 30
pp$data1outmargin <- 30
pp$topmargin <- 20
pp$rightmargin <- 0.02
pp$leftmargin <- 0.18

#===============================================================================
#'  Figure 2b, patient CA071
#===============================================================================

cn.calls <- list("CA071-01"=SO[[1]]$segments, "CA071-02"=SO[[2]]$segments, "CA071-03"=SO[[3]]$segments,
                 "CA071-04"=SO[[4]]$segments, "CA071-05"=SO[[5]]$segments, "CA071-06"=SO[[6]]$segments,
                 "CA071-07"=SO[[7]]$segments, "CA071-08"=SO[[8]]$segments, "CA071-09"=SO[[9]]$segments)

cn.CA071_calls <- lapply(cn.calls, function(x) x<- loadCopyNumberCalls(x, genome = genome_cn) )

pdf("Figure_2b.pdf", 
    height=8, width=6)

kp <- plotKaryotype(chromosomes="chrX", zoom=zoom.region, genome = genome_cn, main = "CA071",plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 5e+05)
kpAddCytobandLabels(kp, cex = 0.7)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")

kpPlotMarkers(kp, chr=markers_chrX$chr, r0=-2.5, x=markers_chrX$pos, labels=markers_chrX$labels,adjust.label.position=TRUE,
              text.orientation = "horizontal",line.color = "red", label.color = "black", y=0.75, pos = 2, offset = 0.5,
              cex=1)

plotCopyNumberCallsAsLines(kp, cn.CA071_calls, r0=0.1, cn.colors = "red_blue", loh.color = "orange", r1=1,
                           style = "line",#add.axis=F, 
                           label.cex = 0.8, add.axis=F, numticks = 2, axis.cex = 5 )
kpText(kp, chr="chrX", x=69.18e6, y=c(0.15,0.24,0.33, 0.42, 0.51,  0.61, 0.70, 0.8, 0.9), cex = 0.8,
       labels=c("max 76.9", "max 132.6", "max 60.5", "max 78.9", "max 97.7","max 70.5","max 74.1","max 110.2", "max 69.4"), col="black", pos=3)

dev.off()
#===============================================================================
#'  Figure 2c, patient CA104
#===============================================================================
grep( "CA104", names(SO))
cn.calls <- list("CA104-01"=SO[[33]]$segments, "CA104-05"=SO[[34]]$segments, "CA104-08"=SO[[35]]$segments,
                 "CA104-09"=SO[[36]]$segments, "CA104-10"=SO[[37]]$segments, "CA104-26"=SO[[38]]$segments,
                 "CA104-27"=SO[[39]]$segments, "CA104-28"=SO[[40]]$segments, "CA104-29"=SO[[42]]$segments,
                 "CA104-33"=SO[[42]]$segments, "CA104-35"=SO[[43]]$segments, "CA104-40"=SO[[44]]$segments,
                 "CA104-41"=SO[[45]]$segments)

cn.CA104_calls <- lapply(cn.calls, function(x) x<- loadCopyNumberCalls(x, genome = genome_cn) )

pdf("Figure_2c.pdf", 
    height=8, width=6)
kp <- plotKaryotype(chromosomes="chrX", zoom=zoom.region, genome = genome_cn, 
                    main = "CA104", plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 5e+05)
kpAddCytobandLabels(kp, cex = 0.7)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")

kpPlotMarkers(kp, chr=markers_chrX$chr, r0=-2.5, x=markers_chrX$pos, labels=markers_chrX$labels,adjust.label.position=TRUE,
              text.orientation = "horizontal",line.color = "red", label.color = "black", y=0.75, pos = 2, offset = 0.5,
              cex=1)

plotCopyNumberCallsAsLines(kp, cn.CA104_calls, r0=0.15, cn.colors = "red_blue", loh.color = "orange", r1=1,
                           style = "line", plot.params = pp,#add.axis=F, 
                           label.cex = 0.8, add.axis=F, numticks = 2, axis.cex = 5 )
kpText(kp, chr="chrX", x=69.18e6, y=c(0.14,0.21,0.27, 0.34, 0.42, 0.48, 0.54, 0.60, 0.66, 0.73, 0.8, 0.88, 0.94), 
       labels=c("max 2.1", "max 1", "max 1", "max 2.1", "max 2.2","max 3.9","max 1.1","max 3.7", "max 1", "max 2","max 1", "max 1.2","max 2"), col="black", pos=3)

dev.off()


#===============================================================================
#'  for Figure 2d. Patient CA113
#===============================================================================

grep( "CA113", names(SO))
cn.calls <- list("CA113-01"=SO[[15]]$segments, "CA113-05"=SO[[16]]$segments, "CA113-10"=SO[[17]]$segments,
                 "CA113-11"=SO[[18]]$segments, "CA113-14"=SO[[19]]$segments, "CA113-18"=SO[[20]]$segments,
                 "CA113-20"=SO[[21]]$segments, "CA113-25"=SO[[22]]$segments )

cn.CA113_calls <- lapply(cn.calls, function(x) x<- loadCopyNumberCalls(x, genome = genome_cn) )

pdf("Figure_2d.pdf", 
    height=8, width=6)

kp <- plotKaryotype(chromosomes="chrX", zoom=zoom.region, genome = genome_cn, main = "CA113",plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 5e+05)
kpAddCytobandLabels(kp, cex = 0.7)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")

kpPlotMarkers(kp, chr=markers_chrX$chr, r0=-2.5, x=markers_chrX$pos, labels=markers_chrX$labels,adjust.label.position=TRUE,
              text.orientation = "horizontal",line.color = "red", label.color = "black", y=0.75, pos = 2, offset = 0.5,
              cex=1)

plotCopyNumberCallsAsLines(kp, cn.CA113_calls, r0=0.15, cn.colors = "red_blue", loh.color = "orange", r1=1,
                           style = "line",#add.axis=F, 
                           label.cex = 0.8, add.axis=F, numticks = 2, axis.cex = 5 )
kpText(kp, chr="chrX", x=69.18e6, y=c(0.18,0.30, 0.40, 0.5,  0.63, 0.70, 0.82, 0.92), cex = 0.8,
       labels=c("max 9.1", "max 2.6", "max 0.9", "max 1.1", "max 0.9","max 2","max 2.6","max 1.7"), col="black", pos=3)

dev.off()

#===============================================================================
#'  for Figure 2e. Patient CA108
#===============================================================================

grep( "CA108", names(SO))

cn.calls <- list("CA108-13"=SO[[46]]$segments, "CA108-39"=SO[[47]]$segments, "CA108-46"=SO[[48]]$segments,
                 "CA108-67"=SO[[49]]$segments, "CA108-74"=SO[[50]]$segments, "CA108-82"=SO[[51]]$segments,
                 "CA108-86"=SO[[52]]$segments, "CA108-94"=SO[[53]]$segments )

cn.CA108_calls <- lapply(cn.calls, function(x) x<- loadCopyNumberCalls(x, genome = genome_cn) )

pdf("Figure_2e.pdf", 
    height=8, width=6)

kp <- plotKaryotype(chromosomes="chrX", zoom=zoom.region, genome = genome_cn, main = "CA108",plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 5e+05)
kpAddCytobandLabels(kp, cex = 0.7)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.2, cex=0.8, gene.name.position = "left")

kpPlotMarkers(kp, chr=markers_chrX$chr, r0=-2.5, x=markers_chrX$pos, labels=markers_chrX$labels,adjust.label.position=TRUE,
              text.orientation = "horizontal",line.color = "red", label.color = "black", y=0.75, pos = 2, offset = 0.5,
              cex=1)

plotCopyNumberCallsAsLines(kp, cn.CA108_calls, r0=0.15, cn.colors = "red_blue", loh.color = "orange", r1=1,
                           style = "line",#add.axis=F, 
                           label.cex = 0.8, add.axis=F, numticks = 2, axis.cex = 5 )
kpText(kp, chr="chrX", x=69.18e6, y=c(0.18,0.28, 0.38, 0.5,  0.58, 0.72, 0.81, 0.92), cex = 0.8,
       labels=c("max 6.6", "max 7.3", "max 8.9", "max 13.2", "max 4.1","max 110.5","max 13","max 13.6"), col="black", pos=3)

dev.off()
