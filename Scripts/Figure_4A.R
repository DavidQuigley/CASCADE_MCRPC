#,===============================================================================
#' Description:
#' Reproduce the Figure 4A
#' Adapted from David's script to plot frequency of ecDNA at the AR locus
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Authors: 
#' - David Quigley & Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Uses Gviz and AA & AC outputs to plot Figure 4A
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. 
#===============================================================================

library( GenomicRanges )
library( Gviz )
library( biomaRt )

#===============================================================================
#'  Creating the tracks to plot Figure 4A 
#===============================================================================

bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# gr_EC_counts_X
positions = c(66000000, 69000000)
biomTrack <- BiomartGeneRegionTrack(genome = "GRCh38", chromosome = "chrX",
                                    start = positions[1], end = positions[2],
                                    filters=list("ensembl_transcript_id" = c("ENST00000374690","ENST00000374719", "ENST00000680612")),
                                    transcriptAnnotation = "symbol",
                                    fontcolor.feature = "black",
                                    name = "",
                                    biomart = bm)

display_properties = list(background.title = "white",
                          col.title = "black",
                          col.axis="black")

gtrack <- GenomeAxisTrack()
dTrack <- DataTrack(gr_EC_counts_X, type=c("mountain", "g"),
                    name = "ecDNA counts", col="black",
                    fill.mountain=c("lightgrey", "lightgrey"))

aTrack <- AnnotationTrack(start = c(66900000), width = 1000,
                          chromosome = "chrX",
                          strand = c("+"), just.group = 'left',
                          fontcolor.feature = "black",
                          id = c("AR enhancer"),
                          groupAnnotation = "id", 
                          fontcolor.feature = 1, 
                          name = "",
                          genome = "hg38")

displayPars(aTrack) <- display_properties
displayPars(biomTrack) = display_properties
displayPars(dTrack) <- display_properties

pdf( 'Figure_4A.pdf', height=5, width=14)
plotTracks( list( dTrack, biomTrack, aTrack, gtrack), cex = 1, #main = "Chr X", 
            from=positions[1], to = positions[2], col.main = "black", cex.main = 0.8)
dev.off()
