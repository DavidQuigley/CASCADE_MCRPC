#,===============================================================================
#' Description:
#' Reproduce the Supplementary Figure 2
#' Adapted from David's script to plot reversion mutations in patient CA071-04
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Authors: 
#' - David Quigley & Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Uses aardvark output and plot_reversion_summary function to plot Supp Fig 2
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. 
#===============================================================================
# install_github("DavidQuigley/aardvark")
library(ggplot2)
library(GenomicAlignments)
library(Rsamtools)
library(Biostrings)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(aardvark)
library(poppy)
library(VariantAnnotation)
library(RColorBrewer)
library(VariantAnnotation)

#===============================================================================
#' Supplementary Figure 2
#===============================================================================
fn_reversion = "Supplementary_Figure_2_CA071_04.pdf"

pdf( fn_reversion, height=12, width=14)
aardvark::plot_reversion_summary(sums_all[[4]]$summary,genome_version = 38,
                                 hsapiens_object = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,biomart_object = bm,
                                 pos_start = 32337500, pos_end=32341000,exclude_blacklist = TRUE, min_freq = 3, text_size_adjustment = 0.8)


dev.off()
