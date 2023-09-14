#,===============================================================================
#' Description:
#' Reproduce Supplementary Figure 3B
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
#' Uses ComplexHeatmap and data from WCDT cohort
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. 
#===============================================================================

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
require(gridExtra)
library(Sushi)
library(Gviz)
library(GenomicFeatures)
library(reshape2)
library(ggrepel)
library(rtracklayer)
require(methods)

#===============================================================================
#'  Function to plot copy number, expression levels and distribution of ecDNA 
#'  and other complex events per selected genes, as well as AR and NE scores for
#'  WCDT cohort ( AR, MYC, and PCAT1 ). I excluded Linear amplifications.
#'  I also removed the legend from this figures, setting show_heatmap_legend to 
#'  FALSE in line 90 and 103, and commented line 143 and 159
#===============================================================================

plot_rna_ecDNA_wcdt_annotations_ht <- function(include_oncoprint = TRUE, add_line_in_annotation = F, anno_str, line_value,
                                      sample_full, include_clinical_info = TRUE, matrix_oncoP, exp_mat, complex_events,
                                      id_order, sv_matrix, mut_sig, extra_ha = c(), extra_lgd = c(),
                                      save = T, fn_output_plot, width = 15, height = 10) {

  library(methods)
  lgd_list <- list()
  clinical_ha <- c()
  if (include_clinical_info) {
    event_col_fun <- c("ecDNA"="#EB8D71","BFB"="#009392","Complex"= "#9CCB86", "Unknown" = "#EBCE8D") 
    AR_score_col_fun <- colorRamp2(c(quantile(sample_full$AR_score, 0), quantile(sample_full$AR_score, 0.5), quantile(sample_full$AR_score, 1)), c("white","#D9D9D9","black"))#colorRamp2(c(quantile(sample_full$AR_score, 0), quantile(sample_full$AR_score, 0.5), quantile(sample_full$AR_score, 1)), c("green","white", "red"))
    NE_score_col_fun <- colorRamp2(c(quantile(sample_full$NE_score, 0), quantile(sample_full$NE_score, 0.5), quantile(sample_full$NE_score, 1)), c("white","#D9D9D9","black"))#brewer.pal(n = 3, name = "Greys"))
    AR_log2_cnv_col_fun <- colorRamp2(c(0, 4), c("#FFFFFF", "#df8640"))
    MYC_log2_cnv_col_fun <- colorRamp2(c(0, 4), c("#FFFFFF", "#Df8640"))
    PCAT1_Log2_cnv_col_fun <- colorRamp2(c(0, 4), c("#FFFFFF", "#Df8640"))
    construct_pch <- function(l) {
      pch <- rep(NA, length(l))
      is_na <- is.na(l)
      pch[is_na] <- ""
      return(pch)
    }

    clinical_ha = HeatmapAnnotation(
      'Complex\nevents'= anno_barplot(complex_events, gp =gpar(color='black', fill=unname(event_col_fun), lwd = 0.7 ),baseline = 0, height = unit(2, "cm")),
      'AR score' = anno_simple(sample_full$AR_score, col = AR_score_col_fun, pch = construct_pch(sample_full$AR_score), na_col = "white", gp = gpar(col = "white")),
      'NE score' = anno_simple(sample_full$NE_score, col = NE_score_col_fun, pch = construct_pch(sample_full$NE_score), na_col = "white", gp = gpar(col = "white")),
      'empty' = anno_empty(border = FALSE, height = unit(0.25, "mm")),
      'AR cn' = anno_simple(sample_full$AR_Log2_cnv, col = AR_log2_cnv_col_fun, pch = construct_pch(sample_full$AR_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'MYC cn' = anno_simple(sample_full$MYC_Log2_cnv, col = MYC_log2_cnv_col_fun, pch = construct_pch(sample_full$MYC_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'PCAT1 cn' = anno_simple(sample_full$PCAT1_Log2_cnv, col = MYC_log2_cnv_col_fun, pch = construct_pch(sample_full$PCAT1_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'empty_2' = anno_empty(border = FALSE, height = unit(0.25, "mm")),
      annotation_name_side = "left"
    )
    
    # add legend for annotations
    lgd_ar_score = Legend(title = "AR/NE\nscore", labels = c("Low","High"), at = c(-6,7), col_fun = AR_score_col_fun)
    lgd_log2_cnv = Legend(title ="log2(cn)", labels = c("0","1","2",">3.3"), at = c(0, 1, 2, 3.3), col_fun = AR_log2_cnv_col_fun)
    lgd_list <- c(lgd_list, 
                  lgd_ar_score,
                  lgd_log2_cnv 
                  )
    
  }
  ht_mat <- as.matrix(exp_mat)
  if (length(extra_ha) == 0) {
    ht <- Heatmap(ht_mat, name = " z-score\nlog2(TPM +1)",  
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F, 
                  rect_gp = gpar(col = "white"), border = TRUE, 
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 8, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 4),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                  },  top_annotation = clinical_ha , column_order = id_order)
  } 
  else {
    ht <- Heatmap(ht_mat, name = " z-score\nlog2(TPM +1)",
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F,
                  rect_gp = gpar(col = "white"), border = TRUE, 
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 8, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 4),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                    }, column_order = id_order,  top_annotation = c(extra_ha, clinical_ha) ) 
  }
  if (include_oncoprint) {
    
    col_oncoP = c( "NULL"= "white","ecDNA" = "#EB8D71","BFB" = "#009392",
                   "unknown" = "#EBCE8D","Complex" =  "#9CCB86")
    alter_fun = list(
      background = function(x, y, w, h) 
        grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#f0f0f0",col = NA, lwd = 1.5)),
      'ecDNA' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9, 
                  gp = gpar(fill = col_oncoP["ecDNA"], col = "black", lwd = 0.5))
      },
      'BFB' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9, 
                  gp = gpar(fill = col_oncoP["BFB"], col = "black", lwd = 0.5))
      },
      'unknown' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9, 
                  gp = gpar(fill = col_oncoP["unknown"], col = "black", lwd = 0.5))
      },
      'Complex' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9,
                  gp = gpar(fill = col_oncoP["Complex"], col = "black", lwd = 0.5))
      }
    )
    oncoP<- oncoPrint(matrix_oncoP, top_annotation = NULL,remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                      column_order = id_order, border = TRUE, row_order = 1:nrow(matrix_oncoP),
                      alter_fun = alter_fun, col = col_oncoP, pct_side = F, alter_fun_is_vectorized = FALSE, row_gap = unit(1, "mm"),
                      row_names_side = "left", show_column_names = TRUE , right_annotation = NULL ) 
    plot_list =  ht %v% oncoP 
    # ComplexHeatmap::draw(plot_list,   annotation_legend_list = lgd_list, merge_legend = T)
    ComplexHeatmap::draw(plot_list, merge_legend = T)
  }
  else {
    ComplexHeatmap::draw(ht, annotation_legend_list = lgd_list)
  }
  
  if (add_line_in_annotation) {
    decorate_annotation(anno_str, {
      grid.lines(c(0, 1), unit(c(line_value, line_value), "native"), gp = gpar(col = "red", lty = 2))
    })
  }
  
  if (save) {
    pdf(fn_output_plot, width = width, height = height)
      plot_list =  ht %v% oncoP 
      # ComplexHeatmap::draw(plot_list,   annotation_legend_list = lgd_list, merge_legend = T)
      ComplexHeatmap::draw(plot_list,  merge_legend = T)
      
    if (add_line_in_annotation) {
      decorate_annotation(anno_str, {
        grid.lines(c(0, 1), unit(c(line_value, line_value), "native"), gp = gpar(col = "red", lty = 2))
      })
    }
    dev.off()
  }
}


#===============================================================================
# To use the function plot_rna_ecDNA_wcdt_annotations_ht, I re-ordered Log2.TPM_wcdt.zc_mod
# matrix_AA_gene_wcdt_freq_sub_new_mod, wcdt_events_mod, and RNA_wcdt dataframes 
# based on AR expression

plot_rna_ecDNA_wcdt_annotations_ht( include_oncoprint = T, exp_mat = Log2.TPM_wcdt.zc_mod[c("AR","MYC","PCAT1"),id_order_ar_wcdt_mod],
                                    sample_full = RNA_wcdt, include_clinical_info = T, id_order = id_order_ar_wcdt_mod, 
                                    matrix_oncoP = matrix_AA_gene_wcdt_freq_sub_new_mod[c("AR","MYC","PCAT1"),id_order_ar_wcdt_mod],
                                    save = T, add_line_in_annotation = F, complex_events =  wcdt_events_mod[id_order_ar_wcdt_mod,3:6],
                                    fn_output_plot = "Supplementary_Figure_3B.pdf", width = 12, height = 4.5) 


