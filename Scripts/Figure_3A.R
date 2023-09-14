#,===============================================================================
#' Description:
#' Reproduce Figure 3A
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
#' Uses ComplexHeatmap and a set of matrices to plot Figure 3A
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
#'  Function to plot ecDNA and other complex events frequency and distribution
#'  in Cascade cohort. I excluded Linear amplification.
#===============================================================================

plot_genes_AA_ht <- function(genes_to_include, include_oncoprint = TRUE, row_split, complex_events,
                             sample_full, include_clinical_info = TRUE, matrix_oncoP,
                             id_order, extra_ha = c(), extra_lgd = c(), add_line_in_annotation = F, anno_str, line_value,
                             save = T, fn_output_plot, width = 15, height = 10) {

    library(methods)
    extra_clinical_annotation = c("psa_category", "gleason", "overall_survival")
    lgd_list <- list()
    clinical_ha <- c()
    if (include_clinical_info) {
      wgd_col_fun <- c("FALSE" ="white", "TRUE" = "black") 
      ploidy_col_fun <- colorRamp2(quantile(sample_full$ploidy_PURPLE), brewer.pal(n = 5, name = "Greys"))
      SVtmbPerMb_col_fun <-  colorRamp2(quantile(sample_full$SVtmbPerMb), brewer.pal(n = 5, name = "Greys"))
      chromo_col_fun <- c("FALSE" = "white", "TRUE" = "black")
      TP53_col_fun <- c("FALSE" = "white", "TRUE" = "black")
      PTEN_col_fun <- c("FALSE" = "white", "TRUE" = "black")
      event_col_fun <- c("ecDNA"="#EB8D71","BFB"="#009392","Complex"= "#9CCB86","Unknown" = "#EBCE8D")
      
      construct_pch <- function(l) {
        pch <- rep(NA, length(l))
        is_na <- is.na(l)
        pch[is_na] <- ""
        return(pch)
      }
      
      
      clinical_ha = HeatmapAnnotation(
        'Complex events'= anno_barplot(complex_events, gp =gpar(color='black', fill=unname(event_col_fun), lwd = 0.7 ),baseline = 0, height = unit(5, "cm")),
        'Whole Genome doubling' = anno_simple(sample_full$WG_doubling, col = wgd_col_fun, pch = construct_pch(sample_full$WG_doubling), na_col = "white", gp = gpar(col = "white")),
        'SV tmbPerMb' = anno_simple(sample_full$SVtmbPerMb, col = SVtmbPerMb_col_fun, pch = construct_pch(sample_full$SVtmbPerMb), na_col = "white", gp = gpar(col = "white")),
        Ploidy = anno_simple(sample_full$ploidy_PURPLE, col = ploidy_col_fun, pch = construct_pch(sample_full$ploidy_PURPLE), na_col = "white", gp = gpar(col = "white")),
        'Chromothripsis' = anno_simple(sample_full$Chromo, col = chromo_col_fun, pch = construct_pch(sample_full$Chromo), na_col = "white", gp = gpar(col = "white")),
        'TP53 alterations' = anno_simple(sample_full$TP53, col = chromo_col_fun, pch = construct_pch(sample_full$TP53), na_col = "white", gp = gpar(col = "white")),
        'PTEN alterations' = anno_simple(sample_full$PTEN, col = chromo_col_fun, pch = construct_pch(sample_full$PTEN), na_col = "white", gp = gpar(col = "white")),
        annotation_name_side = "left"
      )
        
      # add legend for annotations
      lgd_wgd = Legend(title = "WGD", labels = names(wgd_col_fun), legend_gp = gpar(fill = unname(wgd_col_fun)), border = "black")
      lgd_SVtmb = Legend(title = "SV tmbPerMb", labels = c("0", "10","100", "500", "1000"), at = c(0, 10,100, 500, 1000), col_fun = SVtmbPerMb_col_fun)
      lgd_ploidy = Legend(title = "Ploidy", labels = c("1", "2", "3", "4", "5"), at = c(1, 2, 3, 4, 5), col_fun = ploidy_col_fun)
      lgd_chromo = Legend(title = "Chromothripsis", labels = names(chromo_col_fun), legend_gp = gpar(fill = unname(chromo_col_fun)), border = "black")
      lgd_tp53 = Legend(title = "TP53 alt", labels = names(TP53_col_fun), legend_gp = gpar(fill = unname(TP53_col_fun)), border = "black")
      lgd_pten = Legend(title = "PTEN alt", labels = names(PTEN_col_fun), legend_gp = gpar(fill = unname(PTEN_col_fun)), border = "black")
      
        lgd_list <- c(lgd_list,  lgd_wgd,lgd_ploidy, 
                      lgd_SVtmb, lgd_chromo, lgd_tp53, lgd_pten )
    }
    
    
    # plot gene status as a heatmap
    for (gene in genes_to_include) {
      print(gene)
      
      
      sample_full[, gene] <- paste(sample_full[, gene], sample_full[, paste0(gene, "_mut")], 
                                   sample_full[, paste0(gene, "_germline")]
      )
    }
    genes_to_include <- c()
    ht_mat <- as.matrix(sample_full[,genes_to_include])
    rownames(ht_mat) <- sample_full$sample_id
    ht_mat[is.na(ht_mat)] = "NA"
    color_names <- unique(c(ht_mat))
    gene_col_fun <- rep(NA, length(color_names))
    names(gene_col_fun) <- color_names
    
    if (length(extra_ha) == 0) {
        ht <- Heatmap(t(ht_mat),
                      col = gene_col_fun, row_names_side = "left", show_row_names = T,
                      cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F, 
                      rect_gp = gpar(col = "black"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                          if(str_detect(t(ht_mat)[i, j], "TRUE")) {
                              grid.points(x, y, pch = 16, size = unit(2.5, "mm")) 
                          }
                      }, column_order = id_order,  top_annotation = clinical_ha )
    } else {
        ht <- Heatmap(t(ht_mat),
                      col = gene_col_fun, row_names_side = "left", show_row_names = T,
                      cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F, 
                      rect_gp = gpar(col = "black"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                          if(str_detect(t(ht_mat)[i, j], "TRUE")) {
                              grid.points(x, y, pch = 16, size = unit(2.5, "mm")) 
                          }
                      }, column_order = id_order,  top_annotation = c(extra_ha, clinical_ha)) 
    }
    
    if (include_oncoprint) {
 
      col_oncoP = c( "NULL"= "white","ecDNA" = "#EB8D71","BFB" = "#009392",
                     "unknown" = "#EBCE8D","Complex" =  "#9CCB86")
      event_col_fun <- c("ecDNA"="#EB8D71","BFB"="#009392","Complex"= "#9CCB86","Unknown" = "#EBCE8D")
      alter_fun = list(
            background = function(x, y, w, h) 
                grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#f0f0f0", col = NA, lwd = 1.5)),
            
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
        
        oncoP<- oncoPrint(matrix_oncoP, top_annotation = NULL, remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                          right_annotation = NULL, column_order = id_order, alter_fun_is_vectorized = FALSE, row_gap = unit(1.5, "mm"),
                          row_split = row_split,  border = TRUE, row_title_rot = 0, row_title_gp = gpar(fontsize = 10),
                          alter_fun = alter_fun, col = col_oncoP, show_pct = FALSE,row_names_side = "left", show_column_names = TRUE)
        plots_list =  ht %v% oncoP 
        ComplexHeatmap::draw(plots_list, annotation_legend_list = lgd_list, merge_legend = T)
    }
    else {
        plots_list =  ht
        ComplexHeatmap::draw(plots_list, annotation_legend_list = lgd_list, merge_legend = T)
    }
    if (add_line_in_annotation) {
        decorate_annotation(anno_str, {
            grid.lines(c(0, 1), unit(c(line_value, line_value), "native"), gp = gpar(col = "red", lty = 2))
        })
    }
    if (save) {
        pdf(fn_output_plot, width = width, height = height)
        ComplexHeatmap::draw(plots_list,annotation_legend_list = lgd_list, merge_legend = T)
        if (add_line_in_annotation) {
            decorate_annotation(anno_str, {
                grid.lines(c(0, 1), unit(c(line_value, line_value), "native"), gp = gpar(col = "red", lty = 2))
            })
        }
        dev.off()
    }
}


#===============================================================================
# To use the function plot_genes_AA_ht, split_ecDNA_ca_gene_list, id_order_ca,
# cascade_events_mod, SO_summary_sub, and matrix_AA_gene_ca_freq_sub_new_mod 
# DFs are necessary!


plot_genes_AA_ht(genes_to_include,  include_oncoprint = T, row_split = split_ecDNA_ca_gene_list[,1],complex_events = cascade_events_mod[,3:6],
                 sample_full = SO_summary_sub, include_clinical_info = T, matrix_oncoP = matrix_AA_gene_ca_freq_sub_new_mod,
                 id_order = id_order_ca, save = T, add_line_in_annotation = F,
                 fn_output_plot = "Figure_3A.pdf", width = 15, height = 10)

