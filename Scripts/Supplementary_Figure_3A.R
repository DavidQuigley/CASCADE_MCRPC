#,===============================================================================
#' Description:
#' Reproduce Supplementary Figure 3A
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
#' Uses ComplexHeatmap and cascade matrices to plot Supplementary Figure 3A
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
#'  and other complex events per selected genes. In this case AR, MYC, and PCAT1
#'  in Cascade cohort. I excluded Linear amplifications.
#===============================================================================
plot_rna_ecDNA_annotations_ht <- function(include_oncoprint = TRUE, add_line_in_annotation = F, anno_str, line_value,
                                      sample_full, include_clinical_info = TRUE, matrix_oncoP, exp_mat, complex_events,
                                      id_order, sv_matrix, mut_sig, extra_ha = c(), extra_lgd = c(),
                                      save = T, fn_output_plot, width = 15, height = 10) {
  library(methods)
  lgd_list <- list()
  clinical_ha <- c()
  if (include_clinical_info) {
    event_col_fun <- c("ecDNA"="#EB8D71","BFB"="#009392","Complex"= "#9CCB86", "Unknown" = "#EBCE8D")
    patient_id_col_fun <- c("CA071" = "#0D0887FF", "CA093" = "#6A00A8FF", "CA101" = "#B12A90FF", "CA104" = "#E16462FF", "CA108" = "#FCA636FF", "CA113" = "#F0F921FF")
    site_col_fun <- c(rcartocolor::carto_pal (8, "Prism"));names(site_col_fun) <- c("Bone", "Prostate", "Spine", "LN", "Brain", "Liver", "Lung", "Others") 
    asi_drug_col_fun <- c("none" = "white", "abi" = "#D9D9D9", "enz" = "#737373", "both" = "#252525")
    Cons_clust_col_fun <- c( pal_okabe_ito[c(-1, -5, -8)]);names(Cons_clust_col_fun)<- c('AR+/NE+','AR+/NE-','ARL/NE-','AR-/NE+','AR-/NE-')
    AR_score_col_fun <- colorRamp2(c(quantile(sample_full$AR_score, 0), quantile(sample_full$AR_score, 0.5), quantile(sample_full$AR_score, 1)), c("white","#D9D9D9","black"))
    NE_score_col_fun <- colorRamp2(c(quantile(sample_full$NE_score, 0), quantile(sample_full$NE_score, 0.5), quantile(sample_full$NE_score, 1)), c("white","#D9D9D9","black"))
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
      'Patient ID' = anno_simple(sample_full$patient_id, col = patient_id_col_fun, pch = construct_pch(sample_full$patient_id), na_col = "white", gp = gpar(col = "white"), border = TRUE),
      "ASI" = anno_simple(sample_full$ASI_drug, col = asi_drug_col_fun, pch = construct_pch(sample_full$ASI_drug), na_col = "white",gp = gpar(col = "white")),
      'AR score' = anno_simple(sample_full$AR_score, col = AR_score_col_fun, pch = construct_pch(sample_full$AR_score), na_col = "white", gp = gpar(col = "white")),
      'NE score' = anno_simple(sample_full$NE_score, col = NE_score_col_fun, pch = construct_pch(sample_full$NE_score), na_col = "white", gp = gpar(col = "white")),
      'Metastasis site'= anno_simple(sample_full$Tissue, col = site_col_fun, pch = construct_pch(sample_full$Tissue), na_col = "white", gp = gpar(col = "white")),
      'subtypes'= anno_simple(sample_full$Cluster_ID, col = Cons_clust_col_fun, pch = construct_pch(sample_full$Cluster_ID), na_col = "white", gp = gpar(col = "white")),
      'empty' = anno_empty(border = FALSE, height = unit(0.25, "mm")),
      'AR cn' = anno_simple(sample_full$AR_Log2_cnv, col = AR_log2_cnv_col_fun, pch = construct_pch(sample_full$AR_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'MYC cn' = anno_simple(sample_full$MYC_Log2_cnv, col = MYC_log2_cnv_col_fun, pch = construct_pch(sample_full$MYC_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'PCAT1 cn' = anno_simple(sample_full$PCAT1_Log2_cnv, col = MYC_log2_cnv_col_fun, pch = construct_pch(sample_full$PCAT1_Log2_cnv), na_col = "white", gp = gpar(col = "white"), height = unit(5, "mm")),
      'empty_2' = anno_empty(border = FALSE, height = unit(0.25, "mm")),
      annotation_name_side = "left"
    )
    
    # add legend for annotations
    lgd_patient_id = Legend(title = "Patient ID", labels = names(patient_id_col_fun), legend_gp = gpar(fill = unname(patient_id_col_fun)), border = "black")
    lgd_asi_drug = Legend(title = "Androgen signalling\ninhibitors (ASI)", labels = names(asi_drug_col_fun), legend_gp = gpar(fill = unname(asi_drug_col_fun)), border = "black")
    lgd_ar_score = Legend(title = "AR/NE\nscore", labels = c("Low","Medium","High"), at = c(-7,0,7), col_fun = AR_score_col_fun)
    lgd_site = Legend(title ="Metastasis\nsite", labels = names(site_col_fun), legend_gp = gpar(fill = unname(site_col_fun)), border = "black")
    lgd_cluster = Legend(title ="mCRPC subtypes\n(Labrecque et al)", labels = names(Cons_clust_col_fun), legend_gp = gpar(fill = unname(Cons_clust_col_fun)), border = "black")
    lgd_log2_cnv = Legend(title ="log2(cn)", labels = c("0","1","2",">3.3"), at = c(0, 1, 2, 3.3), col_fun = AR_log2_cnv_col_fun)

    lgd_list <- c(lgd_list, 
                  lgd_patient_id,
                  lgd_asi_drug, 
                  lgd_ar_score,
                  lgd_cluster,
                  lgd_log2_cnv,
                  lgd_site)
    
  }

  ht_mat <- as.matrix(exp_mat)
  if (length(extra_ha) == 0) {
    ht <- Heatmap(ht_mat, name = " z-score\nlog2(TPM +1)", 
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = T, border = TRUE, 
                  rect_gp = gpar(col = "white"),
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 8, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 14),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                  },  top_annotation = clinical_ha , column_order = id_order)
  } 
  else {
    ht <- Heatmap(ht_mat, name = " z-score\nlog2(TPM +1)",
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = T,border = TRUE, 
                  rect_gp = gpar(col = "white"), column_split = sample_full$Cluster_ID, column_title_side = "top", 
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 8, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 14),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                    }, column_order = id_order,  top_annotation = c(extra_ha, clinical_ha) ) 
  }

  
  if (include_oncoprint) {
    
    col_oncoP = c( "NULL"= "white","ecDNA" = "#EB8D71","BFB" = "#009392",
                   "unknown" = "#EBCE8D", "Complex" =  "#9CCB86")
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
      
      'Linear' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9,
                  gp = gpar(fill = col_oncoP["Linear"], col = "black", lwd = 0.5))
      },
      'Complex' = function(x, y, w, h) {
        grid.rect(x, y, w*0.9, h*0.9,
                  gp = gpar(fill = col_oncoP["Complex"], col = "black", lwd = 0.5))
      }
      
    )
    oncoP<- oncoPrint(matrix_oncoP, top_annotation = NULL,remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                      column_order = id_order, border = TRUE, row_order = 1:nrow(matrix_oncoP),
                      alter_fun = alter_fun, col = col_oncoP, pct_side = F, alter_fun_is_vectorized = FALSE, row_gap = unit(1, "mm"),
                      row_names_side = "left", show_column_names = TRUE ,
                      right_annotation = NULL ) 
    
    plot_list =  ht %v% oncoP 
    ComplexHeatmap::draw(plot_list,   annotation_legend_list = lgd_list, merge_legend = T)

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
      ComplexHeatmap::draw(plot_list,   annotation_legend_list = lgd_list, merge_legend = T)

    if (add_line_in_annotation) {
      decorate_annotation(anno_str, {
        grid.lines(c(0, 1), unit(c(line_value, line_value), "native"), gp = gpar(col = "red", lty = 2))
      })
    }
    dev.off()
  }
}


#===============================================================================
# To use the function plot_rna_ecDNA_annotations_ht, I re-ordered Log2.TPM_mod.zc
# matrix_AA_gene_ca_freq_sub_new_mod, RNA_Cascade dataframes based on AR expression

plot_rna_ecDNA_annotations_ht( include_oncoprint = T, exp_mat = Log2.TPM_mod.zc_order[c("AR","MYC","PCAT1"),id_order_ar_cascade],
                               sample_full = RNA_Cascade[match(id_order_ar_cascade, RNA_Cascade$sample_id),], include_clinical_info = T, matrix_oncoP = matrix_AA_gene_ca_freq_sub_new_mod[c("AR_enhancer","AR","MYC","PCAT1"),id_order_ar_cascade],
                               id_order = id_order_ar_cascade, save = T, add_line_in_annotation = F, complex_events = cascade_events_mod[id_order_ar_cascade,3:6],
                               fn_output_plot = "Supplementary_Figure_3A.pdf", width = 13.5, height = 5) #height = 9.7)


