#,===============================================================================
#' Description:
#' Reproduce the Figure 2A
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
#' Uses ComplexHeatmap and a set of matrices to plot Figure 2A
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. Function adapted from Figure 1A, and Meng's script.
#===============================================================================
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("circlize")

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
library(rcartocolor)

#===============================================================================
#'  Function modified to plot clinical anotation, RNA expression as well as 
#'  genomic alterations in cascade cohort.
#===============================================================================

plot_rna_annotations_ht <- function(include_oncoprint = TRUE, add_line_in_annotation = F, anno_str, line_value,
                                      sample_full, include_clinical_info = TRUE, matrix_oncoP, exp_mat, 
                                      id_order, sv_matrix, mut_sig, extra_ha = c(), extra_lgd = c(),
                                      save = T, fn_output_plot, width = 15, height = 10) {

  library(methods)
  lgd_list <- list()
  clinical_ha <- c()
  if (include_clinical_info) {
    patient_id_col_fun <- c("CA071" = "#0D0887FF", "CA093" = "#6A00A8FF", "CA101" = "#B12A90FF", "CA104" = "#E16462FF", "CA108" = "#FCA636FF", "CA113" = "#F0F921FF")
    site_col_fun <- c(rcartocolor::carto_pal (8, "Prism"));names(site_col_fun) <- c("Bone", "Prostate", "Spine", "LN", "Brain", "Liver", "Lung", "Others") 
    asi_drug_col_fun <- c("none" = "white", "abi" = "#D9D9D9", "enz" = "#737373", "both" = "#252525")
    Cons_clust_col_fun <- c( pal_okabe_ito[c(-1, -5, -8)]);names(Cons_clust_col_fun)<- c('AR+/NE+','AR+/NE-','ARL/NE-','AR-/NE+','AR-/NE-')
    AR_score_col_fun <- colorRamp2(c(quantile(sample_full$AR_score, 0), quantile(sample_full$AR_score, 0.5), quantile(sample_full$AR_score, 1)), c("white","#D9D9D9","black"))#colorRamp2(c(quantile(sample_full$AR_score, 0), quantile(sample_full$AR_score, 0.5), quantile(sample_full$AR_score, 1)), c("green","white", "red"))
    NE_score_col_fun <- colorRamp2(c(quantile(sample_full$NE_score, 0), quantile(sample_full$NE_score, 0.5), quantile(sample_full$NE_score, 1)), c("white","#D9D9D9","black"))#brewer.pal(n = 3, name = "Greys"))
    construct_pch <- function(l) {
      pch <- rep(NA, length(l))
      is_na <- is.na(l)
      pch[is_na] <- ""
      return(pch)
    }

    clinical_ha = HeatmapAnnotation(
      'Patient ID' = anno_simple(sample_full$patient_id, col = patient_id_col_fun, pch = construct_pch(sample_full$patient_id), na_col = "white", border = TRUE),
      'Metastasis site'= anno_simple(sample_full$Tissue, col = site_col_fun, pch = construct_pch(sample_full$Tissue), na_col = "white", gp = gpar(col = "white")),
      "ASI" = anno_simple(sample_full$ASI_drug, col = asi_drug_col_fun, pch = construct_pch(sample_full$ASI_drug), na_col = "white",gp = gpar(col = "white")),
      'AR score' = anno_simple(sample_full$AR_score, col = AR_score_col_fun, pch = construct_pch(sample_full$AR_score), na_col = "white", gp = gpar(col = "white")),
      'NE score' = anno_simple(sample_full$NE_score, col = NE_score_col_fun, pch = construct_pch(sample_full$NE_score), na_col = "white", gp = gpar(col = "white")),
      'mCRPC subtypes'= anno_simple(sample_full$Cluster_ID, col = Cons_clust_col_fun, pch = construct_pch(sample_full$Cluster_ID), na_col = "white", gp = gpar(col = "white")),
           annotation_name_side = "left"
    )
    
    # add legend for annotations
    lgd_patient_id = Legend(title = "Patient ID", labels = names(patient_id_col_fun), legend_gp = gpar(fill = unname(patient_id_col_fun)), border = "black")
    lgd_asi_drug = Legend(title = "Androgen signalling\ninhibitor (ASI)", labels = names(asi_drug_col_fun), legend_gp = gpar(fill = unname(asi_drug_col_fun)), border = "black")
    lgd_scores = Legend(title = "AR/NE scores", labels = c("Low","Medium","High"), at = c(-7,0,7), col_fun = AR_score_col_fun)
    lgd_site = Legend(title ="Metastasis\nsite", labels = names(site_col_fun), legend_gp = gpar(fill = unname(site_col_fun)), border = "black")
    lgd_cluster = Legend(title ="mCRPC subtypes\n(Labrecque et al)", labels = names(Cons_clust_col_fun), legend_gp = gpar(fill = unname(Cons_clust_col_fun)), border = "black")
    
    lgd_list <- c(lgd_list, 
                  lgd_patient_id,
                  lgd_asi_drug, 
                  lgd_scores,
                  lgd_cluster,
                  lgd_site)
  }
  
  ht_mat <- as.matrix(exp_mat)
  idx_keep = which( P.Nelson_Sig$Gene_ID %in% rownames(ht_mat) )
  col.Nelson <- c(rcartocolor::carto_pal(4, "Temps"))
  names(col.Nelson) <- unique(P.Nelson_Sig[idx_keep,]$GO)
  
  rowAnnot.Nelson = rowAnnotation(
      GO = P.Nelson_Sig[idx_keep,]$GO,
      col = list(GO=col.Nelson)
  )
  
  if (length(extra_ha) == 0) {
    ht <- Heatmap(ht_mat, name = "z-score\nlog2(TPM +1)",
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = T, left_annotation = rowAnnot.Nelson, border = TRUE,
                  rect_gp = gpar(col = "white"), column_split = sample_full$patient_id, row_split  = P.Nelson_Sig[idx_keep,]$GO, 
                  column_names_gp = gpar(col = c("black"), fontsize = 10), column_title = NULL,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                  },  top_annotation = clinical_ha , column_order = id_order)
  } 
  else {
    ht <- Heatmap(ht_mat, name = "z-score\nlog2(TPM +1)",
                  row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"), 
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = T, left_annotation = rowAnnot.Nelson, border = TRUE,
                  rect_gp = gpar(col = "white"), column_split = sample_full$Cluster_ID,  row_split  = P.Nelson_Sig[idx_keep,]$GO, 
                  column_names_gp = gpar(col = c("black"), fontsize = 10), column_title = NULL,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(ht_mat[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                    }, column_order = id_order,  top_annotation = c(extra_ha, clinical_ha) ) 
  }

  if (include_oncoprint) {
    col_oncoP = c( "NULL"= "#D3D3D3",
                   "Copy loss" = "#3288BD",
                   "Copy gain" = "#9E0142",
                   "Stop lost" = "#FDAE61",
                   "Missense" = "#66C2A5",
                   "Missense VUS" ="#66C2A5",
                   "Missense hotspot" = "#3d7463",
                   "Stop gained" = "#F46D43", 
                   "Stop gained hotspot"="#D53E4F", 
                   "Splice region" = "#fa9fb5",
                   "Biallelic"= "black",
                   "Germline" = "#FEE08B",
                   "Frameshift" = "#c0e5bb",
                   "Structural variant" = "#5E4FA2",
                   "Multi-hit" = "#7E7E7E"
    )
    alter_fun = list(
      background = function(x, y, w, h) 
        grid.rect(x, y, w, h, 
                  gp = gpar(fill = "#f0f0f0",
                            col = "white", lwd = 1)),
      'Copy gain' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
          unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
          gp = gpar(fill = col_oncoP["Copy gain"], col = "black"))
      },
      'Copy loss' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
          unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
          gp = gpar(fill = col_oncoP["Copy loss"], col = "black"))
      },
      'Missense VUS' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Missense VUS"], col = "black"))
      },
      'Missense hotspot' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Missense hotspot"], col = "black"))
      }, 
      Missense = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Missense"], col = "black"))
      },
      Frameshift = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Frameshift"], col = "black"))
      },
      'Stop lost' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Stop lost"], col = "black"))
      },       
      'Stop gained hotspot' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Stop gained hotspot"], col = "black"))
      },    
      'Stop gained' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Stop gained"], col = "black"))
      },
      'Splice region' = function(x, y, w, h) {
        grid.polygon(
          unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
          unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
          gp = gpar(fill = col_oncoP["Splice region"], col = "black"))
      },
      'Structural variant' = function(x, y, w, h) 
        grid.rect(x, y, w*0.58, h*0.58, 
                  gp = gpar(fill = col_oncoP["Structural variant"], col = 'black',lwd = 0.6)),
      
      'Multi-hit' = function(x, y, w, h)
        grid.rect(x, y, w*0.38, h*0.38, 
                  gp = gpar(fill = col_oncoP["Multi-hit"], col = 'darkgray',lwd = 0.5)),
      
      Biallelic = function(x, y, w, h) 
        grid.points(x, y, pch = 20, gp=gpar(fill=col_oncoP["Biallelic"],col='black',lwd = 0.8)),
      
      Germline = function(x, y, w, h) 
        grid.points(x, y, pch = 23, gp=gpar(fill=col_oncoP["Germline"],col='black',lwd = 0.8, cex = 0.75))
    )
    
    oncoP<- oncoPrint(matrix_oncoP, top_annotation = NULL,  column_order = id_order, 
                      border = TRUE, row_order = 1:nrow(matrix_oncoP), 
                      alter_fun = alter_fun, col = col_oncoP, pct_side = F,
                      row_names_side = "left", show_column_names = TRUE ,
                      right_annotation = NULL 
                      )
    
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
# To use the function plot_rna_annotations_ht, Log2.TPM_mod.zc.nelson_ft, RNA_Cascade, 
# and mat_onco (same as Fig 1A) matrices are necessary!

plot_rna_annotations_ht( include_oncoprint = T, exp_mat = Log2.TPM_mod.zc.nelson_ft,
                        sample_full = RNA_Cascade, include_clinical_info = T, matrix_oncoP = mat_onco[c("AR", "AR_enhancer", "FOXA1", "CTNNB1", "BRCA2", "BRCA1", "ATM", "CDK12", "MSH2" , "MSH6", "TP53", "RB1","PTEN"),RNA_Cascade$sample_id],
                        id_order = RNA_Cascade$sample_id, save = T, add_line_in_annotation = F, 
                        fn_output_plot = "Figure_2A.pdf", width = 13.5, height = 10.5)


