#,===============================================================================
#' Description:
#' Reproduce the Figure 1A
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
#' Uses ComplexHeatmap and a set of matrices to plot Figure 1A
#****************************************************************************
################################################################################

#===============================================================================
#'  Function modified from Meng's script
#===============================================================================
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("circlize")

library( ggplot2 )
library( ComplexHeatmap )
library( circlize )
library( RColorBrewer )
library( grid )
require( gridExtra )
library( Sushi )
library( Gviz )
library( GenomicFeatures )
library( reshape2 )
library( ggrepel )
library( rtracklayer )
require( methods )
library( stringr )
library( rcartocolor )
library( methods )



# This function plots CASCADE sample annotation with genomic and clinical
# information inlcuding HR-deficiency status. Make sure that you have used the script
# Prepare_cascade_eviroment_2023.R to load and generate all the data.

plot_sample_annotation_sig_ht <- function(status_to_include = F, include_fusion, include_ar_enhancer, include_oncoprint = TRUE,
                                      sample_full, include_clinical_info = TRUE, matrix_oncoP, 
                                      id_order, sv_matrix, extra_ha = c(), extra_lgd = c(), add_line_in_annotation = F, anno_str, line_value,
                                      save = T, fn_output_plot, width = 15, height = 10) {

  lgd_list <- list()
  clinical_ha <- c()
  if (include_clinical_info) {
    # annotate samples
    patient_id_col_fun <- c("CA071" = "#0D0887FF", "CA093" = "#6A00A8FF", "CA101" = "#B12A90FF", "CA104" = "#E16462FF", "CA108" = "#FCA636FF", "CA113" = "#F0F921FF")
    hrd_drug_col_fun <- c("none" = "white", "platinum" = "#D9D9D9", "PARPi" = "#737373", "both" = "#252525")
    wgd_col_fun <- c("FALSE" = "white", "TRUE" = "black")
    purity_col_fun <- colorRamp2(c(0, 100), c("white", "#3F007D"))
    tmbPerMblog2_col_fun <- colorRamp2(quantile(sample_full$tmbPerMb_Log2), brewer.pal(n = 5, name = "Greys"))
    ets_col_fun <- c("FALSE" = "white", "TRUE" = "black")
    HRD_score_col_fun <- colorRamp2(quantile(sample_full$HRD_p_hrd), brewer.pal(n = 5, name = "Greys"))
    chromo_col_fun <- c("FALSE" = "white", "TRUE" = "black")
    site_col_fun <- c(rcartocolor::carto_pal (8, "Prism")); names(site_col_fun) <- c("Bone", "Prostate", "Spine", "LN", "Brain", "Liver", "Lung", "Others") 
    
    construct_pch <- function(l) {
      pch <- rep(NA, length(l))
      is_na <- is.na(l)
      pch[is_na] <- ""
      return(pch)
    }
    align_to = list("BRCA2 Germline" = which(sample_full$HR_status == "BRCA2 Germline"), "BRCA2 Somatic"= which(sample_full$HR_status == "BRCA2 Somatic"),
                    "ATM"= which(sample_full$HR_status == "ATM"), "HR Proficient" = which(sample_full$HR_status == "HR Proficient" ))
    panel_fun = function(index, nm) {
        grid.rect()
        grid.text()
    }
    clinical_ha = HeatmapAnnotation(
      'Patient ID' = anno_simple(sample_full$patient_id, col = patient_id_col_fun, pch = construct_pch(sample_full$patient_id), na_col = "white", border = TRUE, height = unit(0.5, "cm")),
      'Purity' = anno_barplot(sample_full$purity_PURPLE, height = unit(1.5, "cm"), gp = gpar(col = "black", fill = "gray")),
      'Ploidy' = anno_barplot(sample_full$ploidy_PURPLE, height = unit(1.5, "cm"), gp = gpar(col = "black", fill = "gray")),
      'Structural\nVariants\nCounts' = anno_barplot(sv_matrix, height = unit(5, "cm"), gp = gpar(col = "black", fill = SV_colores)),
      'Single Base\nSubtistutions\nSignatures' = anno_barplot(sample_full[, str_detect(colnames(sample_full), "Sig_")], height = unit(5, "cm"), gp = gpar(col = "black", fill = my_hue3_30[c(1,3,4,6,8,15,18,24,25,26,29)])),
      'Metastasis site'= anno_simple(sample_full$Tissue, col = site_col_fun, pch = construct_pch(sample_full$Tissue), na_col = "white", gp = gpar(col = "white")),
      "Platinum/PARPi" = anno_simple(sample_full$HRD_drug, col = hrd_drug_col_fun, pch = construct_pch(sample_full$HRD_drug), na_col = "white",gp = gpar(col = "white")),
      'log2(tmbPerMb)' = anno_simple(sample_full$tmbPerMb_Log2, col = tmbPerMblog2_col_fun, pch = construct_pch(sample_full$tmbPerMb_Log2), na_col = "white", gp = gpar(col = "white")),
      'WGD' = anno_simple(sample_full$WG_doubling, col = wgd_col_fun, pch = construct_pch(sample_full$WG_doubling), na_col = "white", gp = gpar(col = "white")),
      'ETS rearrangements' = anno_simple(sample_full$ETS, col = ets_col_fun, pch = construct_pch(sample_full$ETS), na_col = "white", gp = gpar(col = "white")),
      'Chromothripsis' = anno_simple(sample_full$Chromothripsis, col = chromo_col_fun, pch = construct_pch(sample_full$Chromothripsis), na_col = "white", gp = gpar(col = "white")),
      'CHORD score' = anno_simple(sample_full$HRD_p_hrd, col = HRD_score_col_fun, pch = construct_pch(sample_full$HRD_p_hrd), na_col = "white", gp = gpar(col = "white")),
           annotation_name_side = "left"
    )
    # add legend for annotations
    lgd_patient_id = Legend(title = "Patient ID", labels = names(patient_id_col_fun), legend_gp = gpar(fill = unname(patient_id_col_fun)), border = "black")
    lgd_hrd_drug = Legend(title = "Platinum/PARPi", labels = names(hrd_drug_col_fun), legend_gp = gpar(fill = unname(hrd_drug_col_fun)), border = "black")
    lgd_ets = Legend(title = "ETS\nrearrangement", labels = names(ets_col_fun), legend_gp = gpar(fill = unname(ets_col_fun)), border = "black")
    lgd_wgd = Legend(title = "Whole Genome\nDoubling (WGD)", labels = names(wgd_col_fun), legend_gp = gpar(fill = unname(wgd_col_fun)), border = "black")
    lgd_tmb = Legend(title = "log2(tmbPerMb)", labels = c("-6", "0","2", "4", "8"), at = c(-6, 0, 2, 4, 8), col_fun = tmbPerMblog2_col_fun)
    lgd_hrd_score = Legend(title = "CHORD\nscore", labels = c("0","0.6", "0.7", "0.8","1"), at = c(0, 0.6, 0.7, 0.8, 1), col_fun = HRD_score_col_fun)
    lgd_chromo = Legend(title = "Chromothripsis", labels = names(chromo_col_fun), legend_gp = gpar(fill = unname(chromo_col_fun)), border = "black")
    lgd_site = Legend(title ="Metastasis\nsite", labels = names(site_col_fun), legend_gp = gpar(fill = unname(site_col_fun)), border = "black")
    lgd_sv = Legend(title = "Structural\nVariants", labels = c(colnames(sv_matrix)), legend_gp = gpar(fill = SV_colores), border = "black")
    lgd_sig = Legend(title = "Signatures\nCosmic v2", labels = c(colnames(sample_full[, str_detect(colnames(sample_full), "Sig_")])), legend_gp = gpar(fill = my_hue3_30[c(1,3,4,6,8,15,18,24,25,26,29)]), border = "black")
    
    lgd_list <- c(lgd_list, 
                  lgd_patient_id,
                  lgd_sv,
                  lgd_sig,
                  lgd_site,
                  lgd_tmb,
                  lgd_wgd,
                  lgd_ets,
                  lgd_chromo, 
                  lgd_hrd_score,
                  lgd_hrd_drug)
  }
  
  
  if (status_to_include) {
  # plot extra status as a heatmap
  for (status in status_to_include) {
    print(status)
    
    
    sample_full[, status] <- paste(sample_full[, status], sample_full[, paste0(status, "_status")]
                                 )
  }
  status_to_include <- c()
  ht_mat <- as.matrix(sample_full[,status_to_include])
  rownames(ht_mat) <- sample_full$sample_id
  ht_mat[is.na(ht_mat)] = "NA"
  color_names <- unique(c(ht_mat))
  status_col_fun <- rep(NA, length(color_names))
  names(status_col_fun) <- color_names
  status_col_fun[str_detect(names(status_col_fun), "NA|2")] <- "white"
  status_col_fun[str_detect(names(status_col_fun), "0")] <- "blue"
  status_col_fun[str_detect(names(status_col_fun), "1")] <- "dark gray"
  status_col_fun[str_detect(names(status_col_fun), "3")] <- "red"
  }
  else{
      status_to_include <- c()
      ht_mat <- as.matrix(sample_full[,status_to_include])
      rownames(ht_mat) <- sample_full$sample_id
      ht_mat[is.na(ht_mat)] = "NA"
  }
  
  if (length(extra_ha) == 0) {
    ht <- Heatmap(t(ht_mat),
                  col = status_col_fun, row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"),
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F, 
                  rect_gp = gpar(col = "black"), column_split = sample_full$HR_status, column_title_side = "top",
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 10, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 12),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(t(ht_mat)[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) # "."
                    }
                  },  top_annotation = clinical_ha , column_order = id_order)
  } else {
    ht <- Heatmap(t(ht_mat),
                  col = status_col_fun, row_names_side = "left", show_row_names = T,column_gap = unit(1.5, "mm"),
                  cluster_rows = F, show_column_dend = FALSE, show_heatmap_legend = F, 
                  rect_gp = gpar(col = "black"), column_split = sample_full$HR_status, column_title_side = "top",
                  column_title_gp = gpar(fill = c("white"), just = "bottom",fontsize = 10, fontface = "bold"),show_column_names = T,column_title_rot = 0,
                  column_names_gp = gpar(col = c("black"), fontsize = 12),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(str_detect(t(ht_mat)[i, j], "TRUE")) {
                      grid.points(x, y, pch = 16, size = unit(2.5, "mm")) 
                    }
                    }, column_order = id_order,  top_annotation = c(extra_ha, clinical_ha)) 
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
                     "Multi-hit" = "#7E7E7E" )
      
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
      'Missense' = function(x, y, w, h) {
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
                      border = TRUE, 
                      alter_fun = alter_fun, col = col_oncoP, pct_side = F,#"right",
                      row_names_side = "left", show_column_names = TRUE,
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
# To use the function plot_sample_annotation_sig_ht, SV_gridss_summary, SO_summary, 
# and mat_onco matrices are necessary!


plot_sample_annotation_sig_ht(status_to_include = F, sv_matrix = SV_gridss_summary,
                     include_fusion = F, include_ar_enhancer = F, include_oncoprint = T,
                     sample_full = SO_summary, include_clinical_info = T, matrix_oncoP = mat_onco[c("BRCA2", "BRCA1", "ATM", "CDK12", "MSH2" , "MSH6", "TP53", "AR", "AR_enhancer","RB1", "FOXA1", "PTEN", "CTNNB1"),],
                     id_order = SO_summary$sample_id, save = T, add_line_in_annotation = F, 
                     fn_output_plot = "Figure_1a.pdf", width = 17, height = 10.5) 

