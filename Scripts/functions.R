#,===============================================================================
#' Description:
#' Pre-processing of poppy objects to generate a RData and reproduce the figures
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
#' Functions necessary for the analysis 
#****************************************************************************
################################################################################

#===============================================================================
#'  Function to calculate chromothripsis from David's Cell paper
#===============================================================================
chromo_score = function( sample_id="", chrom="", 
                         window_density=10000, window_width=2500000 ){
  
  total_length=chrom_lengths$V2[ which(rownames(chrom_lengths)==chrom) ]
  half_window=(window_width/2)
  window_starts = seq(from=1, to=total_length+window_width-1, by=window_density)
  interval_start = window_starts - half_window
  interval_end = window_starts+half_window
  
  n_inv = rep(0, length(window_starts))
  n_del = rep(0, length(window_starts))
  n_cna_switch = rep(0, length(window_starts))
  med_CNA = rep(0, length(window_starts)) 
  
  list_sv_m_cur = list_sv_gridss[list_sv_gridss$sample_id==sample_id ,]
  bed = (SO[[sample_id]]$segments )
  bed_chrom = bed[bed$chrom==chrom,]
  
  for( i in 1:length(window_starts) ){
    n_inv[i] = sum( list_sv_m_cur$svtype=="INV" & 
                      list_sv_m_cur$chrom_1==chrom & list_sv_m_cur$chrom_2 == chrom &
                      (
                        (list_sv_m_cur$start_1>interval_start[i] &
                           list_sv_m_cur$start_1<interval_end[i]) | 
                          (list_sv_m_cur$end_2>interval_start[i] &
                             list_sv_m_cur$end_2<interval_end[i])
                      )
    )
    n_del[i] = sum( list_sv_m_cur$svtype=="DEL" & 
                      list_sv_m_cur$chrom_1==chrom & list_sv_m_cur$chrom_2 == chrom &
                      (
                        (list_sv_m_cur$start_1>interval_start[i] &
                           list_sv_m_cur$start_1<interval_end[i]) | 
                          (list_sv_m_cur$end_2>interval_start[i] &
                             list_sv_m_cur$end_2<interval_end[i])
                      )
    )
    n_bounces=0
    bed_idx = which( bed_chrom$start>interval_start[i] &
                       bed_chrom$start<interval_end[i] )
    if( length(bed_idx) > 1 ){
      cnas = bed_chrom$copies[bed_idx]
      cur_dir = "greater"
      for(j in 1:(length(cnas)-1)){
        if( cur_dir=="greater" ){
          if( cnas[j+1] > cnas[j] ){
            n_cna_switch[i] = n_cna_switch[i]+1
            cur_dir = "lesser"            
          }
        }else{
          if( cnas[j+1] < cnas[j] ){
            n_cna_switch[i] = n_cna_switch[i]+1
            cur_dir = "greater"            
          }
        }
      }
    }
    med_CNA[i] = median( (bed_chrom$copies[bed_idx] * (bed_chrom$end[bed_idx]-bed_chrom$start[bed_idx]))/window_width, na.rm=TRUE )
  }
  
  data.frame(sample_id=rep(sample_id, length(window_starts)),
             chrom=rep(chrom, length(window_starts)),
             window_starts, 
             n_inv, 
             n_del, 
             n_cna_switch,
             med_CNA)
}

#===============================================================================
#'  Match.idx and get.split.col functions from David
#===============================================================================

match.idx <-function(A, B, allow.multiple.B=F){
  if( allow.multiple.B ){
    idx.B = which(B %in% A)
    idx.A = match(B[idx.B], A)
  }
  else{
    in.both = intersect(A,B)
    idx.A = match(in.both, A)
    idx.B = match(in.both, B)
  }
  C= data.frame(idx.A, idx.B)
  if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
    stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
  C
}


get.split.col = function(v, string, col=0, last=F, first=F){
  if( last & first )
    stop("Cannot request both last and first column")
  if( col==0 & !last & !first)
    stop("Must request either a column by index, first, or last")
  for(i in 1:length(v)){
    x = strsplit( v[i], string, fixed=T)[[1]]
    if(last){
      v[i] = x[length(x)]
    }
    else if(first){
      v[i] = x[1]
    }
    else{
      v[i] = x[col]
    }
  }
  v
}

#===============================================================================
#'  Function to annotate the copy number from a non-coding region (e.g AR enhancer)
#===============================================================================

somatic_add_CNA_by_region_from_segments = function( somatic, regions ){
  
  N = dim(regions)[1]
  CN_chrom = rep(NA, N )
  CN_weighted_mean = rep(NA, N )
  CN_largest_segment = rep(NA, N )
  gr_segments <- GenomicRanges::makeGRangesFromDataFrame(somatic$segments[,c("chrom","start","end","copyNumber")],
                                                         keep.extra.columns = T)
  # identify segments that intersect the start and end of each gene
  # trim intersecting segments to gene bounds
  # calculate 1) weighted average segment CN and 2) CN of largest segment
  for(i in 1:N ){
    gr_region = GenomicRanges::GRanges(seqnames=regions$chrom[i],
                                       ranges = IRanges(start=regions$start[i],
                                                        end=regions$end[i]))
    df_sv_overlapping_gene = data.frame( subsetByOverlaps( gr_segments, gr_region) )
    df_sv_overlapping_gene$start[ df_sv_overlapping_gene$start < regions$start[i] ] = regions$start[i]
    df_sv_overlapping_gene$end[ df_sv_overlapping_gene$end > regions$end[i] ] = regions$end[i]
    CN_weighted_mean[i] = weighted.mean( df_sv_overlapping_gene$copyNumber, df_sv_overlapping_gene$width)
    CN_largest_segment[i] = df_sv_overlapping_gene$copyNumber[ which( df_sv_overlapping_gene$width == max(df_sv_overlapping_gene$width)[1] ) ]
    CN_chrom[i] = regions$chrom[i]
  }
  df_regions=data.frame(chrom = CN_chrom,
                        weighted_mean=CN_weighted_mean, 
                        largest_segment = CN_largest_segment,
                        region=regions$subset.CNA )
  df_regions$weighted_mean = round( df_regions$weighted_mean, 3)
  df_regions$largest_segment = round( df_regions$largest_segment, 3)
  df_regions
}
#===============================================================================
#'  Function to generate a matrix of coding mutations
#===============================================================================

generate_coding_matrix_mod = function( idx,min_depth=40 ){
  symbols_gene = c()
  for(i in idx){
    me = SO[[i]]$mutect_exome
    idx_effect=grep( "frameshift|stop_gain|splice_acceptor|splice_donor|missense", me$effect)
    idx_depth = which(me$DP_SNV_alt_T >= min_depth)
    
    symbols_gene = c(symbols_gene, unique(
      paste(me[ intersect( idx_effect, idx_depth ) ,]$gene,
            me[ intersect( idx_effect, idx_depth ) ,]$pos,
            me[ intersect( idx_effect, idx_depth ) ,]$ref,
            me[ intersect( idx_effect, idx_depth ) ,]$alt, sep="_")
    ))
  }
  symbols_gene = sort( unique(symbols_gene) )
  M_var = matrix( 0, nrow=length(symbols_gene), ncol=length(idx))
  dimnames(M_var)[[1]] = symbols_gene
  dimnames(M_var)[[2]] = get.split.col( names(SO)[ idx ], "_", first=TRUE)
  
  mut_phyl = matrix(nrow=length(idx) , ncol = 2)

  for(i in 1:length(idx)){
    me = SO[[ idx[i] ]]$mutect_exome
    idx_effect=grep( "frameshift|stop_gain|splice_acceptor|splice_donor|missense", me$effect)
    idx_depth = which(me$DP_SNV_alt_T>= min_depth)
    symbols_sample = paste(me[ intersect( idx_effect, idx_depth ) ,]$gene,
                           me[ intersect( idx_effect, idx_depth ) ,]$pos, 
                           me[ intersect( idx_effect, idx_depth ) ,]$ref, 
                           me[ intersect( idx_effect, idx_depth ) ,]$alt, sep="_")
    
    m = match.idx( dimnames(M_var)[[1]], symbols_sample )
    M_var[m$idx.A, i] = 1
    mut_app <- c()
    for (j in 1:length(symbols_gene) ){
      
      if(  symbols_gene[j]  %in% symbols_sample  ) {
        mut_app = append(mut_app, 1)
      }
      else{
        mut_app = append(mut_app, 0)
        next
      }
      mut_app
    }
    mut_phyl[i,1] <-do.call('rbind',strsplit(as.character(SO[[ idx[i] ]]$sample_base), '-' ,fixed=TRUE))[,2]
    mut_phyl[i,2] <- paste0(mut_app,collapse="")
  }
  M_var = list( M_var, mut_phyl )
  M_var
}

#===============================================================================
#'  Function to annotate enhancer amplified by complex events  (including ecDNA, BFB, etc)
#'  e.g.  AR enhancer 
#===============================================================================

annotate_enhancer_interval <- function(  coord_intervals, anno_region ){
  
  region_GR <- GenomicRanges::makeGRangesFromDataFrame(anno_region,keep.extra.columns = T)
  int_gr <- GenomicRanges::makeGRangesFromDataFrame(coord_intervals,keep.extra.columns = T)
  hits <- GenomicRanges::findOverlaps( region_GR, int_gr )
  int_gr_sub <- int_gr[subjectHits(hits)]
  return(int_gr_sub) 
}

#===============================================================================
#'  Other functions
#===============================================================================

# The oposite of "%in%" function
# Taken from https://stackoverflow.com/questions/38351820/negation-of-in-in-r

'%out%' <- function(a,b) ! a %in% b

# Get upper and lower  triangle of the correlation matrix
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
# https://stackoverflow.com/questions/47475897/correlation-matrix-tidyr-gather-v-reshape2-melt

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


#===============================================================================
#'  My colors palette collection
#===============================================================================

my_hue3_30 <- c("#ee6e96", "#e38d95", "#ffab96", "#eb825c", "#ff9f54", "#d3a474",
                         "#d29e3f", "#ffd26e", "#f3eb66", "#e3dd95", "#9ea838", "#9cad5b",
                         "#93ad6d", "#bce2ab", "#82b07f", "#84e796", "#719a7e", "#55ce97",
                         "#4fae8e", "#43f2ca", "#7de7d3", "#5bbbb9", "#55d1e9", "#59b0d6",
                         "#4fa2e8", "#a9bdf1", "#8c97d3", "#928fee", "#d9afea", "#cf91bd")
                         

# SV_colores<- c("#E41A1C","#f2855d","#377EB8", "#4DAF4A", "#984EA3", "#FF7F00") #old colors

SV_colores<- c("#0072B2","#56B4E9","#E69F00","#009E73","#984EA3", "#E41A1C")

# Barrier-free color palette
# Source: Okabe & Ito (2008): Color Universal Design (CUD):
#         Fig. 16 of <https://jfly.uni-koeln.de/color/>:
# https://bookdown.org/hneth/ds4psy/D-5-apx-colors-define-use-custom.html#apx:defining-color-palettes
# (a) Vector of colors (as RGB values):
o_i_colors <- c(rgb(  0,   0,   0, maxColorValue = 255),  # black
                rgb(230, 159,   0, maxColorValue = 255),  # orange
                rgb( 86, 180, 233, maxColorValue = 255),  # skyblue
                rgb(  0, 158, 115, maxColorValue = 255),  # green
                rgb(240, 228,  66, maxColorValue = 255),  # yellow
                rgb(  0, 114, 178, maxColorValue = 255),  # blue
                rgb(213,  94,   0, maxColorValue = 255),  # vermillion
                rgb(204, 121, 167, maxColorValue = 255)   # purple
)

# (b) Vector of color names:
o_i_names <- c("black", "orange", "skyblue", "green", "yellow", "blue", "vermillion", "purple")

library(unikn)
# (c) Use newpal() to combine colors and names:
pal_okabe_ito <- newpal(col = o_i_colors,
                        names = o_i_names)

