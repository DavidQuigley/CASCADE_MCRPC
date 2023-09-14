#,===============================================================================
#' Description:
#' Reproduce the Figure 4A
#' Adapted from David's script to plot frequency of reversion mutations
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Authors: 
#' - David Quigley & Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Uses aardvark output to plot Figure 1B
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
#' Figure 1B
#===============================================================================
fn_heatmap = "Figure_1B.pdf"

#===============================================================================
# heatmap of reversion mutation frequency

MIN_OBS=3
hm = data.frame()
for(i in 1:length(sums_all)){
    idx = which( sums_all[[i]]$summary$N >= MIN_OBS &
                     sums_all[[i]]$summary$evidence != "overlaps_blacklist")
    h = data.frame( sample = rep( names(sums_all)[i], length(idx) ),
                    slug = dimnames( sums_all[[i]]$summary )[[1]][ idx ],
                    N = sums_all[[i]]$summary$N[idx],
                    stringsAsFactors = FALSE)
    hm = rbind(hm, h)
}

hm$sample = factor( hm$sample, levels = names(sums_all)[ length(names(sums_all)):1] )

rev_share_filter = rev_share
rev_share_filter[rev_share_filter<3]=0
rev_share_filter=rev_share_filter[which(rowSums(rev_share_filter>0)>1),]
h=hclust( dist( t(rev_share_filter) ) )
#plot(h)

h_rev=hclust( dist( rev_share_filter ) )

n_seen = rep(0, dim(hm)[1])
for(i in 1:length(n_seen)){
    n_seen[i] = sum( hm$slug == hm$slug[i] )
}

#===============================================================================
# test whether samples observed in a single tumor have different frequency from
# those seen in multiple tumors
name_singleton = dimnames(rev_share)[[1]][ rowSums( rev_share>0 )==1  ]
m = match.idx( name_singleton, hm$slug)
# boxplot( log2( hm$N[ m$idx.B ]), log2( hm$N[ setdiff( (1:dim(hm)[1]), m$idx.B ) ] ) )
wilcox.test( hm$N[ m$idx.B ], hm$N[ setdiff( (1:dim(hm)[1]), m$idx.B ) ] )

#W = 3205.5, p-value = 0.0008186

t.test( hm$N[ m$idx.B ], hm$N[ setdiff( (1:dim(hm)[1]), m$idx.B ) ] )

#===============================================================================
# creating levels

hm$sample = factor( hm$sample, levels = dimnames(rev_share_filter)[[2]][ h$order ] )
hm$slug = factor( hm$slug, levels = names(sort( table(hm$slug), decreasing = TRUE)) )

# trim slugs
slugs = as.character(hm$slug )
for(i in 1:length(slugs) ){
    newslug = c()
    for( token in strsplit( slugs[i], "_")[[1]]  ){
        newslug = c(newslug, paste( strsplit( token, ":")[[1]][1:3], sep=":", collapse=":" ) )
    }
    if( length(newslug)==1 ){
        slugs[i] = newslug[1]
    }else{
        slugs[i] = paste( newslug, collapse=" ")
    }
}
hm$slug = slugs
hm$slug = factor( hm$slug, levels = names( sort( table(hm$slug), decreasing = TRUE)) )

#===============================================================================
# Plot

pdf( fn_heatmap, height=6.5, width=12)

(p <- ggplot( hm, aes(slug, sample) ) +
        labs(x = "", y = "") +
        geom_tile(aes(fill = N), colour = "white") +
        theme(panel.background = element_rect(fill = "#eeeeee", size=1, colour = "black"),
              axis.ticks.y = element_blank(), axis.ticks.x = element_blank() ,
              axis.text.y = element_text( size=14),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8),
              panel.border = element_rect(colour = "black", fill=NA, size=0.75) ) +
        scale_fill_gradient(name="Number of\nreads", low = "lightblue1", high = "black",
                            trans="log10") )

dev.off()


#===============================================================================
# Total number of individual reversions
dim(rev_share)[1]
# 129

# number of reversions shared in multiple samples
sum( rowSums( rev_share>=MIN_OBS) > 1 )
# 39

# number of reversions shared in at least 3 samples
sum( rowSums( rev_share>=MIN_OBS) >= 3)
# 20

# most widely shared reversion
max( rowSums( rev_share>=MIN_OBS)  )
# 6
which(  rowSums( rev_share>=MIN_OBS) ==6 )
# D:30:32339629:32339658 D:51:32339610:32339660

#===============================================================================
# Write out supplementary table
tbl = rev_share
tbl[tbl<MIN_OBS] = 0
which(rowSums(tbl>0)==0)
write.table(tbl, paste0(dir_cascade_results,'supplementary_table_reversions.txt'), sep='\t', quote=FALSE)
