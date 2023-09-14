#,===============================================================================
#' Description:
#' Reproduce Supplementary Figures 1 A-F
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
#' Uses phangorn and generate_coding_matrix_mod function to create the 
#' phylogenetic trees and plot Supplementary Figure 1 and the individual trees
#****************************************************************************
################################################################################
# We used the Neighbor-Joining (NJ) algorithm because allows for unequal rates 
# of evolution, so that branch lengths are proportional to amount of change. 
# If rates on different branches are not markedly unequal, the branching orders 
# produced by others methods will not differ.

#===============================================================================
#'  Loading libraries. 
#===============================================================================

library( phangorn )
library( ape )
library( readr )
library( ggtree )

#===============================================================================
#' To Draw rotated plot
#===============================================================================

vp_phylo <- viewport(width = 0.5,
                     height = 0.6,
                     x = unit(0.5, "npc"), 
                     y = unit(0.5, "npc"),
                     angle = 270)

#===============================================================================
#' Supplementary Figure 1A CA071
#===============================================================================

idx=grep("CA071", names(SO))
M_var_CA071 = generate_coding_matrix_mod( idx, min_depth = 10 )

CA071_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA071[[1]])[1] ), M_var_CA071[[2]] ) )
readr::write_delim(CA071_mut_phyl, "CA071_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA071_mut_phyl <- read.phyDat("CA071_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA071_mut_phyl  <- dist.ml(CA071_mut_phyl)
CA071_treeNJ  <- nj(dm_CA071_mut_phyl)
CA071_treeNJ_mp <- midpoint(CA071_treeNJ)

pdf("Suppl_Figure_1A_CA071_treeNJ.pdf", width = 4, height = 3)
p_CA071 <- ggtree(CA071_treeNJ,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.5, hjust=0.5, offset = 0.005) 
print(p_CA071,vp=vp_phylo)
dev.off()

#===============================================================================
#' Supplementary Figure 1B CA104
#===============================================================================

idx=grep("CA104", names(SO))
M_var_CA104 = generate_coding_matrix_mod( idx, min_depth = 15 )

CA104_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA104[[1]])[1] ), M_var_CA104[[2]] ) )
readr::write_delim(CA104_mut_phyl, "CA104_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA104_mut_phyl <- read.phyDat("CA104_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA104_mut_phyl  <- dist.ml(CA104_mut_phyl)
# Neighbor joining tree 
CA104_treeNJ  <- nj(dm_CA104_mut_phyl)
CA104_treeNJ_mp <- midpoint(CA104_treeNJ)

pdf("Suppl_Figure_1B_CA104_treeNJ.pdf", width = 4, height = 3)
p_CA104 <- ggtree(CA104_treeNJ_mp,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.5, hjust=0.5, offset = 0.005) 
print(p_CA104,vp=vp_phylo)
dev.off()

#===============================================================================
#' Supplementary Figure 1C CA113
#===============================================================================

idx=grep("CA113", names(SO))
M_var_CA113 = generate_coding_matrix_mod(idx, min_depth=5)

CA113_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA113[[1]])[1] ), M_var_CA113[[2]] ) )
write_delim(CA113_mut_phyl, "CA113_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA113_mut_phyl <- read.phyDat("CA113_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA113_mut_phyl  <- dist.ml(CA113_mut_phyl)
dm_CA113_mut_phyl[is.na(dm_CA113_mut_phyl)] = 0.89

CA113_treeNJ  <- nj(dm_CA113_mut_phyl)
CA113_treeNJ_mp <- midpoint(CA113_treeNJ)

pdf("Suppl_Figure_1B_CA113_treeNJ.pdf", width = 4, height = 3)
p_CA113 <- ggtree(CA113_treeNJ,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.5, hjust=0.5, 
                                              offset = 0.005, ignore.negative.edge=TRUE ) 
print(p_CA113,vp=vp_phylo)
dev.off()

#===============================================================================
#' Supplementary Figure 1D CA101
#===============================================================================
idx=grep("CA101", names(SO))

M_var_CA101 = generate_coding_matrix_mod( idx, min_depth = 10 )

CA101_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA101[[1]])[1] ), M_var_CA101[[2]] ) )
write_delim(CA101_mut_phyl, "CA101_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA101_mut_phyl <- read.phyDat("CA101_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA101_mut_phyl  <- dist.ml(CA101_mut_phyl)

CA101_treeNJ  <- nj(dm_CA101_mut_phyl)

pdf("Suppl_Figure_1B_CA101_treeNJ.pdf", width = 4, height = 3)
p_CA101 <- ggtree(CA101_treeNJ,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.8, hjust=0.5, 
                                              offset = 0.005, ignore.negative.edge=TRUE ) 
print(p_CA101,vp=vp_phylo)
dev.off()


#===============================================================================
#' Supplementary Figure 1E CA108
#===============================================================================

idx=grep("CA108", names(SO))
M_var_CA108 = generate_coding_matrix_mod( idx, min_depth = 10 )

CA108_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA108[[1]])[1] ), M_var_CA108[[2]] ) )
write_delim(CA108_mut_phyl, "CA108_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA108_mut_phyl <- read.phyDat("CA108_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA108_mut_phyl  <- dist.ml(CA108_mut_phyl)
dm_CA108_mut_phyl[is.na(dm_CA108_mut_phyl)] = 0.2

CA108_treeNJ  <- nj(dm_CA108_mut_phyl)

pdf("Suppl_Figure_1B_CA108_treeNJ.pdf", width = 4, height = 3)
p_CA108 <- ggtree(CA108_treeNJ,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.8, hjust=0.5, 
                                              offset = 0.005, ignore.negative.edge=TRUE ) 
print(p_CA108,vp=vp_phylo)
dev.off()

#===============================================================================
#' Supplementary Figure 1F CA093
#===============================================================================

idx=grep("CA093", names(SO) )
M_var_CA093 = generate_coding_matrix_mod(idx, min_depth=5)

CA093_mut_phyl <- as.data.frame( rbind( c(length(idx), dim(M_var_CA093[[1]])[1] ), M_var_CA093[[2]] ) )
write_delim(CA093_mut_phyl, "CA093_Mut_mod_phyDat.txt", delim = "\t", col_names = F)
CA093_mut_phyl <- read.phyDat("CA093_Mut_mod_phyDat.txt",  type="USER", levels = c(0, 1))
dm_CA093_mut_phyl  <- dist.ml(CA093_mut_phyl)
CA093_treeNJ  <- nj(dm_CA093_mut_phyl)

pdf("Suppl_Figure_1B_CA093_treeNJ.pdf", width = 4, height = 3)
p_CA093 <- ggtree(CA093_treeNJ,options(ignore.negative.edge=TRUE)) + geom_tiplab(angle=90, fontface='plain', size=1.8, hjust=0.5, 
                                              offset = 0.005, ignore.negative.edge=TRUE ) 
print(p_CA093,vp=vp_phylo)
dev.off()

#===============================================================================
#' Plotting of phylogenetic tree in one pdf Supplementary Figure 1
#===============================================================================
pdf("Suppl_Figure_1.pdf", width = 10, height = 10)
pushViewport(viewport(angle = 270, layout = grid.layout(nrow = 3, ncol = 2))) # invert rows and cols due the angle
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p_CA071, vp = vplayout(3, 1))
print(p_CA104, vp = vplayout(2, 1))
print(p_CA113, vp = vplayout(1, 1))
print(p_CA101, vp = vplayout(3, 2))
print(p_CA108, vp = vplayout(2, 2))
print(p_CA093, vp = vplayout(1, 2))
dev.off()
