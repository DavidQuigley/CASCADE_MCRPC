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
#' Uses corrplot.mixed function from corrplot to plot Suppl Figures 4A-B
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. 
#===============================================================================

library( corrplot )
mar=c(0,0,1,0)
# par(mfrow=c(1,2))
par(ask = TRUE)

#===============================================================================
# Plotting patient CA071 in Supp Figure 4A

pdf(height=4, width=4.5, file="Supplementary_Figure_4a.pdf")
corrplot::corrplot.mixed(CA071_correlation_similarity, 
                         lower.col = "black", 
                         upper.col = COL2('PiYG', 10),
                         addCoef.col = 'black', 
                         col.lim=c(0, 1), 
                         tl.pos = c("d"), 
                         lower = "number",
                         number.cex = 0.6,
                         # cl.pos= 'full',
                         mar=c(0,0,2,0),
                         upper = "shade",
                         is.corr = FALSE, 
                         cl.cex = 0.5,
                         tl.cex = 0.65 , 
                         addrect = 2,
                         tl.col = "black") %>%
  corrRect(namesMat = c('CA071-01', 'CA071-08', 'CA071-08', 'CA071-01' ) )
mtext(text = "CA071 AR amplicons Jaccard Breakpoint Index", side = 2, line = 3, las = 3 , cex=1)
mtext(text = 'CA071 AR amplicons Similarity Scores', side = 3, line = 3, las = 1 , cex=1) 
dev.off()

pdf(height=2.3, width=2.6, file="Supplementary_Figure_4b.pdf")

#===============================================================================
# Plotting patient CA108 in Supp Figure 4B

corrplot.mixed(CA108_correlation_similarity, 
               lower.col = "black", 
               upper.col = COL2('PiYG', 10),
               addCoef.col = 'black', 
               col.lim=c(0, 1), 
               tl.pos = c("d"), 
               lower = "number",
               number.cex = 0.6,
               mar=c(0,0,2,0),
               upper = "shade",
               is.corr = FALSE, 
               cl.cex = 0.5,
               tl.cex = 0.65 , 
               addrect = 2,
               tl.col = "black")
mtext(text = "CA108 AR amplicons\nJaccard Breakpoint Index", 
      side = 2, line = 1.5, las = 3, cex=1 )
mtext(text = 'CA108 AR amplicons\nSimilarity Scores', 
      side = 3, line = 2, las = 1 , cex=1)

dev.off()




