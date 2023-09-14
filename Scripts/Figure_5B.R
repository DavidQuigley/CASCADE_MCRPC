#,===============================================================================
#' Description:
#' Reproduce the Figure 5B
#===============================================================================
#' Article: Cascade 
#===============================================================================
#' Author: 
#' - Thaidy Moreno-Rodriguez
#===============================================================================
#' Date:  May 2023
#===============================================================================
#****************************************************************************
#' Uses ggstatsplot to plot Figure 5B using WCDT data
#****************************************************************************
################################################################################

#===============================================================================
#'  Loading libraries. 
#===============================================================================

library( rstatix )
library( ggpubr )
library( ggprism )
library( ggstatsplot ) 
library( gapminder )
library( car )
library( cowplot )



#===============================================================================
# First I check the distribution, normality of the data.
set.seed(12345)

#Test each group for normality
ecDNA_stats_wcdt%>%
    group_by(ecDNA_cat) %>%
    summarise(`W Stat` = shapiro.test(ploidy_PURPLE)$statistic,
              p.value = shapiro.test(ploidy_PURPLE)$p.value)
 
# # A tibble: 2 × 3
# ecDNA_cat `W Stat`  p.value
# <chr>        <dbl>    <dbl>
#   1 ecDNA(+)    0.925 0.000429    
#   2 ecDNA(-)    0.785 0.0000000228


# To check for variance:
leveneTest(data = ecDNA_stats_wcdt, ploidy_PURPLE ~ ecDNA_cat) 
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group   1  0.0014 0.9701
#       133               

# To check for normal distribution:
shapiro.test(ecDNA_stats_wcdt$ploidy_PURPLE[ecDNA_stats_wcdt$ecDNA_cat == "ecDNA(-)"]) 
# Shapiro-Wilk normality test
# 
# data:  ecDNA_stats_wcdt_sub$ploidy_PURPLE[ecDNA_stats_wcdt_sub$ecDNA_cat == "ecDNA_Neg"]
# W = 0.7845, p-value = 2.276e-08
shapiro.test(ecDNA_stats_wcdt$ploidy_PURPLE[ecDNA_stats_wcdt$ecDNA_cat == "ecDNA(+)"]) 
# Shapiro-Wilk normality test
# 
# data:  ecDNA_stats_wcdt$ploidy_PURPLE[ecDNA_stats_wcdt$ecDNA_cat == "ecDNA_Pos"]
# W = 0.92475, p-value = 0.0004288

#===============================================================================
# Perform manually the Mann-Whitney U test
m1<-wilcox.test(ploidy_PURPLE ~ ecDNA_cat, data=ecDNA_stats_wcdt, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m1)
# Wilcoxon rank sum test with continuity correction
# 
# data:  ploidy_PURPLE by ecDNA_cat
# W = 1477, p-value = 0.0004424

#===============================================================================
# Plotting Ploidy 
pw_1 <- ggbetweenstats(
  data             = ecDNA_stats_wcdt,
  x                = ecDNA_cat,
  y                = ploidy_PURPLE,
  var.equal = T, 
  type = "np",
  centrality.plotting = F,
  title            = "wcdt",
  xlab             = "",
  ylab = "Ploidy",
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.8, size = 2, stroke = 0, na.rm = TRUE), 
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  facet.proptest = T
) +
  ggplot2::scale_color_manual(values = c("#9CCB86", "#EB8D71" )) + 
  theme(text = element_text(size = 5),
        plot.subtitle = element_text(size = 3),
        plot.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.1), 
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.1),
        panel.grid.minor.y = element_line(colour = "lightgrey", linetype = "dashed", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) 
ggsave(plot = pw_1,"Figure_5b_Ploidy.pdf", width=2.5,height=2.8)

#===============================================================================
# Statistis for SVtmbPerMb
#===============================================================================

# To check for variance:
leveneTest(data = ecDNA_stats_wcdt, SVtmbPerMb ~ ecDNA_cat) 

# To test each group for normality
ecDNA_stats_wcdt%>%
  group_by(ecDNA_cat) %>%
  summarise(`W Stat` = shapiro.test(SVtmbPerMb)$statistic,
            p.value = shapiro.test(SVtmbPerMb)$p.value)
# # A tibble: 2 × 3
# ecDNA_cat `W Stat`   p.value
# <chr>        <dbl>     <dbl>
# 1 ecDNA(+)     0.889 0.0000149 
# 2 ecDNA(-)     0.856 0.00000210

# To check for distribution
shapiro.test(ecDNA_stats_wcdt$SVtmbPerMb[ecDNA_stats_wcdt$ecDNA_cat == "ecDNA(-)"]) 
shapiro.test(ecDNA_stats_wcdt$SVtmbPerMb[ecDNA_stats_wcdt$ecDNA_cat == "ecDNA(+)"]) 

#===============================================================================
# Perform the Mann-Whitney U test manually

m2<-wilcox.test(SVtmbPerMb ~ ecDNA_cat, data=ecDNA_stats_wcdt, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m2)
# Wilcoxon rank sum test with continuity correction
# 
# data:  SVtmbPerMb by ecDNA_cat
# W = 1577, p-value = 0.002129
#===============================================================================

pw_2 <- ggbetweenstats(
  data             = ecDNA_stats_wcdt,
  x                = ecDNA_cat,
  y                = SVtmbPerMb,
  var.equal = T, 
  type = "np",
  centrality.plotting = F,
  title            = "wcdt",
  xlab             = "",
  ylab = "SV TMBPerMb",
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.8, size = 2, stroke = 0, na.rm = TRUE), 
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  facet.proptest = T
) +
  ggplot2::scale_color_manual(values = c("#9CCB86", "#EB8D71" )) + 
  theme(text = element_text(size = 5),
        plot.subtitle = element_text(size = 3),
        plot.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.1), 
        panel.grid.major.x = element_line(colour = "lightgrey", size = 0.1),
        panel.grid.minor.y = element_line(colour = "lightgrey", linetype = "dashed", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75)) 
ggsave(plot = pw_2, "Figure_5b_SVtmbPerMb.pdf", width=2.5,height=2.8)


#'==============================================================================
#' Chromothripsis 
#'==============================================================================

pw_3 <- ggbarstats(
  data             = ecDNA_stats_wcdt,
  x                = Chromo,
  y                = 'ecDNA_cat',
  title            = "Chromothripsis",
  label.fill.alpha = 0.5,
  xlab             = "",
  ylab             = NULL,
  label.fill.alpha = 0.5,
  label = "counts",
  label.text.size = 20,
  bf.message = FALSE,
  fixed.margin = rows,
  legend.title    = "Alterations", 
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  facet.proptest = TRUE
) + 
  scale_fill_manual(values = alpha(c("#E6A024", "#5BB4E5"), 0.75) ) +
  theme(text = element_text(size = 5),
        plot.subtitle = element_text(size = 3),
        plot.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(colour = "lightgrey", linetype = "dashed", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        legend.position = "none") 

extract_stats(pw_3)
ggsave(plot = pw_3, "Figure_5b_Chromothripsis.pdf", width=1.5,height=2.8)

 
#'==============================================================================
#' TP53 
#'==============================================================================
 
pw_4 <- ggbarstats(
  data             = ecDNA_stats_wcdt,
  x                = TP53,
  y                = 'ecDNA_cat',
  title            = "TP53",
  label.fill.alpha = 0.5,
  xlab             = "",
  ylab             = NULL,
  label.fill.alpha = 0.5,
  label = "counts",
  label.text.size = 20,
  bf.message = FALSE,
  fixed.margin = rows,
  legend.title = NULL, 
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  facet.proptest = TRUE
) + 
  scale_fill_manual(values = alpha(c("#E6A024", "#5BB4E5"), 0.75) ) +
  theme(text = element_text(size = 5),
        plot.subtitle = element_text(size = 3),
        plot.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(colour = "lightgrey", linetype = "dashed", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        legend.position = "none") 

extract_stats(pw_4)
ggsave(plot = pw_4, "Figure_5b_TP53.pdf", width=1.5,height=2.8)

#'==============================================================================
#' PTEN 
#'==============================================================================

pw_5 <- ggbarstats(
  data             = ecDNA_stats_wcdt,
  x                = PTEN,
  y                = 'ecDNA_cat',
  title            = "PTEN",
  label.fill.alpha = 0.5,
  xlab             = "",
  ylab             = NULL,
  label.fill.alpha = 0.5,
  label = "counts",
  label.text.size = 20,
  bf.message = FALSE,
  fixed.margin = rows,
  legend.title = NULL, 
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  facet.proptest = TRUE
) + 
  scale_fill_manual(values = alpha(c("#E6A024", "#5BB4E5"), 0.75) ) +
  theme(text = element_text(size = 5),
        plot.subtitle = element_text(size = 3),
        plot.title = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(colour = "lightgrey", linetype = "dashed", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        legend.position = "none") 

extract_stats(pw_5)
ggsave(plot = pw_5, "Figure_5b_PTEN.pdf", width=1.5,height=2.8)

#'==============================================================================
#' List of plots  
#'==============================================================================

# create a list with the plots
p_list_wcdt<- list(pw_1, pw_2, pw_3, 
                      pw_4, pw_5) 

# extract a legend that is laid out horizontally
legend_w_c_d <- get_legend(
  pw_3 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") )

# combining plots
prow <- cowplot::plot_grid(
  plotlist = p_list_wcdt,   nrow = 1, ncol = 5,
  align = "h",  rel_widths = c(2, 2, 1.3, 1.3,1.3), axis = "b")
# add the legend underneath the row we made earlier resized at 10%

prow2 <- cowplot::plot_grid(prow, legend_c_d, ncol = 1, rel_heights = c(1, .1))
ggsave(plot =prow2, "Figures_5B.pdf", width=9,height=3)
