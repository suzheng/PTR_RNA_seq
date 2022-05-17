library(ggplot2)
library(dplyr)
library(ggstatsplot)
#library(tidyverse)
library(grid)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(ggthemes)

library("ggsci")
library("scales")
library(ggplot2)
#show_col(c(pal_npg("nrc")(10), "red", "green"))
#col_sz <- pal_npg("nrc")(10)
#con_sz <- c(col_sz[4],rev(brewer.pal(11,'Spectral')[3:10]), col_sz[1])
lwd_box <- 0.9
#Maximum text size for all other text should be 7 pt; minimum text size should be 5 pt
pt_sz <- 7
family_sz <- "sans"
par(family=family_sz)
#theme_set(family=family_sz, size=pt_sz)

col_sz <- pal_jama("default")(7)[c(3,6,2,4,7,5,1)]
con_sz <- c(pal_npg("nrc")(10)[4], rev(brewer.pal(11,'Spectral')[3:10]), pal_npg("nrc")(10)[1])

theme_Publication <- function(base_size=pt_sz, base_family=family_sz) {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "top",
           legend.direction = "horizontal",
           legend.key.size= unit(0.35, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(1,1,1,1),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}