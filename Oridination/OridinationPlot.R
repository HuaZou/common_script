#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
#   This R program is using to Generate DimensionReduction Result.           #
#                                                                            #   
# Args :  5 args                                                             #
#   RData:  Data with all result                                             #
#   TYPE:   type of RData                                                    #
#   TOTAL:  All result or exclude model list(T | F)                          #
#   Title:  the title of plot                                                # 
#   out:    result with director and prefix")                                #
#                                                                            #
# Output:     plots                                                          #
#   plot1:   show by score                                                   #
#   plot2:   show by model list                                              #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   dplyr; tibble; ggplot2; imputeTS; factoextra                             #                                  #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
library(dplyr)
library(tibble)
library(ggplot2)
library(imputeTS)
library(factoextra)

args <- commandArgs(T)

if (length(args) < 5) {
    stop("Usage:
        Rscript OridinationPlot.R
          RData:  Data with all result                                             
          TYPE:   type of RData                                                    
          TOTAL:  All result or exclude model list(T | F)                          
          Title:  the title of plot                                                
          out:    result with director and prefix")                                
}

# prepare for function 
load("demo/PCA.F.RData")
TYPE <- "PCA"
TOTAL <- T
Title <- "PCA"
out <- "./"
PlotFun <- function(x, type=TYPE, total=TOTAL, title=Title){
# Display the plots of different methods of oridination  
#
# Args:
#   RData:  Data with all result                                             
#   TYPE:   type of RData                                                    
#   TOTAL:  All result or exclude model list(T | F)                          
#   Title:  the title of plot  
#
# Returns:
#   The list of plots 
  
    PlotScoreFun <- function(x, title=Title){
    # the plots of  score
    #
    # Args:
    #   x:  RData with all result                                 
    #   Title:  the title of plot 
    #
    # Returns:
    #   Score plot   
      
      dat <- dat.res$score
      colnames(dat)[1:2] <- c("x", "y")
      xlab <- dat.res$explains[1]
      ylab <- dat.res$explains[2]
      
      pl <- ggplot(dat, aes(x=x, y=y))+
        geom_point(aes(color=group), alpha=1, size=2, shape=16)+
        guides(fill=F, color=guide_legend(title = NULL,keywidth=.6,keyheight=.6))+
        theme_bw()+
        labs(title=Title, x=xlab, y=ylab)+
        theme(plot.title = element_text(size = 12, hjust = .5, face = "bold"), 
              axis.title = element_text(size = 10, face = "bold"),
              panel.grid = element_blank(),  
              axis.text = element_text(size = 8),
              legend.text = element_text(size = 8),
              legend.position = c(1,1),
              legend.justification = c(1,1),
              legend.background = element_rect(fill="white", color = "black"))
      return(pl)
    }


    # which methods of multidimensional scaling to apply
    if (type == "PCA") {
        res <- PCAFun(datphe, datprf, scale, number, total)
    } else if (type == "PCoA") {
        res <- PCoAFun(datphe, datprf, scale, number, total)
    } else if (type == "NMDS") {
        res <- NMDSFun(datphe, datprf, scale, number, total)
    } else if (type == "CA") {
        res <- CAFun(datprf, number, total)
    } else if (type == "DA") {
        res <- DAFun(datphe, datprf, total)
    } else if (type == "Rtsne") {
        res <- RtsneFun(datphe, datprf, scale, number, total)
    } else if (type == "RDA") {
        res <- RDAFun(x, y, scale, total)
    } else if (type == "CCA") {
        res <- CCAFun(x, y, scale, total)
    } else if (type == "dbRDA") {
        res <- dbRDAFun(x, y, scale, total)
    } else if (type == "CAP") {
        res <- CAPFun(x, y, scale, total)
    }

    return(res)
}

dat.res <- OrdinationFun(phen, prof, TYPE, SCALE, Number, Total)
file <- paste0(out,"/", paste(TYPE, Total, "RData", sep="."))
save(dat.res, file=file)
