#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program checks significant difference in Oridination                #
# Args :  5 args                                                             #
#   phen:   phenotype with sampleID, ID for repeated participants and group  #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   TYPE:   the kind of profile                                              #
#   out:    result with director")                                           #
#                                                                            #
# Output: 14 cols result                                                     #
#   TYPE:           kind of data                                             #
#   Estimate:       correlation coefficient                                  #
#   F.Model:        F test                                                   #
#   P.value:        p value                                                  #
#   FDR:            adjusted P value by BH									                 #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#   export PATH=/ldfssz1/ST_META/share/User/zouhua/Miniconda3/bin:$PATH      #
#                                                                            #
# Packages:                                                                  #
#   vegan; dplyr; tibble; imputeTS; readr                                    #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
library(vegan)
library(dplyr)
library(tibble)
library(imputeTS)
library(readr)
library(varhandle)

args <- commandArgs(T)

if (length(args) < 5) {
  stop("Usage:
       Rscript compare_Result.R 
       phen:   phenotype with sampleID, ID for repeated participants and group(csv file) 
       prof:   profile table rownames->taxonomy; colnames->sampleID
       DNAID:  names of sampleID to connect phen and prof
       TYPE:   the kind of profile 
       out:    result with director and prefix")
}

# prepare for function 
phen <- read_csv("../Result/Phenotype/n50.phenotype.merge.csv")                               
prof <- read_delim("../Result/Profile/Species.profile", delim = "\t", col_types = cols()) %>% 
        column_to_rownames("X1")
DNAID <- "SampleID"	 
TYPE <- "ANOSIM"		
out <- "BASE"		 


colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"

TestFun <- function(x, y, type=TYPE) {
  # take a test for oridination result
  #
  # Args:
  #   x:       profile
  #   y:       phenotype
  #   TYPE:    the kind of profile
  #
  # Returns:   
  #     P value of test 
  
  # filter profile
  prf <- y %>% select(colnames(.)[colnames(.) %in% x$SampleID]) %>%
    # resevred the rownames
    rownames_to_column("Type") %>%
    # occurence of rows more than 0.2 
    filter(apply(select(., -one_of("Type")), 1, function(x){sum(x>0)/length(x)}) > 0.2) %>%
    column_to_rownames("Type")  %>% 
    t(.) %>% data.frame(.)
  
  # order phenotype and profile by SampleID
  phe.ord <- x[order(x$SampleID), ]
  prf.ord <- prf[order(rownames(prf)), ]
  
  if (length(!rownames(prf.ord) == phe.ord$SampleID) < nrow(prf.ord)) {
    stop("
         The order of Phenotype and profile are `inconsistant` in SampleID")
  }
  
  PERMANOVAFun <- function(a, b){
    # permanova test for profile and phenotype
    # Args:
    #   a:       profile
    #   b:       phenotype
    #
    # Returns:
    #   PERMANOVA test result
    
    per <- apply(a %>% select(-one_of("SampleID")), 2, function(x, pf){
      dat <- data.frame(y = x, pf) %>% na.omit()
      datphe <- dat[, 1] %>% unfactor()
      if (length(datphe) == 0 | unique(datphe) == 1) {
        res <- data.frame(length(datphe), rep(NA, 6))
        next
      }
      
      if (length(unique(datphe)) < 10) {
        datphe <- as.factor(datphe)
      } 
      
      # distance 
      datprf <- dat[, -1, F] %>% scale(center = F, scale = T)
      dis <- vegdist(datprf, method = "bray")
      
      set.seed(123)
      ad <- adonis(dis ~ datphe, permutations = 999)
      tmp <- as.data.frame(ad$aov.tab) %>% slice(1)
      res <- c(length(datphe), tmp[, 1], tmp[, 2], tmp[, 3],
               tmp[, 4], tmp[, 5], tmp[, 6])
      return(res)
    }, b) %>% t() %>% data.frame()
    
    colnames(per) <- c("SumsOfSample", "Df", "SumsOfSqs", 
                       "MeanSqs", "F.Model", "R2", "Pr(>F)")
    per$FDR <- p.adjust(per$`Pr(>F)`, method = "BH")
    return(per)  
  }
  
  # ANOSIM : Analysis of similarities (ANOSIM) provides a way to test statistically 
  #whether there is a significant difference between two or more groups of sampling units.
  
  ANOSIMFun <- function(a, b){
    # permanova test for profile and phenotype
    # Args:
    #   a:       profile
    #   b:       phenotype
    #
    # Returns:
    #   PERMANOVA test result
    ano <- apply(a %>% select(-one_of("SampleID")), 2, function(x, pf){
      dat <- data.frame(y = x, pf) %>% na.omit()
      datphe <- dat[, 1] 
      if (length(levels(datphe)) > 6 | length(datphe) == 0 | length(unique(datphe)) == 1) {
        res <- c(length(datphe), rep(NA, 3))
      } else {
        datphe <- as.factor(datphe)
        # distance 
        datprf <- dat[, -1, F] %>% scale(center = F, scale = T)
        dis <- vegdist(datprf, method = "bray")
        
        set.seed(123)
        ad <- anosim(dis, datphe, permutations = 999)
        
        res <- c(length(datphe), length(levels(datphe))-1, 
                 as.numeric(ad$statistic), 
                 as.numeric(ad$signif))
      } 
      
      return(res)
    }, b) %>% t() %>% data.frame() %>% na.omit()
    
    colnames(ano) <- c("SumsOfSample", "Df",  "R2", "P.value")
    ano$FDR <- p.adjust(ano$P.value, method = "BH")
    return(ano)  
  }
  
  # which methods of testing to apply
  if (type == "PERMANOVA") {
    res <- PERMANOVAFun(phe.ord, prf.ord)
  } else if (type == "ANOSIM") {
    res <- ANOSIMFun(phe.ord, prf.ord)
  } 
  
  return(res)
  }

res <- TestFun(phen, prof, TYPE)

file <- paste(out, TYPE, "csv", sep = ".")
write.csv(res, file, row.names = T)