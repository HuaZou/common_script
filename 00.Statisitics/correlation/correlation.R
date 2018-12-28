#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program Calculate Association by GEE for longitude data             #
# Args :  7 args                                                             #
#   phen:   phenotype with sampleID, ID for repeated participants and group  #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       #
#   TYPE:   the kind of profile                                              #
#   level:    one of groups                                                  #
#   out:    result with director")                                           #
#                                                                            #
# Output: 14 cols result                                                     #
#   TYPE:        kind of data                                                #
#   Estimate:    correlation coefficient                                     #
#   Std.err:                                                                 #
#   Wald:        wald test value                                             #
#   Pr(>|W|):    wald test p value                                           #
#   Phenotype:   phenotype name                                              #
#   FDR:            adjusted P value by BH								     #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   dplyr; tibble; imputeTS; varhandle; readr; argparser                     #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
suppressPackageStartupMessages(library(dplyr))
library(tibble)
library(imputeTS)
suppressPackageStartupMessages(library(readr))
library(varhandle)
suppressPackageStartupMessages(library(argparser))

# parameter input
parser <- arg_parser("Spearman correlation coefficient") %>%
    add_argument("-p", "--phen", 
        help = "phenotype with sampleID, ID for repeated participants and group(csv file)") %>%
    add_argument("-f", "--prof", 
        help = "profile table rownames->taxonomy; colnames->sampleID") %>%
    add_argument("-d", "--DNAID", 
        help = "names of sampleID to connect phen and prof") %>%
    add_argument("-s", "--GROUP", 
        help = "names of group information")  %>%
    add_argument("-t", "--TYPE", 
        help = "the kind of profile")  %>%
    add_argument("-l", "--level", 
        help = "one of groups") %>%
    add_argument("-o", "--out", 
        help = "result with director", default = "./")

args <- parse_args(parser)

# prepare for function 
phen <- read.csv(args$p)                             
prof <- read_delim(args$f, delim = "\t", col_types = cols()) %>% 
        column_to_rownames("X1")
DNAID <- args$d	 
GROUP <- args$s  	 
TYPE <- args$t		
level <- args$l	 		 	
out <- args$o		 

# Judging phenotype with two cols and names are corret
colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"

# phen less than 3 will stop
if (length(which(colnames(phen) %in% c("SampleID", "ID", "Stage"))) != 3) {
    	stop("`phen` without 3 cols.\n",
		"please check your `GROUP` levels or the number of columns of `phen`.\n")
}

corFun <- function(tag, x, y, Type=TYPE, level="BASE"){
  # Calculate the corretion coefficient between two stage by PGEE
  #
  # Args:
  #   tag:     one of responses in phenotype	
  #   x:       profile
  #   y:       phenotype
  #   TYPE:    the kind of profile
  #   level:   the stage for arranging to be 0 with default "BASE"
  #
  # Returns:   
  #     a table of estimate and p.value less than 0.05
  
  # choose data 
  colnames(x)[which(colnames(x)==tag)] <- "group" 
  id <- c("SampleID", "Stage", "ID", "Age", "Gender", "group")
  # get intersect by SampleID
  sid <- intersect(as.character(x$SampleID), colnames(y))

  phe <- x %>% filter(Stage %in% level) %>%
    filter(SampleID %in% sid) %>%
    arrange(SampleID, Stage) %>%
    select(one_of(id)) %>% 
    filter(!is.na(group)) %>%
    droplevels(.)   
  
  prf <- y %>% select(colnames(.)[colnames(.) %in% phe$SampleID]) %>%
    # resevred the rownames
    rownames_to_column(Type) %>%
    # occurence of rows more than 2 
    filter(apply(select(., -one_of(Type)), 1, function(x){sum(x>0)/length(x)}) > 0.2) %>%
    column_to_rownames(Type)  %>% 
    t(.) %>% data.frame(.)
  
  # merge data by sampleID
  dat <- left_join(phe, prf %>% rownames_to_column("SampleID"), by="SampleID") %>%
    select(-one_of("SampleID")) %>% rename("id"="ID") %>% 
    na.replace(., 0) %>% 
    mutate(group = scale(log(group), scale = T, center = T))
  
  # choose dat 
  datprf <- dat[, c(6:ncol(dat))] 
  datphe <- dat[, c(1:5)]  
  
  # cor calculate
  res <- apply(t(datprf), 1, function(x, phs){
    d <- cbind(datphe, y = as.numeric(x))
    # median;mean;occurence
    med <- median(d$y)
    occ <- round(sum(d$y > 0)/length(d$y), 4)
    men <- mean(d$y)
    d$y[which(d$y == 0)] <- 1e-10
    #print(d)
    d$y <- scale(log(d$y), center = T, scale = T)
    cor <- cor.test(as.numeric(d$y), 
                    as.numeric(d$group), method = "spearman")
    
    res <- c(as.numeric(cor$estimate), 
             as.numeric(cor$p.value), 
             occ, med,  men)
    return(res)
  }, datphe) %>% t(.) %>% data.frame(.) %>%
    rownames_to_column(Type) %>% 
    mutate(Phentype = tag)
    
  colnames(res)[2:6] <- c("Estimate", "P.value", 
                          paste(level, "\n Occurence"),
                          paste(level, "\n Abundance median"),
                          paste(level, "\n Abundance mean"))

  res$FDR <- p.adjust(res$P.value, method = "BH")
  return(res)
}

# phenotype index to run spearman correlation coefficient
idx_names <- colnames(phen)[c(10, 14:33)]
res <- lapply(idx_names, corFun, phen, prof, TYPE, level)
list.res <- NULL
for(i in 1:length(res)){
  list.res <- rbind(list.res, res[[i]])
}

file <- paste0(out,"/", paste(TYPE, level, "csv", sep="."))
write.csv(list.res, file, row.names = F)
