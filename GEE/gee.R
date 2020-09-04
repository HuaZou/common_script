#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program Calculate Association by GEE for longitude data             #
# Args :  8 args                                                             #
#   phen:   phenotype with sampleID, ID for repeated participants and group  #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       #
#   TYPE:   the kind of profile                                              #
#   grp1:   one of groups to be converted into 0                             #
#   grp2:   one of groups to be converted into 1                             #
#   out:    result with director")                                           #
#                                                                            #
# Output: 14 cols result                                                     #
#   TYPE:        kind of data                                                #
#   Estimate:    correlation coefficient                                     #
#   Std.err:                                                                 #
#   Wald:        wald test value                                             #
#   Pr(>|W|):    wald test p value                                           #
#   Phenotype:   phenotype name                                              #
#   FDR:            adjusted P value by BH									                 #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   geepack; dplyr; tibble; imputeTS; varhandle; readr; argparser            #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
suppressPackageStartupMessages(library(geepack))
suppressPackageStartupMessages(library(dplyr))
library(tibble)
library(imputeTS)
suppressPackageStartupMessages(library(readr))
library(varhandle)
suppressPackageStartupMessages(library(argparser))

# parameter input
parser <- arg_parser("Generalized Estimating Equation") %>%
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
    add_argument("-s1", "--group1", 
        help = "one of groups to be converted into 0") %>%
    add_argument("-s2", "--group2", 
        help = "one of groups to be converted into 1")  %>%
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
grp1 <- args$s1	 	
grp2 <- args$s2	 	
out <- args$o		 

# change colnames
colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"

# phen less than 3 will stop (ID for repeat measure)
if (length(which(colnames(phen) %in% c("SampleID", "ID", "Stage"))) != 3) {
    	stop("`phen` without 3 cols.\n",
		"please check your `GROUP` levels or the number of columns of `phen`.\n")
}

geeFun <- function(tag, x, y, Type=TYPE, grp1="BASE", grp2="LOW"){
  # Calculate the corretion coefficient between two stage by PGEE
  #
  # Args:
  #   tag:     one of responses in phenotype	
  #   x:       profile
  #   y:       phenotype
  #   TYPE:    the kind of profile
  #   grp1:    the stage for arranging to be 0 with default "BASE"
  #   grp2:    the stage for arranging to be 1 with default "LOW" 
  #
  # Returns:   
  #     a table of estimate and p.value less than 0.05
  
  # choose data 
  colnames(x)[which(colnames(x)==tag)] <- "group" 
  id <- c("SampleID", "Stage", "ID", "Age", "Gender", "group")
   
  # intersect
  intersectFun <- function(data){
    data <- data %>% mutate(Stage = factor(Stage))
    id <- unique(as.character(data$ID))
    for (i in 1:length(levels(data$Stage))) {
      id <- intersect(id,
        as.character(unlist(data %>% filter(Stage == levels(Stage)[i]) %>% select(ID))))
    }
    return(id)
  }

  phe <- x %>% filter(Stage %in% c(grp1, grp2)) %>%
    filter(ID %in% intersectFun(.)) %>%
    arrange(SampleID, Stage) %>%
    select(one_of(id)) %>% 
    filter(!is.na(group)) %>%
    droplevels(.) %>% 
    # convert stage into numeric
    mutate(Stage = case_when( 
      Stage == grp1 ~ 0,
      Stage == grp2 ~ 1))   
  
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
    mutate(group = scale(group, scale = T, center = T))
  
  #scale dat 
  datprf <- dat[, c(6:ncol(dat))]
  datphe <- dat[, c(1:5)]  
  datphe$Stage <- factor(datphe$Stage)
  
  # gee calculate
  res <- apply(t(datprf), 1, function(x, datphe){
    d <- cbind(datphe, y = as.numeric(x))
    # median;mean;enrichment;occurence
    med <- median(d$y)
    md <- tapply(d$y, datphe$Stage, median)
    if ( md[1] > md[2]) {
      enrich <- levels(datphe$Stage)[1]
    } else {
      enrich <- levels(datphe$Stage)[2]
    }
    occ <- tapply(d$y, datphe$Stage, function(x){
      round(sum(x > 0)/length(x), 4)})
    men <- mean(d$y)
    mn <- tapply(d$y, datphe$Stage, mean)
    # GEE
    d$y <- scale(d$y, scale = T, center = T)
    fm <- formula(group ~ y+Age+Gender)
    fit <- geeglm(fm, family=gaussian, id=id, data=d,
                  corstr="exchangeable")
    res.fit <- data.frame(coef(summary(fit))) %>% 
      rename("P.value" = "Pr...W..") %>% 
      slice(2) %>% data.frame(.)
    
    res <- c(enrich, occ, med, md, men, mn,
              res.fit[, 1], res.fit[, 2], res.fit[, 3], 
              res.fit[, 4], tag)
    return(res)
  }, datphe) %>% t(.) %>% data.frame(.) %>%
    rownames_to_column(Type)
  
  colnames(res)[2:15] <- c(paste0("Enrichment \n(0,", grp1,";","1,",grp2,")"),
        paste(c(grp1, grp2), "\n occurence"), "Abundance median \n in all", 
        paste("Abundance median \n", c(grp1, grp2)), "Abundance mean \n in all",
        paste("Abundance mean \n", c(grp1, grp2)), 
        "Estimate", "Std.err", "Wald", "Pr(>|W|)", "Phenotype")
  
  res$FDR <- p.adjust(res$`Pr(>|W|)`, method = "BH")
  return(res)
}

# phenotype index to run GEE
idx_names <- colnames(phen)[c(10, 14:44)]
res <- lapply(idx_names, geeFun, phen, prof, TYPE, grp1, grp2)
gee.res <- NULL
for(i in 1:length(res)){
  gee.res <- rbind(gee.res, res[[i]])
}

file <- paste0(out,"/", paste(TYPE, grp1, grp2, "csv", sep="."))
write.csv(gee.res, file, row.names = F)