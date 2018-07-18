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
#   TYPE:           kind of data                                             #
#   Estimate:       correlation coefficient                                  #
#   Naive.S.E:                                                               #
#   Naive.z:        occurence of two groups                                  #
#   Robust.S.E.:    both or each group                                       #
#   Robust.z:       obth or each group                                       #
#   P.value:        p value                                                  #
#   FDR:            adjusted P value by BH									                 #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   pacman; geepack; dplyr; tibble; imputeTS; varhandle                      #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
if(!require(pacman)){
  install.packages("pacman", dependencies = T)
}
pacman::p_load(geepack, dplyr, tibble, imputeTS, varhandle)

args <- commandArgs(T)

if (length(args) < 8) {
  stop("Usage:
       Rscript compare_Result.R 
       phen:   phenotype with sampleID, ID for repeated participants and group(csv file) 
       prof:   profile table rownames->taxonomy; colnames->sampleID
       DNAID:  names of sampleID to connect phen and prof
       GROUP:  names of group information
       TYPE:   the kind of profile 
       grp1:   one of groups to be converted into 0
       grp2:   one of groups to be converted into 1
       out:    result with director and prefix")
}

# prepare for function 
phen <- read.csv(args[1])                               # necessary
prof <- read.table(args[2], header=T, row.names=1)      # necessary
DNAID <- args[3]	 # necessary
GROUP <- args[4]  	 # necessary
TYPE <- args[5]		 # necessary
grp1 <- args[6]		 # necessary	
grp2 <- args[7]		 # necessary	
out <- args[8]		 # necessary

# Judging phenotype with two cols and names are corret
if (!(length(which(colnames(phen) == "SampleID")) > 0)) {
  warning("colnames of Your Phenotype doesn't cmontain SampleID")
  colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
}

if (!(length(which(colnames(phen)=="Stage")) > 0)) {
  warning("colnames of Your Phenotype doesn't contain Stage")
  colnames(phen)[which(colnames(phen)==GROUP)] <- "Stage"
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
  
  phe <- x %>% filter(Stage %in% c(grp1, grp2)) %>%
    arrange(SampleID, Stage) %>%
    select(one_of(id)) %>% 
    filter(!is.na(group)) %>%
    droplevels(.) %>% 
    # convert stage into numeric
    mutate(Stage = case_when( 
      Stage == grp1 ~ 0,
      Stage == grp2 ~ 1))   
  
  prf <- y %>% select(colnames(.)[phe$SampleID]) %>%
    # resevred the rownames
    rownames_to_column(Type) %>%
    # occurence of rows more than 2 
    filter(apply(select(., -one_of(Type)), 1, function(x){sum(x>0)/length(x)}) > 0.2) %>%
    column_to_rownames(Type)  %>% 
    t(.) %>% data.frame(.)
  
  # merge data by sampleID
  dat <- left_join(phe, prf %>% rownames_to_column("SampleID"), by="SampleID") %>%
    select(-one_of("SampleID")) %>% rename("id"="ID") %>% 
    na.replace(., 0)
  
  #scale dat 
  dat$group <- scale(dat$group, scale = T, center = T)
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
    
    res <- cbind(enrich, occ[1], occ[2], med,
              md[1], md[2], men, mn[1], mn[2],
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
idx_names <- colnames(phen)[c(10, 14:33)]
res <- lapply(idx_names, geeFun, phen, prof, TYPE, grp1, grp2)
gee.res <- NULL
for(i in 1:length(res)){
  gee.res <- rbind(gee.res, res[[i]])
}

file <- paste0(out,"/", paste(TYPE, grp1, grp2, "csv", sep="."))
write.csv(gee.res, file, row.names = T)