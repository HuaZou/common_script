#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program Calculate Association by PGEE for longitude data            #
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
#   pacman; PGEE; dplyr; tibble; imputeTS                                    #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
if(!require(pacman)){
  install.packages("pacman", dependencies = T)
}
pacman::p_load(PGEE, dplyr, tibble, imputeTS)

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

PgeeFun <- function(tag, x, y, Type=TYPE, grp1="BASE", grp2="LOW"){
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
  id <- c("SampleID", "Stage", "ID", "group")
  
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
    select(-one_of("SampleID")) %>% rename("id"="ID", "y"="group") %>% 
    na.replace(., 0)
  
  #scale dat 
  dat$id <- as.numeric(sub("P","", dat$id))
  dat$y <- scale(dat$y, scale = T, center = T)
  dat[, 4:ncol(dat)] <- scale(dat[, 4:ncol(dat)], scale = T, center = T)
  
  formula <- "y ~. -id"
  family <- gaussian(link = "identity")
  lambda.vec <- seq(0.01, 0.25, 0.01)
  
  # best cv 
  cv <- CVfit(formula = formula, id = id, data = dat, family = family, scale.fix = TRUE,
              scale.value = 1, fold = 5, lambda.vec = lambda.vec, pindex = c(1, 2),
              eps = 10^-6, maxiter = 30, tol = 10^-6)
  
  # fit PGEE model 
  fit <- PGEE(formula = formula, id = id, data = dat, na.action = NULL,
              family = family, corstr = "independence", Mv = NULL,
              beta_int = c(rep(0,dim(dat)[2]-1)), R = NULL, scale.fix = TRUE,
              scale.value = 1, lambda = cv$lam.opt, pindex = c(1,2), eps = 10^-6,
              maxiter = 30, tol = 10^-6, silent = TRUE)
  
  # result 
  ## see a portion of the results returned by coef(summary(fit))
  res <- data.frame(coef(summary(fit))) %>% 
    rownames_to_column(Type) %>%
    mutate(p.value = 2*pnorm(-abs(Robust.z))) %>%
    filter(p.value < 0.05) %>%
    mutate(FDR = p.adjust(as.numeric(p.value))) %>%
    slice(-c(1:2)) %>% 
    mutate(Phenotype = tag)
  
  return(res)
}

# phenotype index to run PGEE
idx_names <- colnames(phen)[c(10, 14:33)]
res <- lapply(idx_names, PgeeFun, phen, prof, TYPE, grp1, grp2)
pgee.res <- NULL
for(i in 1:length(res)){
  pgee.res <- rbind(pgee.res, res[[i]])
}

file <- paste0(out,"/", paste(TYPE, grp1, grp2, "csv", sep="."))
write.csv(pgee.res, file, row.names = T)