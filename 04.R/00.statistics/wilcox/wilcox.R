#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program is using to Generate the result of table of comparison      #
# Args :  9 args                                                             #
#   phen:   phenotype with sampleID, ID for repeated participants and group  #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       # 
#   TYPE:   the kind of profile                                              #
#   PAIRED: wilcox test paired or not default                                #
#   grp1:   one of groups to be converted into 0                             #
#   grp2:   one of groups to be converted into 1                             #
#   out:    result with director and prefix")                                #
#                                                                            #
# Output: 14 cols result                                                     #
#   TYPE:       kind of data                                                 #
#   P-value:    P by wilcox test                                             #
#   Enrichment: directory of enrichment                                      #
#   Occurence:  occurence of two groups                                      #
#   median:     both or each group                                           #
#   mean:       obth or each group                                           #
#   FDR:        adjusted P value by BH                                       #
#   95% CI:     glm result for +/- 1.96					     #
#   compare: compare information                                             #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   pacman; dplyr; tibble; varhandle; readr                                  #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
if(!require(pacman)){
    install.packages("pacman", dependencies = T)
}
pacman::p_load(dplyr, tibble, varhandle, readr)

args <- commandArgs(T)

if (length(args) < 9) {
    stop("Usage:
        Rscript compare_Result.R 
            phen:   phenotype with sampleID, ID for repeated participants and group(csv file) 
            prof:   profile table rownames->taxonomy; colnames->sampleID
            DNAID:  names of sampleID to connect phen and prof
            GROUP:  names of group information
            TYPE:   the kind of profile 
            PAIRED: wilcox test paired or not (T or F)
            grp1:   one of groups to be converted into 0
            grp2:   one of groups to be converted into 1
            out:    result with director")
}

# prepare for function 
phen <- read_csv(args[1],  col_types = cols())                               
prof <- read_delim(args[2], col_types = cols(), delim =  "\t") %>% 
			column_to_rownames("X1")      
DNAID <- args[3]	 
GROUP <- args[4]  	
TYPE <- args[5]		 
PAIRED <- args[6]	 
grp1 <- args[7]		 
grp2 <- args[8]		 
out <- args[9]		

# Judging phenotype with two cols and names are corret
if (!(length(which(colnames(phen) == "SampleID")) > 0)) {
  warning("colnames of Your Phenotype doesn't contain SampleID")
  colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
}

if (!(length(which(colnames(phen)=="Stage")) > 0)) {
  warning("colnames of Your Phenotype doesn't contain Stage")
  colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"
}

if (length(which(colnames(phen) %in% c("SampleID", "ID", "Stage"))) != 3){
    waring("phenotype without 3 cols: DNAID, ID, GROUP")
}


TestFun <- function(x, y, Type=TYPE, paired=PAIRED, grp1="BASE", grp2="LOW"){
# calculate the p value and mean abundance of profile 
#
# Args:
#   x:  phenotype which contains stage and sampleID for selecting profile
#   y:  profile, a table which cols are sampleID and rows are taxonomy etc
#   TYPE: the kind of profile
#   PAIRED: wilcox test paired or not default
#   grp1: the stage for arranging to be 0 with default "BASE"
#   grp2: the stage for arranging to be 1 with default "LOW"
#
# Returns:
#   The p value between two stages
   
    # select two stages phenotype and profile
    if (paired) {
        phe <- tbl_df(x) %>% filter(Stage %in% c(grp1, grp2)) %>%
        select(ID, SampleID, Stage) %>%
        droplevels() %>%
        arrange(SampleID, Stage) %>%
        # convert stage into numeric
        mutate(group = case_when( 
            Stage == grp1 ~ 0,
            Stage == grp2 ~ 1)) %>%
        data.frame()
    } else {
        phe <- tbl_df(x) %>% filter(Stage %in% c(grp1, grp2)) %>%
        select(SampleID, Stage) %>%
        droplevels() %>%
        arrange(SampleID, Stage) %>%
        # convert stage into numeric
        mutate(group = case_when( 
            Stage == grp1 ~ 0,
            Stage == grp2 ~ 1)) %>%
        data.frame()
    }

    prf <- tbl_df(y) %>% select(colnames(.)[colnames(.) %in% phe$SampleID]) %>%
            select(colnames(.)[order(colnames(.))]) %>%
            # remove missing value in cols (NA)
            #filter(complete.cases(.)) %>% 
            # resevred the rownames
            rownames_to_column(Type) %>%
            # occurence of rows more than 0.1 
            filter(apply(select(., -one_of(Type)), 1, function(x){sum(x>0)/length(x)}) > 0.1) %>%
            data.frame(.) %>% 
            column_to_rownames(Type)
    # scale profile
    prf.cln <- prf %>% t(.) %>% scale(., center = T, scale = T) %>% t(.)
    idx <- which(colnames(phe) == "group")
    
    # glm result for odd ratios 95%CI
    glmFun <- function(m, n){
    # calculate the glm between profile and group information
    #
    # Args:
    #   m:  result of group information which must to be numeric
    #   n:  taxonomy to be glm
    #
    # Returns:
    #   the glm result of between taxonomy group 
        group <- m
        marker <- n
        model <- glm(group ~ marker, family = binomial(link = "logit"))
        cof <- coef(summary(model))["marker", c("Estimate", "Std. Error", "Pr(>|z|)")]
        return(cof) 
    }
    
    glm_res <- t(apply(prf.cln, 1, function(x, group){
        res <- glmFun(group, as.numeric(x))
        return(res)
    }, group = phe[, idx])) %>% data.frame() %>%
        rename("Std.Error" = "Std..Error", "Pr(>|z|)" = "Pr...z..")

    OR_low_up <- glm_res %>% transmute(
        OR = round(exp(as.numeric(Estimate)), 2), 
        lower = round(exp(as.numeric(Estimate)) - 1.96*as.numeric(Std.Error), 2),
        upper = round(exp(as.numeric(Estimate)) + 1.96*as.numeric(Std.Error), 2)) 
    colnames(OR_low_up) = c("OR", "lower", "upper")
    OR_lu_res <- OR_low_up %>% mutate("Odds Ratio (95% CI)" = paste0(OR, " (", lower, ";", upper, ")"))
    rownames(OR_lu_res) <- rownames(glm_res)
    OR_lu_res <- OR_lu_res %>% rownames_to_column(Type)

    WilcoxTest <- function(phew, prfw){
    # paired wilcox test by arranging ID
    #
    # Args:
    #   phew:  phenotype
    #   prfw:  profile 
    #
    # Returns:
    #   wilcox test result
        
        if (paired) {
        datphe <- phew %>% arrange(ID, SampleID) %>%
            mutate_all(funs(as.factor))
        } else {
        datphe <- phew %>% arrange(SampleID, group) %>%
            mutate_all(funs(as.factor))
        }

        datprf <- prfw %>% select(colnames(.)[datphe$SampleID])
        
        if (length(levels(datphe$group)) > 2) {
        stop("The levels of group are more than 2")
        }

        grp <- datphe[, idx]
        res <- apply(datprf, 1, function(x, grp, paired){
            dat <- as.numeric(x)
			# p value;median;mean;enrichment;occurence;FDR
            if (paired){
				p <- wilcox.test(dat ~ grp, paired = T)$p.value	
			} else {
				p <- wilcox.test(dat ~ grp, paired = F)$p.value	
			}
		    med <- median(dat)
            md <- tapply(dat, grp, median)
            if ( md[1] > md[2]) {
                enrich <- levels(grp)[1]
            } else {
                enrich <- levels(grp)[2]
            }
            occ <- tapply(dat, grp, function(x){
                round(sum(x > 0)/length(x), 4)})
            men <- mean(dat)
            mn <- tapply(dat, grp, mean)
            res <- c(p, enrich, occ, med, md, men, mn)
            return(res)
        }, grp, paired) %>% t(.) %>% data.frame(.) %>%
        rownames_to_column(Type) %>% unfactor(.) %>%
        mutate(FDR = p.adjust(as.numeric(X1)))
        
    colnames(res)[2:12] <- c("P-value", paste0("Enrichment \n(0,", grp1,";","1,",grp2,")"),
        paste(c(grp1, grp2), "\n occurence"), "Abundance median \n in all", 
        paste("Abundance median \n", c(grp1, grp2)), "Abundance mean \n in all",
        paste("Abundance mean \n", c(grp1, grp2)), "Adjusted P-value")
        
        return(res)
    }
    test.res <- WilcoxTest(phe, prf)
    
    # cbind wilcox test and Odd ratio 95%CI
    wilcox.CI <- left_join(test.res, OR_lu_res %>% select(c(1, 5)), 
                            by = Type)
    wilcox.CI$Compare <- paste(grp1, "vs", grp2)
    return(wilcox.CI)
}

res <- TestFun(phen, prof, TYPE, PAIRED, grp1, grp2) 
file <- paste0(out,"/", paste(TYPE, grp1, grp2, "csv", sep="."))
# output
write.csv(res, file, row.names=F)