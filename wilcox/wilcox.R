#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
# This R program is using to Generate the result of table of comparison      #
# Args :  11 args                                                            #
#   phen:   phenotype with sampleID, ID for repeated participants and group  #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       #
#   ID:     id for paired test                                               # 
#   TYPE:   the kind of profile                                              #
#   Method: wilcox or ttest                                                  #
#   PAIRED: paired or not default                                            #
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
#   95% CI:     glm result for +/- 1.96									     #
#   compare: compare information                                      		 #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   dplyr; tibble; varhandle; readr                                          #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
suppressPackageStartupMessages(library(dplyr))
library(tibble)
suppressPackageStartupMessages(library(readr))
library(varhandle)
suppressPackageStartupMessages(library(argparser))

# parameter input
parser <- arg_parser("Wilcox.test") %>%
    add_argument("-p", "--phen", 
        help = "phenotype with sampleID, ID for repeated participants and group(csv file)") %>%
    add_argument("-f", "--prof", 
        help = "profile table rownames->taxonomy; colnames->sampleID") %>%
    add_argument("-d", "--DNAID", 
        help = "names of sampleID to connect phen and prof") %>%
    add_argument("-s", "--GROUP", 
        help = "names of group information")  %>%
    add_argument("-i", "--ID", 
        help = "id for paired test")  %>%
    add_argument("-t", "--TYPE", 
        help = "the kind of profile")  %>%
    add_argument("-m", "--METHOD", 
        help = "wilcox or ttest")  %>%
    add_argument("-r", "--PAIRED", 
        help = "paired or not (T or F)")  %>%
    add_argument("-s1", "--group1", 
        help = "one of groups to be converted into 0") %>%
    add_argument("-s2", "--group2", 
        help = "one of groups to be converted into 1")  %>%
    add_argument("-o", "--out", 
        help = "result with director", default = "./")

args <- parse_args(parser)

# prepare for function 
phen <- read.csv(args$p)                              
prof <- read_delim(args$f, col_types = cols(), delim =  "\t") %>% 
			column_to_rownames("X1")      
DNAID <- args$d	 
GROUP <- args$s
ID <- args$i  	
TYPE <- args$t
METHOD <- args$m		 
PAIRED <- args$r	 
grp1 <- args$s1		 
grp2 <- args$s2		 
out <- args$o		

# Judging phenotype with two cols and names are corret
colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"
colnames(phen)[which(colnames(phen) == ID)] <- "ID"

if (length(which(colnames(phen) %in% c("SampleID", "ID", "Stage"))) != 3){
    waring("phenotype without 3 cols: DNAID, ID, GROUP")
}


TestFun <- function(x, y, Type=TYPE, method=METHOD, paired=PAIRED, grp1="BASE", grp2="LOW"){
# calculate the p value and mean abundance of profile 
#
# Args:
#   x:      phenotype which contains stage and sampleID for selecting profile
#   y:      profile, a table which cols are sampleID and rows are taxonomy etc
#   TYPE:   the kind of profile
#   METHOD: wilcox or ttests
#   PAIRED: wilcox test paired or not default
#   grp1:   the stage for arranging to be 0 with default "BASE"
#   grp2:   the stage for arranging to be 1 with default "LOW"
#
# Returns:
#   The p value between two stages
		
    intersectFun <- function(data){
        data <- data %>% mutate(Stage = factor(Stage))
        id <- unique(as.character(data$ID))
        for (i in 1:length(levels(data$Stage))) {
            id <- intersect(id,
            unlist(data %>% filter(Stage == levels(Stage)[i]) %>% select(ID)))
        }
        return(id)
    }

    # get intersect by SampleID
    sid <- intersect(as.character(x$SampleID), colnames(y))  

#    # y filter cols with NAs(more than one)
#    prf <- tbl_df(y) %>% select(colnames(.)[apply(., 2, function(x){!(length(x[is.na(x)]) > 0)})])
    
    # select two stages phenotype and profile
    if (paired) {
        phe <- tbl_df(x) %>% filter(Stage %in% c(grp1, grp2)) %>%
        filter(SampleID %in% sid) %>% 
        select(ID, SampleID, Stage) %>% 
        droplevels() %>% 
        # select intersect ID
		filter(ID %in% intersectFun(.)) %>% 
        arrange(SampleID, Stage) %>%
        # convert stage into numeric
        mutate(group = case_when( 
            Stage == grp1 ~ 0,
            Stage == grp2 ~ 1)) %>%
        data.frame()
    } else {
        phe <- tbl_df(x) %>% filter(Stage %in% c(grp1, grp2)) %>%
        filter(SampleID %in% sid) %>%
        select(SampleID, Stage) %>%
        droplevels() %>%
        arrange(SampleID, Stage) %>%
        # convert stage into numeric
        mutate(group = case_when( 
            Stage == grp1 ~ 0,
            Stage == grp2 ~ 1)) %>%
        data.frame()
    }
    # prf filter by phenotype
    prf <- tbl_df(y) %>% select(colnames(.)[colnames(.) %in% as.character(phe$SampleID)]) %>%
            select(as.character(phe$SampleID)) %>%
            # resevred the rownames
            rownames_to_column(Type) %>% 
            # occurence of rows more than 0.1 
            filter(apply(select(., -one_of(Type)), 1, function(x){sum(x != 0)/length(x)}) > 0.1) %>%
            data.frame(.) %>% 
            column_to_rownames(Type)
	
	# judge no row of profile filter
	if (nrow(prf) == 0) {
		stop("No row of profile to be choosed\n")
	}

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
        setNames(c("Estimate", "Std.Error", "Pr(>|z|)"))

    OR_low_up <- glm_res %>% transmute(
        OR = round(exp(as.numeric(Estimate)), 2), 
        lower = round(exp(as.numeric(Estimate)) - 1.96*as.numeric(Std.Error), 2),
        upper = round(exp(as.numeric(Estimate)) + 1.96*as.numeric(Std.Error), 2)) 
    colnames(OR_low_up) <- c("OR", "lower", "upper")
    OR_lu_res <- OR_low_up %>% mutate("Odds Ratio (95% CI)" = paste0(OR, " (", lower, ";", upper, ")"))
    rownames(OR_lu_res) <- rownames(glm_res)
    OR_lu_res <- OR_lu_res %>% rownames_to_column(Type)

    TestRes <- function(phew, prfw, meth){
    # test by arranging ID
    #
    # Args:
    #   phew:  phenotype
    #   prfw:  profile
    #   meth:  test type 
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

        # order colnames by datphe sampleID
        datprf <- prf %>% select(colnames(.)[colnames(.) %in% as.character(datphe$SampleID)]) %>%
                select(as.character(datphe$SampleID))
        
        # determine the right order and group levels 
        for(i in 1:ncol(datprf)){ 
          if (!(colnames(datprf) == datphe$SampleID)[i]) {
            stop(paste0(i, " Wrong"))
          }
        }
        if (length(levels(datphe$group)) > 2) {
        stop("The levels of `group` are more than 2")
        }
        grp <- datphe[, idx]

        res <- apply(datprf, 1, function(x, grp, paired, method){
            dat <- as.numeric(x)
			# p value;median;mean;enrichment;occurence;FDR

            # t.test or wilcox.test
            if (paired){
                if (method == "wilcox") {
                    p <- wilcox.test(dat ~ grp, paired = T)$p.value	
                } else if (method == "ttest"){
                    p <- t.test(dat ~ grp, paired = T)$p.value	
                }
			} else {
                if (method == "wilcox") {
                    p <- wilcox.test(dat ~ grp, paired = T)$p.value	
                } else if (method == "ttest"){
                    p <- t.test(dat ~ grp, paired = T)$p.value	
                }
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
        }, grp, paired, meth) %>% t(.) %>% data.frame(.) %>%
        rownames_to_column(Type) %>% unfactor(.) 
    res$FDR <- p.adjust(as.numeric(res[, 2]), method="BH")    
    colnames(res)[2:12] <- c("P-value", paste0("Enrichment \n(0,", grp1,";","1,",grp2,")"),
        paste(c(grp1, grp2), "\n occurence"), "Abundance median \n in all", 
        paste("Abundance median \n", c(grp1, grp2)), "Abundance mean \n in all",
        paste("Abundance mean \n", c(grp1, grp2)), "Adjusted P-value")
        
        return(res)
    }
    
    test.res <- TestRes(phe, prf, method)
    
    # cbind wilcox test and Odd ratio 95%CI
    wilcox.CI <- left_join(test.res, OR_lu_res %>% select(c(1, 5)), 
                            by = Type)
    wilcox.CI$Compare <- paste(grp1, "vs", grp2)
    return(wilcox.CI)
}

res <- TestFun(phen, prof, TYPE, METHOD, PAIRED, grp1, grp2) 
file <- paste0(out,"/", paste(TYPE, METHOD, grp1, grp2, "csv", sep="."))
# output
write.csv(res, file, row.names=F)
