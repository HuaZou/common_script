#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
#   This R program is using to Generate the result of table of comparison.   #
# first, friedman.test is used to get pvalue; second, if the tax of p value  #
# which are less than 0.05 are used to do next test by PMCMRplus             #
#                                                                            #   
# Args :  7 args                                                             #
#   phen:   phenotype with sampleID, ID and group more than 2 leves          #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       # 
#   TYPE:   the kind of profile                                              #
#   Filter: post test after filtering with p < 0.05 or not                   #
#   out:    result with director and prefix")                                #
#                                                                            #
# Output:     TYPE; p.value for FD&posthc & mean+/-sd                        #              
#   TYPE:      type of names                                                 #
#   Mean+/-sd: mean+/-sd for each levels                                     #
#   P.value: kruskal & wilcox test                                           #
#                                                                            #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   argparser; dplyr; tibble; varhandle; PMCMR; readr                        #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
suppressPackageStartupMessages(library(PMCMR))
suppressPackageStartupMessages(library(dplyr))
library(tibble)
library(imputeTS)
suppressPackageStartupMessages(library(readr))
library(varhandle)
suppressPackageStartupMessages(library(argparser))

# parameter input
parser <- arg_parser("Test for PERMANOVA and etc") %>%
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
    add_argument("-l", "--Filter", 
        help = "post test after filtering with p < 0.05 or not  (TRUE | FALSE)") %>%
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
Filter <- args$l	 	
out <- args$o		

# replace phenotype with two cols and names 
colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
colnames(phen)[which(colnames(phen) == GROUP)] <- "Stage"

# Stage levels and cols of phen less than 3 will stop
if (length(levels(factor(phen$Stage))) < 3 & 
	length(which(colnames(phen) %in% c("SampleID", "ID", "Stage"))) != 3) {
    	stop("
		`GROUP` no more than 3 levels or `phen` without 3 cols.\n",
		"please check your `GROUP` levels or the number of columns of `phen`.\n")
}

TestFun <- function(x, y, Type=TYPE, FILTER=Filter){
  # calculate the p value and mean±sd abundance of profile 
  #
  # Args:
  #   x:  phenotype which contains stage and sampleID for selecting profile
  #   y:  profile, a table which cols are sampleID and rows are taxonomy etc
  #   TYPE: the kind of profile
  #   FILTER: post test after filtering with p < 0.05 or not  (TRUE | FALSE) 
  #
  # Returns:
  #   The p value among more than 2 levels
  
  intersectFun <- function(data){
    data <- data %>% mutate(group = factor(group))
    id <- unique(as.character(data$ID))
    for (i in 1:length(levels(data$group))) {
      id <- intersect(id,
      unlist(data%>% filter(group == levels(group)[i]) %>% select(ID)))
    }
    return(id)
  }

  # intersect of x,y in sampleID
  sid <- intersect(as.character(x$SampleID), colnames(y))

  phe <-  x %>% mutate(group = as.factor(Stage)) %>%
    filter(SampleID %in% sid) %>% 
    select(SampleID, ID, group)
    
  prf <- y %>% select(colnames(.)[colnames(.) %in% as.character(phe$SampleID)]) %>%
    # resevred the rownames
    rownames_to_column("tmp") %>%
    # occurence of rows more than 0.1 
    filter(apply(select(., -one_of("tmp")), 1, 
		function(x){sum(x[!is.na(x)]>0)/length(x[!is.na(x)])}) > 0.1) %>%
    data.frame(.) %>% 
    column_to_rownames("tmp")
  
  # join phe & prf by sampleID 
  datTol <- left_join(phe, 
  		prf %>% t() %>% data.frame() %>% rownames_to_column("SampleID"),
		by = "SampleID") %>% arrange(ID, group) %>%
    # get intersect ID for each group
    filter(ID %in% intersectFun(.))

  idx <- which(colnames(datTol) == "group")
  datprf <- datTol[, 4:(nrow(prf)+3)] %>% t() %>% data.frame()

  fd.res <- apply(datprf, 1, function(x, grp){
    dat <- data.frame(y=as.numeric(x), group=grp) %>% na.omit()
    # p value; mean±sd    
    mn <- tapply(dat$y, dat$group, function(x){
      num = paste(mean(x), "+/-", sd(x))
      return(num)
    })

    dat2 <- data.frame(y=as.numeric(x), group=grp) %>% 
      unstack() %>% t() %>% t() %>% na.omit()
    p <- friedman.test(dat2)$p.value

    res <- c(p, mn)
    return(res)
  }, datTol[, idx]) %>% t(.) %>% data.frame(.) %>%
    rownames_to_column("tmp") %>% unfactor(.) 
  
  fr <- phe$group 
  cl <- unlist(lapply(levels(fr), function(a){paste0(a, "\nMean+/-Sd")}))
  colnames(fd.res)[2:ncol(fd.res)] <- c("P.value.fd", cl)
  
  # filter by p value < 0.05 or run all 
  if (FILTER) {
      fd <- fd.res %>% filter(P.value.fd < 0.05)
      if (nrow(fd) == 0){
          res <- fd.res
          return(res)
      } else {
          fd.prf <- datprf %>% filter(rownames(.) %in% fd$tmp)
      }
  } else {
      fd <- fd.res
      fd.prf <- datprf
  }
  
  post.res <- apply(fd.prf, 1, function(x, grp){
    dat <- data.frame(y=as.numeric(x), group=grp) %>% 
      unstack() %>% t() %>% t() %>% na.omit()
    # p value; mean+/-sd
    p <- posthoc.durbin.test(dat, p.adj="none")$p.value
    p.val <- p[!upper.tri(p)]
    
    # mn <- tapply(as.numeric(x), grp, function(x){
    #   num = paste(mean(x), "+/-", sd(x))
    #   return(num)
    # })
    # res <- c(p.val, mn)
    return(c(p.val, nrow(dat)))
  }, datTol[, idx]) %>% t(.) %>% data.frame(.) %>% 
    unfactor(.)
  
  # names 
  rownames(post.res) <- fd$tmp
  n <- length(levels(fr))
  cl <- NULL
  for(i in 1:(n-1)){
    for(j in (i+1):n){
       cl <- c(cl, paste("P.value\n", levels(fr)[i], "vs", levels(fr)[j]))
    }
  }
  colnames(post.res)[1:(n+1)] <- c(cl, "Times")

  # cbind results
  res <- left_join(fd, post.res %>% rownames_to_column("tmp"), by = "tmp")
  colnames(res)[which(colnames(res)=="tmp")] <- Type

  return(res)
}

res <- TestFun(phen, prof, TYPE, Filter) 
file <- paste0(out,"/", paste(TYPE, Filter, "csv", sep="."))
# output
write.csv(res, file, row.names=F)