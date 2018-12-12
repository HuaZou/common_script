#----------------------------------------------------------------------------#
# Copyright (c) 2018 Hua Zou (BGI-shenzhen). Allrights reserved              #
#   This R program is using to Generate DimensionReduction Result.           #
#                                                                            #   
# Args :  9 args                                                             #
#   phen:   phenotype with sampleID, group more than 2 leves                 #
#   prof:   profile table rownames->taxonomy; colnames->sampleID             #
#   DNAID:  names of sampleID to connect phen and prof                       #
#   GROUP:  names of group information                                       # 
#   TYPE:   the kind of Dimension Reduction methods                          #
#   SCALE:  the variables individually scaled or not                         #   
#   Number: the count of Dimension                                           #
#   Total:  All result or exclude model list(T | F)                          #
#   out:    result with director and prefix")                                #
#                                                                            #
# Output:     RData with all result                                          #
#   RData:     All result in RData                                           #
#                                                                            #
#                                                                            #
# R-version:                                                                 #
#   version.string R version 3.5.1 (2018-07-02)                              #
#                                                                            #
# Packages:                                                                  #
#   dplyr; tibble; readr; imputeTS; factoextra                               #
#   vegan; ape; ca; Rtsne; factoextra; MASS                                  #
#----------------------------------------------------------------------------#

# clear all vectors
rm(list = ls())

# library function
library(dplyr)
library(tibble)
library(imputeTS)
library(readr)
library(factoextra)
library(MASS)
library(vegan)
library(ape)
library(ca)
library(Rtsne)

args <- commandArgs(T)

if (length(args) < 8) {
    stop("Usage:
        Rscript oridination.R 
            phen:   phenotype with sampleID, group(csv file) 
            prof:   profile table rownames->taxonomy; colnames->sampleID
            DNAID:  names of sampleID to connect phen and prof
            GROUP:  names of group information
            TYPE:   the kind of Dimension Reduction methods 
            SCALE:  the variables individually scaled or not(T | F)                        
            Number: the counts of Dimension                                           
            Total:  All result or exclude model list(T | F) 
            TOP:    the toppest eigenvectors in variables                
            out:    result with director")
}

# prepare for function 
phen <- read.csv(args[1])                               
prof <- read_delim(args[2], delim = "\t", col_types = cols()) %>% 
    column_to_rownames("X1")      
DNAID <- args[3]	 
GROUP <- args[4]  	
TYPE <- args[5]
SCALE <- args[6]
Number <- args[7]
Total <- args[8]		  	 
out <- args[9]	

# rename cols
colnames(phen)[which(colnames(phen) == DNAID)] <- "SampleID"
colnames(phen)[which(colnames(phen) == GROUP)] <- "group"


# group levels are less than 2 will stop
if (length(levels(factor(phen$group))) < 2 ) {
    stop("
		`GROUP` no more than 2 levels.\n",
		"please check your `GROUP` levels.\n")
}


# Unconstrained: only based on the species matrix table
#   1. Principal Component Analysis (PCA) 
#   2. Principal Coordinates Analysis (PCOA, MDS)
#   3. Non-Metric Multidimensional Scaling (NMDS)
#   4. Correspondence Analysis (CA)
#   5. Discriminant Analysis (DA)
#   6. t-Distributed Stochastic Neighbor Embedding (Rtsne)
#
# Constrained: Use information from both the species and the environmental matrices.
# the aim of constrained ordination is to find the variability in composition that can
# be explained by the measured environmental variables
#   1. Redundancy Discriminant Analysis (RDA)
#   2. Canonical Correspondence Analysis (CCA)
#   3. Distance-based redundancy analysis (dbRDA)
#   4. Canonical analysis of principal coordinates (CAP)
#
# conclusions:
# db-RDA "can be applied to determine how well additional environmental parameters 
#can explain the variation among objects in the matrix". In contrast, CAP do not try 
#to increase the explanation of the variance with other parameters, but only perform correlations.

OrdinationFun <- function(x, y, type=TYPE, scale=SCALE, number=Number, total=TOTAL){
# Reduce the dimension of profile table with methods  
#
# Args:
#   x:      phenotype which contains group and sampleID for selecting profile
#   y:      profile, a table which cols are sampleID and rows are taxonomy etc
#   type:   the kind of Dimension Reduction methods
#   scale:  the variables individually scaled or not(T | F)
#   number: the counts of Dimension
#   total:  All result or exclude model list(T | F)
#
# Returns:
#   The result of Reduce Dimension 

    phe <- x %>% dplyr::select(SampleID, group) %>% 
        arrange(SampleID, group)
    
    prf <- y %>% dplyr::select(colnames(.)[colnames(.) %in% phe$SampleID]) %>%
        # resevred the rownames
        rownames_to_column("tmp") %>%
        # occurence of rows more than 0.2 
        filter(apply(dplyr::select(., -one_of("tmp")), 1, 
		    function(x){sum(x[!is.na(x)]>0)/length(x[!is.na(x)])}) > 0.2) %>%
        data.frame(.) %>% 
        column_to_rownames("tmp")

    # join phe & prf by sampleID 
    datTol <- left_join(phe, 
  		prf %>% t() %>% data.frame() %>% rownames_to_column("SampleID"),
		by = "SampleID") %>% arrange(SampleID, group) %>% 
        na.replace(., 0) 
    
    idx <- c("SampleID", "group")
    datphe <- datTol %>% dplyr::select(one_of(idx))
    datprf <- datTol %>% dplyr::select(-one_of(idx))

    ScaleFun <- function(dat, scale){
    # scale variables or not   
    #
    # Args:
    #   dat:    table rows -> variables; cols -> samples
    #   scale:  scale or not (T | F)
    #
    # Returns:
    #   A table of scale or not 
        if (scale) {
            scale.res <- apply(dat, 2, function(x){
                scale(x, center=F, scale=T)}) %>% data.frame()
            rownames(scale.res) <- rownames(dat)
            return(scale.res)
        } else {
            return(dat)
        }
    }

    ExtractFun <- function(a, b){
    # filter phen and prof for contrain ordination   
    #
    # Args:
    #   a:  phenotype
    #   b:  profile
    #
    # Returns:
    #   A list of phen and prof
    
        id <- c("SampleID", "group", "BMI", "Tyrosine")
        if ( length(which(colnames(a) %in% id)) < 4) {
    	    stop("
		    `phen` without 4 cols.\n",
		    "please check the number of columns of `phen`.\n")
        }

        # remove missing value in datphen and filter datprof
        datphen <- a %>% mutate(group = as.numeric(group)) %>%
                na.omit() %>% arrange(SampleID) 
        datprof <- b %>% dplyr::select(colnames(.)[colnames(.) %in% datphen$SampleID]) %>%
         # resevred the rownames
        rownames_to_column("tmp") %>%        
        # occurence of rows more than 0.2 
        filter(apply(dplyr::select(., -one_of("tmp")), 1, 
		    function(x){sum(x[!is.na(x)]>0)/length(x[!is.na(x)])}) > 0.2) %>%
        data.frame(.) %>% column_to_rownames("tmp")

        datprof <- data.frame(t(datprof[, order(colnames(datprof))]))
        datphen <- datphen %>% arrange(SampleID)
        ph <- datphen %>% dplyr::select(-one_of("SampleID")) 

        return(list(datprof = datprof, datphen=ph))
    }

    ## unconstrain ordination methods 
    # 1. PCA:  a statistical procedure that uses an orthogonal transformation to 
    # convert a set of observations of possibly correlated variables into a set 
    # of values of linearly uncorrelated variables called principal components.
    PCAFun <- function(m, n, Sca, Num, Tot, top=5){
    # Reduce the dimension of profile table with PCA  
    #
    # Args:
    #   m:    phenotype which contains group and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F) 
    #   Num:  the counts of principal components
    #   Tot:  All result or exclude model list(T | F)
    #   top:  the toppest eigenvectors in variables
    #
    # Returns:
    #   a list of pca result: pca;contribution;explian;score   

        # PCA modeling 
        dat <- ScaleFun(n, Sca)
        pca <- prcomp(dat)

        # each eigenvectors contributions of PCA
        # rotation matrix 
        #   column: contains the principal component loading vector
        #   row:    contains the variables 
        loading <- pca$rotation[, c(1:Num)]
		contribution <- apply(loading, 1, function(x){sum(sapply(x, function(x){x^2}))}) %>% 
            data.frame() %>% rename("value" = ".") %>%  
            rownames_to_column("tmp")%>% 
            arrange(desc(value)) %>% slice(1:top)
        
        # compute standard deviation of each principal component 
        # find the components which explain the maximum variance to retain as 
        # much information as possible using these components
        std_dev <- pca$sdev
        pr_var <- (std_dev^2)[1:Num]
        pr_var_explain <- round(pr_var/sum(std_dev^2), 4) * 100
        explains <- paste0(paste0("PC", seq(Num)), " (", paste0(pr_var_explain, "%"), ")")

        # principal component score of each sample
        score <- data.frame(pca$x[, seq(Num)], group = m[, "group", F])
        colnames(score) <- c(paste0("PCA", seq(Num)), "group")

        # return PCA result 
        if (Tot) {
            pca.res <- list(PCA=pca, contribution=contribution,
                        explains=explains, score=score)
        } else {
            pca.res <- list(contribution=contribution,
                        explains=explains, score=score)
        }
        return(pca.res)
    }
    
    # 2. PCoA: It takes an input matrix giving dissimilarities between pairs of items 
    # and outputs a coordinate matrix whose configuration minimizes a loss function called strain
    PCoAFun <- function(m, n, Sca, Num, Tot, method="euclidean", correction="none"){
    # Reduce the dimension of profile table with PCoA  
    #
    # Args:
    #   m:    phenotype which contains group and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F) 
    #   Num:  the counts of principal coordinate decomposition
    #   Tot:  All result or exclude model list(T | F)
    #   method: methods of calculate the distance between samples
    #
    # Returns:
    #   a list of pcoa result: pcoa/scores/explains
    
        # scale calculate 
        dat <- ScaleFun(n, Sca)

        # table with 0 can't using pcoa 
        # detemine to scale table 
        if (length(dat[dat==0])) {
            dat.no.zero <- ScaleFun(n, TRUE)
        } else {
            dat.no.zero <- dat
        }

        # PCoA modeling
        dis <- vegdist(dat.no.zero, method=method)
        pcoa <- pcoa(dis, correction=correction)
        
        # compute eigenvalues of each principal coordinate decomposition 
        # find the eigenvalues which explain the maximum variance to retain as 
        # much information as possible using these eigenvalues
        eig <- pcoa$values[, "Eigenvalues"]
        eig_var <- eig[1:Num]
        eig_var_explain <- round(eig_var/sum(eig), 4) * 100 
        explains <- paste0(paste0("PCoA", seq(Num)), " (", paste0(eig_var_explain, "%"), ")")

        # principal coordinate decomposition score of each sample
        score <- data.frame(pcoa$vectors[, c(1:Num)], group = m[, "group", F])
        colnames(score) <- c(paste0("PCoA", seq(Num)), "group")

        # return PCoA result 
        if (Tot) {
            pcoa.res <- list(PCoA=pcoa,
                        explains=explains, score=score)
        } else {
            pcoa.res <- list(explains=explains, score=score)
        }
        return(pcoa.res)
    }

    # 3. NMDS: NMDS uses rank orders(take more time to run) not rely on distance, and thus is an 
    # extremely flexible technique that can accommodate a `variety of different kinds` of data
    # NMDS does not use the absolute abundances of species in communities, but rather their rank orders. 
    NMDSFun <- function(m, n, Sca, Num, Tot, method="bray"){
    # Reduce the dimension of profile table with non-metric multidimensional scaling  
    #
    # Args:
    #   m:    phenotype which contains group and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Num:  the counts of principal coordinate decomposition 
    #   Tot:  All result or exclude model list(T | F)
    #   method: methods of calculate the distance between samples
    #
    # Returns:
    #   a list of NMDS results:   

        # scale or not
        # NMDS modeling 
        dat <- ScaleFun(n, Sca)
        mds <- metaMDS(dat, k=as.numeric(Num), trymax=100, distance=method)
    
        # multidimensional's score of each sample
        score <- data.frame(mds$points, group = m[, "group", F])
        colnames(score) <- c(paste0("MDS", seq(Num)), "group")

        # return NMDS result 
        if (Tot) {
            NMDS.res <- list(MDS=mds, score=score)
        } else {
            NMDS.res <- list(score=score)
        }
        return(NMDS.res)
    }

    # 4. CA: Correspondence Analysis 
    CAFun <- function(n, Num, Tot){
    # Reduce the dimension of profile table with CA 
    #
    # Args:
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Num:  the counts of principal coordinate decomposition 
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of CA results:   

        # CA modeling 
        ca <- ca(n, Num)

        # Eigenvalues
        # The proportion of the principal dimensions can inertia explain the variability 
        eigenvalues <- get_eigenvalue(ca)
        eigen_proportion <- round(eigenvalues$variance.percent, 2)[seq(Num)]
        explains <- paste0(paste0("CA", seq(Num)), " (", paste0(eigen_proportion, "%"), ")")

        # the standard coordinates of variable: rows and cols
        score.row <- data.frame(ca$rowcoord)
        colnames(score.row) <- paste0("CA", seq(Num))
        score.col <- data.frame(ca$colcoord)
        colnames(score.col) <- paste0("CA", seq(Num))        

        # return CA result 
        if (Tot) {
            CA.res <- list(CA=ca, score.row=score.row,
                           score.col=score.col)
        } else {
            CA.res <- list(score.row=score.row,
                           score.col=score.col)
        }
        return(CA.res)
    }

    # 5. DA: Discriminant Analysis->canonical discriminant analysis->linear discriminant analysis (LDA)
    #  The purpose of discriminant analysis (LDA) is to find the linear combinations of the original 
    # variables that gives the best possible separation between the groups. 
    DAFun <- function(m, n, Tot){
    # Reduce the dimension of profile table with DA 
    #
    # Args:
    #   m:    phenotype which contains group and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of DA results:   

        # join phe & prf by sampleID 
        dat <- cbind(m, n) %>% 
                dplyr::select(-one_of("SampleID"))
        
        # lda modeling 
        lda <- lda(group~., data=dat)
        lda.pre <- predict(lda, dat)
    
        # multidimensional's score of each sample
        score <- data.frame(lda.pre$x, group = m[, "group", F])

        # return DA result 
        if (Tot) {
            DA.res <- list(DA=lda, score=score)
        } else {
            DA.res <- list(score=score)
        }
        return(DA.res)
    }

    # 6. Rtsne: t-SNE stands for t-Distributed Stochastic Neighbor Embedding  
    RtsneFun <- function(m, n, Sca, Num, Tot, method="bray"){
    # Reduce the dimension of profile table with Rtsne  
    #
    # Args:
    #   m:    phenotype which contains group and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Num:  the counts of Rtsne 
    #   Tot:  All result or exclude model list(T | F)
    #   method: methods of calculate the distance between samples
    #
    # Returns:
    #   a list of Rtsne results:   

        # scale or not
        # Rtsne modeling 
        dat <- ScaleFun(n, Sca)
        dis <- vegdist(dat, method=method)
		Rtsne <- Rtsne(dis, dim=as.numeric(Num), preplexity=30, max_iter=100)
        
		# Rtsne's score of each sample
        score <- data.frame(Rtsne$Y, group = m[, "group", F])
        colnames(score) <- c(paste0("Tsne", seq(Num)), "group")
		
		# return Rtsne result 
        if (Tot) {
            Rtsne.res <- list(Rtsne=Rtsne, score=score)
        } else {
            Rtsne.res <- list(score=score)
        }
        return(Rtsne.res)
    }

    ## constrain ordination methods
    # 1. RDA: Redundancy Discriminant Analysis, based on principal components analysis.
    RDAFun <- function(m, n, Sca, Tot){
    # Attempt to explain difference in species composition between sites by
    # differences in environmental variables(BMI;Amino acids etc)  
    #
    # Args:
    #   m:    phenotype which contains group&environmetal variables and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of RDA result: RDA;explians;score   

        dat <- ExtractFun(m, n)

        # rda modeling
        fm <- formula(dat$datprof ~ BMI + Tyrosine + group)
        rda <- rda(fm, dat$datphen)

        # Eigenvalues
        Eig <- as.numeric(rda$CA$eig)
        Eig_explain <- (round(Eig/sum(Eig), 4) * 100)[1:2]
        explains <- paste0(paste0("CAP", seq(2)), " (", paste0(Eig_explain, "%"), ")")

        # constrain analysis of principal component score of each sample
        score <- data.frame(rda$CA$u[, c(1:2)], group=dat$datphen[, "group", F])
        colnames(score) <- c(paste0("RDA", seq(2)), "group")

        # return RDA result 
        if (Tot) {
            rda.res <- list(RDA=rda,
                        explains=explains, score=score)
        } else {
            rda.res <- list(explains=explains, score=score)
        }
        return(rda.res)
    }

    # 2. CCA: Canonical Correspondence Analysis, based on correspondence analysis.
    CCAFun <- function(m, n, Sca, Tot){
    # Attempt to explain difference in species composition between sites by
    # differences in environmental variables(BMI;Amino acids etc)  
    #
    # Args:
    #   m:    phenotype which contains group&environmetal variables and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of CCA result: CCA;explians;score   

        dat <- ExtractFun(m, n)

        fm <- formula(dat$datprof ~ BMI + Tyrosine + group)
        cca <- cca(fm, dat$datphen)
        
        # Eigenvalues
        Eig <- as.numeric(cca$CCA$eig)
        Eig_explain <- round(Eig/sum(Eig), 4) * 100
        explains <- paste0(paste0("CAP", seq(2)), " (", paste0(Eig_explain[1:2], "%"), ")")

        # constrain analysis of principal component score of each sample
        score <- data.frame(cca$CCA$u[, 1:2], group=dat$datphen[, "group", F])
        colnames(score) <- c(paste0("CCA", seq(2)), "group")

        # return CCA result 
        if (Tot) {
            cca.res <- list(CCA=cca,
                        explains=explains, score=score)
        } else {
            cca.res <- list(explains=explains, score=score)
        }
        return(cca.res)
    }    

    # 3. dbRDA: Distance-based redundancy analysis, based on non-Euclidean distance measures.
    dbRDAFun <- function(m, n, Sca, Tot){
    # Attempt to explain difference in species composition between sites by
    # differences in environmental variables(BMI;Amino acids etc)  
    #
    # Args:
    #   m:    phenotype which contains group&environmetal variables and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of dbRDA result: dbRDA;explians;score   
        dat <- ExtractFun(m, n)

        fm <- formula(dat$datprof ~ BMI + Tyrosine + group)
        dbrda <- capscale(fm, dat$datphen, dist = "bray", dfun = vegdist, sqrt.dist = TRUE)
        
        # Eigenvalues
        Eig <- as.numeric(dbrda$CCA$eig)
        Eig_explain <- round(Eig/sum(Eig), 4) * 100
        explains <- paste0(paste0("dbRDA", seq(3)), " (", paste0(Eig_explain, "%"), ")")

        # constrain analysis of principal component score of each sample
        score <- data.frame(dbrda$CCA$u)

        # return dbRDA result 
        if (Tot) {
            dbrda.res <- list(dbRDA=dbrda,
                        explains=explains, score=score)
        } else {
            dbrda.res <- list(explains=explains, score=score)
        }
        return(dbrda.res)
    }

    # 4. CAP: Constrained Analysis of Principal Coordinates (CAP) is strictly linear and metric.
    # distance: Euclidean(RDA->Redundancy Analysis)/Manhattan/Bray-curtis distance
    # formula: table1->response; table2->vector (including condition:partialled out)
    CAPFun <- function(m, n, Sca, Tot){
    # Attempt to explain difference in species composition between sites by
    # differences in environmental variables(BMI;Amino acids etc)  
    #
    # Args:
    #   m:    phenotype which contains group&environmetal variables and sampleID for selecting profile
    #   n:    profile, a table which cols are sampleID and rows are taxonomy etc
    #   Sca:  the variables individually scaled or not(T | F)
    #   Tot:  All result or exclude model list(T | F)
    #
    # Returns:
    #   a list of CAP result: CPA;explians;score   

        dat <- ExtractFun(m, n)

        fm <- formula(dat$datprof ~ BMI + Tyrosine + Condition(group))
        cap <- capscale(fm, dat$datphen, dist = "bray", dfun = vegdist
                        , add = TRUE, sqrt.dist = TRUE)
        
        # Eigenvalues
        Eig <- as.numeric(cap$CCA$eig)
        Eig_explain <- round(Eig/sum(Eig), 4) * 100
        explains <- paste0(paste0("CAP", seq(2)), " (", paste0(Eig_explain[1:2], "%"), ")")

        # constrain analysis of principal component score of each sample
        score <- data.frame(cap$CCA$u)

        # return CAP result 
        if (Tot) {
            cap.res <- list(CAP=cap,
                        explains=explains, score=score)
        } else {
            cap.res <- list(explains=explains, score=score)
        }
        return(cap.res)
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
