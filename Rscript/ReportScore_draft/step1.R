## script contain 5 steps
### 1. filter ko by all_map_ko.list
### 2. split ko profile&phenotype by  config
### 3. filter ko by number :default 6 
### 4. calculate each ko wilcox.test in two groups 
### 5. calculate pathway and module coverage

##### module1 All_map_ko.list filter ###########################
filterFun <- function(x, y){
  res <- NULL
  for(i in 1:nrow(x)){
    dat <- as.character(x[i, 2])
    tmp <- unlist(strsplit(dat, ","))
    res <- unique(c(res, tmp))
  }
  ko <- intersect(rownames(y), res)
  ko.prf <- y[rownames(y)%in%ko, ]
  return(ko.prf)
}

###### module2 split by configure & filter data by number#########
splitFun <- function(phen, prof, stage=cfg1, n=6, id=cfg2, cg=cfg){
  ## phen : phenotype 
  ## prof : profile
  ## stage : Stage
  ## cfg : factor col
  ## n : filter  number
  ## id : phenotype colname

  group <- c(id, cg)
  datphe <- subset(phen, (Stage==stage), select = group)
  datphe[, 1] <- droplevels(datphe[, 1])
  datpro <- prof[, colnames(prof)%in%datphe[, 1]]
  datphe <- datphe[order(datphe$SampleID), ]
  ko.prf <- datpro[, order(colnames(datpro))]
  ko.prf <- ko.prf[which(rowSums(ko.prf) > 0), ]
  res <- data.frame()
  for(i in 1:nrow(ko.prf)){
    x <- ko.prf[i, ]
    datphe[, 2] <- factor(droplevels(datphe[, 2])) 
    dat1 <- x[which(datphe[, 2]==levels(datphe[, 2])[1])]
    dat2 <- x[which(datphe[, 2]==levels(datphe[, 2])[2])]
    if(length(which(dat1>0)) > n & length(which(dat2>0)) > n){
      res <- rbind(res, x)
    }
  }
  return(list(ko.pf=res, ko.pe=datphe))
}

###### module3 wilcox.test #####################################
testFun <- function(x, y){
  dat.prf <- t(scale(t(x)))
  dat.phe <- y
  logFun <- function(marker, result){
    model <- glm(result ~ marker, family=binomial())
    res <- coef(summary(model))["marker",c("Estimate","Std. Error","Pr(>|z|)")]
    return(res)
  }
  ## glm result
  dat1 <- t(apply(dat.prf, 1, function(x, cfg){
    logFun(as.numeric(x), cfg)}, dat.phe[, 2]))
  OR <- round(exp(as.numeric(dat1[, 1])), 4)
  lower <- round(exp(as.numeric(dat1[,1]) - 1.96*as.numeric(dat1[, 2])), 4)
  uper <- round(exp(as.numeric(dat1[,1]) + 1.96*as.numeric(dat1[, 2])), 4)
  dat2 <- data.frame("OR"=OR, "lower"=lower, "uper"=uper,
                     "com"=paste0(OR, " (", lower, ";", uper, ")"))
  rownames(dat2) <- rownames(dat.prf)
  
  ## wilcox test 
  res <- matrix(0, nrow(dat.prf), 9)
  rownames(res) <- rownames(dat.prf)
  fr <- factor(droplevels(dat.phe[, 2]))
  for(i in 1:nrow(dat.prf)){
    dat <- as.numeric(dat.prf[i, ])
    res[i, 1] <- wilcox.test(dat ~ fr)$p.value
    
    dat.rank <- rank(dat)  ## rank's mean
    tmp1 <- tapply(dat.rank, fr, mean)
    res[i, 2] <- round(tmp1[1])
    res[i, 3] <- round(tmp1[2])
    
    res[i, 4] <- levels(dat.phe[, 2])[1] ## enrichment dir
    if(res[i, 2] > res[i, 3]){ 
      res[i, 5] <- round(qnorm(1- (wilcox.test(dat ~ fr, exact=F)$p.value)/2), 2)   
    }else{
      res[i, 4] <- levels(dat.phe[, 2])[2]
      res[i, 5] <- round(qnorm((wilcox.test(dat ~ fr, exact=F)$p.value)/2), 2)
    }
    
    tmp2 <- tapply(dat, fr, function(x){ ## value > 0 mean
      round(sum(x[x>0])/sum(x), 4)})
    res[i, 6] <- tmp2[1]
    res[i, 7] <- tmp2[2] 
    
    tmp3 <- tapply(dat, fr, mean) ## mean 
    res[i, 8] <- tmp3[1]
    res[i, 9] <- tmp3[2]
  }
  q.value <- p.adjust(res[, 1], method = "BH")
  dat.res <- cbind(res, dat2, q.value)
  colnames(dat.res) <- c("p.value", paste0(levels(dat.phe[, 2]), "_rank"), "Enrichment", "qnorm",
      paste0(levels(dat.phe[, 2]), "_Mean1"), paste0(levels(dat.phe[, 2]), "_Mean2"),
      colnames(dat2), "FDR")
  return(dat.res)
}

##### module4 calculate pathway and module coverage ######## 
### ko occurance in each factor must more than 10%
### calculate each factor occurance rate
coverFun <- function(phen, prof, mk){
  # pro : profile
  # phe : configure
  #data1 : map_KO_info
  #data2 : Module_KO_info
  fr <- factor(droplevels(phen[, 2]))
  ind <- which(phen[, 2]==levels(fr)[1])
  p0 <- prof[, ind]
  p1 <- prof[, -ind]
  
  cutoff <- function(x){
    res <- data.frame()
    for(i in 1:nrow(x)){
      dat <- x[i, ]
      rate <- length(which(dat>0))/ncol(dat)
      if( rate > 0.1 ){
        res <- rbind(res, dat)
      }
    }
    return(res)
  }
  p0.res <- cutoff(p0)
  p1.res <- cutoff(p1)
  
  # compute coverage of pathway and module
  cov <- function(x){
    # x : mpk or mdk
    res <- data.frame()
    for(i in 1:nrow(x)){
      dat <- as.character(x[i, 3])
      ko <- unlist(strsplit(dat, ","))
      n0 <- length(intersect(ko, rownames(p0.res)))
      n1 <- length(intersect(ko, rownames(p1.res)))
      c0 <- n0/x[i, 2]
      c1 <- n1/x[i, 2]
      res <- rbind(res, cbind("cov1"=c0, "cov2"=c1))
    }
    colnames(res) <- c(paste0(levels(fr), "_coverage"))
    return (res)
  }
  cov.path <- cov(mk)
  rownames(cov.path) <- rownames(mk)
  return(cov.pth=cov.path)
}

#########use function ##################
if(type == "path"){
  dat.filter <- filterFun(ampk, pro)
  dat.splt <- splitFun(phe, dat.filter)
  dat.test <- testFun(dat.splt$ko.pf, dat.splt$ko.pe)
  dat.cover <- coverFun(dat.splt$ko.pe, dat.splt$ko.pf, mpki)
  
}else{
  dat.filter <- filterFun(amdk, pro)
  dat.splt <- splitFun(phe, dat.filter)
  dat.test <- testFun(dat.splt$ko.pf, dat.splt$ko.pe)
  dat.cover <- coverFun(dat.splt$ko.pe, dat.splt$ko.pf, mdki)
}

#### output: profile/phenotype/coverage/wilcox.test ###########################################
if(!file.exists(dir)){
  dir.create(dir)
}
name1 <- paste0(dir, "/", prefix,".step1.", type)

#1 ko profile 
file11 <- paste0(name1, ".pro.fiter.csv")
write.csv(dat.splt$ko.pf, file11, quote = F, row.names = T)

#2 phenotype
file12 <- paste0(name1, ".phe.fiter.csv")
write.csv(dat.splt$ko.pe, file12, quote = F, row.names = F)

#3 wilcox.test result
file13 <- paste0(name1, ".wilcox.csv")
write.csv(dat.test, file13, quote = F, row.names = T)

#4 coverage
file14 <- paste0(name1,  ".cover.csv")
write.csv(dat.cover, file14, quote = F, row.names = T)