## script contain 5 steps
### 1. test result wilcox/ttest
### 2. count how many ko to one map or module
### 3. permutation : based on map count / based on map 

########### module1 test result #############################
testFun <- function(x, y, test=method, pair=paired){
  phen <- x[order(rownames(x)), , F]
  prof <- y[, order(colnames(y))]
  fr <- factor(phen[, 1])
  
  # result:
  ### 1:two side  2:less side  3:greater side 4/5/6:FDR 6/7/8: z_score
  ### 10/11: mean   12/13:sd 14/15: no-zore num 16:significant difference(1/2, 0)
  res <- matrix(0, nrow(prof), 16)
  rownames(res) <- rownames(prof)
  
  for(i in 1:nrow(prof)){
    dat <- as.numeric(prof[i, ])
    if(test == "wilcox"){
      res[i, 1] <- wilcox.test(dat ~ fr, paired=as.logical(pair), alternative="two.sided")$p.value
      res[i, 2] <- wilcox.test(dat ~ fr, paired=as.logical(pair), alternative="less")$p.value
      res[i, 3] <- wilcox.test(dat ~ fr, paired=as.logical(pair), alternative="greater")$p.value
    }else{
      res[i, 1] <- t.test(dat ~ fr, paired=as.logical(pair), alternative="two.sided")$p.value
      res[i, 2] <- t.test(dat ~ fr, paired=as.logical(pair), alternative="less")$p.value
      res[i, 3] <- t.test(dat ~ fr, paired=as.logical(pair), alternative="greater")$p.value
    }
    
    tmp1 <- tapply(dat, fr, mean)
    res[i, 10] <- tmp1[1]
    res[i, 11] <- tmp1[2]
    
    tmp2 <- tapply(dat, fr, sd)
    res[i, 12] <- tmp2[1]
    res[i, 13] <- tmp2[2]
    
    tmp3 <- tapply(dat, fr, function(x){sum(x>0)})
    res[i, 14] <- tmp3[1]
    res[i, 15] <- tmp3[2]
    
    if(res[i, 1] < 0.05){
      if(res[i, 2] < 0.05){
        res[i, 16] <- levels(fr)[1] 
      }else{
        res[i, 16] <- levels(fr)[2] 
      }
    }else{
      res[i, 16] <- "NS"
    }
  }
  
  ### p.adjust-->q.value and get z-score
  zmax <- 8.20953
  zmin <- -8.20953
  for(i in 1:3){
    x <- as.numeric(as.character(res[, i]))
    x <- x[!is.na(x)]
    res <- res[!is.na(res[, 1]), ]
    x.a <- p.adjust(x, method = "BH")
    res[, 3+i] <- x.a
    z.a <- qnorm(1-x)
    for(j in 1:length(z.a))
    {
      if((z.a[j] > 0) && is.infinite(z.a[j])){
        z.a[j] <- zmax
      }
      if((z.a[j]< 0) && is.infinite(z.a[j])){
        z.a[j] <- zmin
      }
    }
    res[, 6+i] <- z.a
  }
  colnames(res) <- c("p_value_twosides", "p_value_less", "p_value_greater",
                     "p_adjust_twosides", "p_adjust_less", "p_adjust_greater", 
                     "z_score(twosides)", "z_score(less)", "z_score(greater)", 
                     paste0("mean_", levels(fr)), paste0("sd_", levels(fr)), 
                     paste0("Occ_", levels(fr)),"enrichment_dir")
  return(data.frame(res))
}

##### permutation: based number/based map#################
permFun <- function(x, y, phen=phe, prof=pro, times=1000){
  
  phen <- phen[order(rownames(phen)), ,F]
  prof <- prof[, order(colnames(prof))] 
  ## each map has how many kos
  map_ko <- data.frame()
  for(i in 1:nrow(y)){
    dat <- as.character(y[i, 2])
    tmp1 <- unlist(strsplit(dat, ","))
    array <- intersect(tmp1, rownames(x))
    tmp2 <- paste0(array, collapse = ";")
    num <- length(array)
    if(num > 0){
      map_ko <- rbind(map_ko, cbind("ID"=rownames(y)[i], "Times"=num, "KO"=tmp2))
    }
  }
  ## permutition by ko numbers
  col <- unique(as.character(map_ko$Times))
  datn <- x[, c(7:9)]
  dat.cln <- apply(datn, 2, as.numeric)
  rownames(dat.cln) <- rownames(datn)
  perm_num <- matrix(NA, length(col), ncol(dat.cln)*2+1)
  for(i in 1:length(col)){
    value <- as.numeric(col[i])
    tmp <- rep(0, times)
    perm_num[i, 1] <- value
    for(j in 1:ncol(datn)){
      for(k in 1:times){
        perm <- dat.cln[sample(nrow(dat.cln), value), j]
        tmp[k] <- sum(perm)/sqrt(value)
      }
      perm_num[i, j*2] <- mean(tmp)
      perm_num[i, j*2+1] <- sd(tmp)
    }
  }
  colnames(perm_num) <- c("Map_KO_Number", paste0("two",c("_Mean", "_Sd")), 
                      paste0("less",c("_Mean", "_Sd")), paste0("greater",c("_Mean", "_Sd")))
  perm_num <- perm_num[order(perm_num[, 1]), ] 
  
  ##### permutation by map ######################################
  # p0 p1 : 通过检验P值判断KO的富集方向和数目  p.ajust  permutation
  fr <- factor(droplevels(phen[, 1]))
  perm_map <- matrix(NA, nrow(y), 11)
  
  for(m in 1:nrow(y)){
    dat1.cln <- as.character(y[m, 2])
    tmp3 <- unlist(strsplit(dat1.cln, ","))
    ko.arr <- intersect(tmp3, rownames(pro))
    value2 <- length(ko.arr) 
    perm_map[m, 1] <- rownames(y)[m]
    perm_map[m, 2] <- value2
    if(value2==0){
      perm_map[m, 3:11] <- NA
    }else{
      occ <- x[rownames(x)%in%ko.arr, , F]
      occ.cln <- sapply(occ[, -16, F], function(x){as.numeric(as.character(x))})
    if(value2==1){
        perm_map[m, 3] <- length(occ.cln[[14]]>0)
        perm_map[m, 4] <- length(occ.cln[[15]]>0)
      }else{
        rownames(occ.cln) <- rownames(occ)
        num <- apply(occ.cln[, c(14:15)], 2, function(x){length(x>0)})
        perm_map[m, 3] <- num[[1]]
        perm_map[m, 4] <- num[[2]]
      }
      NS <- length(which(occ[, 16] == "NS"))
      p0 <- length(which(occ[, 16] == levels(fr)[1]))
      p1 <- length(which(occ[, 16] == levels(fr)[2]))
      rich <- paste(NS, p0, p1, sep = "; ")
      perm_map[m, 5] <- rich
      
      for(n in 1:3){
        if(value2 == 1){tmp4 <- sum(occ.cln[[n+6]])
        }else{
          tmp4 <- sum(occ.cln[, n+6])/sqrt(value2)}
          perm_map[m, n+5] <- tmp4
          perm_num_val <- perm_num[which(perm_num[, 1]==value2), ]
          tmp5 <- (tmp4 - perm_num_val[[n*2]])/(perm_num_val[[n*2+1]]) 
          perm_map[m, n+8] <- tmp5
         }
     }
    }
    perm_map_res <- data.frame(perm_map[which(perm_map[, 2]!=0), ], stringsAsFactors = FALSE)
    colnames(perm_map_res) <- c("ID", "KO_sum", paste0(levels(fr), "_KO_Num"),
          paste("NS", levels(fr)[1], levels(fr)[2], sep = ";"), 
          paste0("zscore_", c("two", "less", "greater")), 
          paste0("zscore_adj_", c("two", "less", "greater")))
    
    return(list(num=perm_num, map=perm_map_res))
}


#### use function############
phe <- dat.splt$ko.pe
rownames(phe) <- phe[, 1]
phe <- phe[, -1, F]
pro <- dat.splt$ko.pf
test.res <- testFun(phe, pro)
if(type == "path"){
  dat.perm <- permFun(test.res, ampk)
}else{
  dat.perm <- permFun(test.res, amdk)
}
##### output ######################################################
name2 <- paste0(dir, "/", prefix,".step2.", type)

###1.test.result#######
file21 <- paste0(name2,".test.csv")
write.csv(test.res, file21, row.names=T, quote=F)

###2. p.value.hist
test.res.cln <- sapply(test.res, function(x){suppressWarnings(as.numeric(as.character(x)))})
rownames(test.res.cln) <- rownames(test.res)
dat.p <- test.res.cln[, 1]
lab <- hist(dat.p , plot=F, breaks=20)
Peak <- lab$mids[which.max(lab$density)]
sd <- sd(dat.p)
main <- paste0("p.value.Peak=",Peak," Sd=",round(sd,3))
pdf.out <- paste0(name2, ".P.value.hist.pdf")
pdf(pdf.out)
hist(dat.p, freq=F, breaks=20, xlab="p.value", main=main)
lines(density(dat.p))
dev.off()

###3. FDR  ####
dat.fdr <- test.res.cln[, c("p_value_twosides", "p_adjust_twosides")]
dat.fdr.cln <- dat.fdr[order(dat.fdr[, c("p_value_twosides")]), ]
fdr.1 <- dat.fdr.cln[max(which(dat.fdr.cln[, 1] < 0.01)), ]
fdr.5 <- dat.fdr.cln[max(which(dat.fdr.cln[, 1] < 0.05)), ]
fdr <- data.frame("Significant_level"=c(0.01,0.05), "FDR"=c(fdr.1[2], fdr.5[2]))
fdr.out <- paste0(name2,".FDR.csv")
write.csv(fdr, fdr.out, row.names=F, quote=F)

###4. permutation result ####
file24 <- paste0(name2, ".permutation.csv")
write.csv(dat.perm$num, file24, quote = F, row.names = F)

###5. reportscore result ####
file5 <- paste0(name2, ".reportscore.csv")
write.csv(dat.perm$map, file5, quote = F, row.names = F)