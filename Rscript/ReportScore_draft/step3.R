## script contain 2 steps
### 1. prepare database
### 2. choose remain/delete

### map/module -> levels#####################################
preFun <- function(x, y){
  res1 <- NULL
  for(i in 1:nrow(x)){
    dat1 <- as.character(x[i, 1])
    map1 <- rownames(x)[i]
    tmp1 <- unlist(strsplit(dat1, ";"))
    res1 <- data.frame(rbind(res1, cbind(Map=map1, level1=tmp1[1],
              level2=tmp1[2],Description_map=tmp1[3])))
  }
  rownames(res1) <- res1$Map
  res2 <- NULL
  for(i in 1:nrow(info)){
    dat2 <- as.character(y[i, 3])
    tmp2 <- substr(dat2, 0, 8)
    chr1 <- as.character(y[i, 1])
    chr2 <- as.character(y[i, 2])
    res2 <- rbind(res2, cbind(Module=chr1, Description_module=chr2, Map=tmp2))
  }
  mdu <- merge(res1, res2, by="Map")
  rownames(mdu) <- mdu$Module
  return(list(path=res1[, -1], modu=mdu[, -c(1, 5)]))
}

###intersection###########
inteFun <- function(x, y, z, m){
  name <- intersect(x$ID, intersect(rownames(y), z$V1))
  por.cln <- x[x$ID%in%name, ]
  pre.cln <- y[rownames(y)%in%name, ]
  cov.cln <- m[rownames(m)%in%name, ]
  por.cln <- por.cln[order(por.cln$ID), ]
  pre.cln <- pre.cln[order(rownames(pre.cln)), ]
  cov.cln <- cov.cln[order(rownames(cov.cln)), ]
  dat <- cbind(por.cln, pre.cln, cov.cln)
  return(dat)
}

###select reportscore########################################
selFun <- function(x){
  res <- data.frame()
  for(i in 1:nrow(x)){
    dat1 <- as.numeric(x[i, 7])
    dat2 <- as.numeric(x[i, 8])
    dat3 <- as.numeric(x[i, 2])
    dat4 <- as.numeric(x[i, 3])
    dat5 <- as.numeric(x[i, 4])
    if(dat1 > 0){
      score <- dat1
    }else{
      score <- -dat2
    }
    tmp1 <- round((dat4/dat3), 2)
    tmp2 <- round((dat5/dat3), 2)
    if(tmp1>0.4 & tmp2 >0.4){
      tmp1.res <- "Remain"
    }else{
      tmp1.res <- "Delete"
    }
    chr1 <- as.character(x[i, 1])
    tmp <- cbind("ID"=chr1, x[i, 12:17], "Reportscore"=score, "State"=tmp1.res)
    res <- rbind(res, tmp)
  }
  return(res)
}

#### use function ####################
por <- dat.perm$map
cov <- dat.cover
pre <- preFun(mpki, info)
if(type=="path"){
  int.res <- inteFun(por, pre$path, mrp, cov)
}else{
  int.res <- inteFun(por, pre$modu, mrm, cov)
}
dat.score <- selFun(int.res)

####   output  ###########
name3 <- paste0(dir, "/", prefix,".step3.", type)
file31 <- paste0(name3, ".result.txt")
write.table(dat.score, file31, quote=F, row.names=F, sep = "\t")