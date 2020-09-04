mypca <-  function(data, config, scale=F, top=5,pc1.var=1, pc2.var=2){
  # 是否scale
  data <- data[,colSums(data)!=0]
  if(scale==F){pc <- prcomp(data)}else{pc <- prcomp(data, scale. = T)}
  # 
  loading <- pc$rotation[,c(pc1.var, pc2.var)]
  roatL <- apply(loading,1,function(x){(x[1]^2+x[2]^2)})
  var <- names(sort(roatL, decreasing = T)[1:top])
  loading <- loading[var, ]
  datapc <- data.frame(varnames=rownames(loading), loading)
  # 解释度
  eig <- pc$sdev^2
  pc1 <- eig[pc1.var]/sum(eig)*100
  pc2 <- eig[pc2.var]/sum(eig)*100
  pc1 <- paste0("PC", pc1.var,"(",round(pc1,2),"%)")
  pc2 <- paste0("PC", pc2.var,"(",round(pc2,2),"%)")
  dat2 <- data.frame(pc$x[,c(pc1.var, pc2.var)])
  # 
  mult <- min(
    (max(dat2[,2]) - min(dat2[,2])/(max(datapc[,3])-min(datapc[,3]))),
    (max(dat2[,1]) - min(dat2[,1])/(max(datapc[,2])-min(datapc[,2])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * datapc[,2],
                      v2 = .7 * mult * datapc[,3]
  )
  #
  dat2$group <- as.factor(config[rownames(dat2),1])
  
  plot <- ggplot(dat2,aes(dat2[,1], dat2[,2],color=group))+geom_point()+theme_bw()+
    ggtitle("pca")+theme(plot.title = element_text(hjust = 0.5, size =14),
                        axis.title = element_text(size = 10),
                        axis.text = element_text(size = 13),
                        panel.grid = element_blank(),
                        legend.title = element_text(size = 15),
                         legend.text = element_text(size = 13))+xlab(pc1)+ylab(pc2)
  plot <- plot+geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), 
                         size = 4, vjust=1, color="black")
  plot <- plot+geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2),
                        arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black")
  
  return(plot)
  
}

### pcoa
mypcoa <- function(data, config, scale=F, method = "bray"){
  data <- data[,colSums(data)!=0]
  if(scale==F){
    dis <- vegdist(data, method = method)
    pco <- pcoa(dis)
  }else{
    data <- scale(data, center=F, scale=T)
    dis <- vegdist(data, method = method)
    pco <- pcoa(dis)
  }
  eig <- pco$value[,1]
  pc1 <- eig[1]/sum(eig)*100
  pc2 <- eig[2]/sum(eig)*100
  pc1 <- paste0("pcoa1(",round(pc1,2),"%)")
  pc2 <- paste0("pcoa2(",round(pc2,2),"%)")
  dat2 <- data.frame(pco$vector[,1:2])
  colnames(dat2) <- c("PC1", "PC2")
  dat2$group <- as.factor(config[rownames(dat2),1])
  plot <- ggplot(dat2,aes(PC1,PC2,color=group))+geom_point()+theme_bw()+ggtitle("pcoa")+theme(plot.title = element_text(hjust = 0.5,
                                                    size = 14),
                            axis.title = element_text(size = 10),
                            axis.text = element_text(size = 13),
                            panel.grid = element_blank(),
                           legend.title = element_text(size = 15),
                           legend.text = element_text(size = 13))+xlab(pc1)+ylab(pc2)
  return(plot)
}


# NMDS
myNMDS <- function(data, config, scale=F, method = "bray"){
  data <- data[,colSums(data)!=0]
  if(scale==F){
    metaMDS(data,k=2,trymax=100,  distance = method ) -> MDSfit
  }else{
    data <- scale(data, center=F, scale=T)
    metaMDS(data,k=2,trymax=100, distfun=betadiver, distance = method) -> MDSfit
  }
  
  dat2 <- as.data.frame(MDSfit$points)
  dat2$group <- as.factor(config[rownames(dat2),1])
  plot <- ggplot(dat2,aes(MDS1,MDS2,color=group))+geom_point()+theme_bw()+ggtitle("NMDS")+
                    theme(plot.title = element_text(hjust = 0.5,size = 14),
                    axis.title = element_text(size = 10),
                    axis.text = element_text(size = 13),
                    panel.grid = element_blank(),
                    legend.title = element_text(size = 15),
                    legend.text = element_text(size = 13))
  return(plot)
}




# t-SNE

mytsne <- function(data, config, scale=F, method="bray"){
  data <- data[,colSums(data)!=0]
  if(scale==F){
    dis <- vegdist(data, method = method)
    tsne <- Rtsne(dis, dims = 2, preplexity=30, verbose=T, max_iter = 500)
  }else{
    data <- scale(data, center=F, scale=T)
    dis <- vegdist(data, method = method)
    tsne <- Rtsne(dis, dims = 2, preplexity=30, verbose=T, max_iter = 500)
  }
  
  pc1 <- "t-SNE1"
  pc2 <- "t-SNE2"
  dat2 <- data.frame(tsne$Y)
  colnames(dat2) <- c("y1", "y2")
  rownames(dat2) <- rownames(data)
  dat2$group <- as.factor(config[rownames(dat2),1])
  plot <- ggplot(dat2,aes(y1,y2,color=group))+geom_point()+theme_bw()+ggtitle("tsne")+
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 13),
          panel.grid = element_blank(),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13))+xlab(pc1)+ylab(pc2)
  return(plot)
  
  
}

## get the core genus

core <- function(dat){
  sum <- apply(dat, 1, sum)
  index <- 1
  for (i in 1:length(sum)) {
    sum.sort <- sort(sum, decreasing = T) 
    if(sum(sum.sort[1:i])/sum(sum)>0.99){
      index <- i
      break
    }
  }
  return(names(sum.sort)[1:index])
}


### myCCA
mypca <-  function(data, config, scale=F, top=5,pc1.var=1, pc2.var=2){
  # 是否scale
  data <- data[,colSums(data)!=0]
  if(scale==F){pc <- prcomp(data)}else{pc <- prcomp(data, scale. = T)}
  # 
  loading <- pc$rotation[,c(pc1.var, pc2.var)]
  roatL <- apply(loading,1,function(x){(x[1]^2+x[2]^2)})
  var <- names(sort(roatL, decreasing = T)[1:top])
  loading <- loading[var, ]
  datapc <- data.frame(varnames=rownames(loading), loading)
  # 解释度
  eig <- pc$sdev^2
  pc1 <- eig[pc1.var]/sum(eig)*100
  pc2 <- eig[pc2.var]/sum(eig)*100
  pc1 <- paste0("PC", pc1.var,"(",round(pc1,2),"%)")
  pc2 <- paste0("PC", pc2.var,"(",round(pc2,2),"%)")
  dat2 <- data.frame(pc$x[,c(pc1.var, pc2.var)])
  # 
  mult <- min(
    (max(dat2[,2]) - min(dat2[,2])/(max(datapc[,3])-min(datapc[,3]))),
    (max(dat2[,1]) - min(dat2[,1])/(max(datapc[,2])-min(datapc[,2])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * datapc[,2],
                      v2 = .7 * mult * datapc[,3]
  )
  #
  dat2$group <- as.factor(config[rownames(dat2),1])
  
  plot <- ggplot(dat2,aes(dat2[,1], dat2[,2],color=group))+geom_point()+theme_bw()+
    ggtitle("pca")+theme(plot.title = element_text(hjust = 0.5, size =14),
                         axis.title = element_text(size = 10),
                         axis.text = element_text(size = 13),
                         panel.grid = element_blank(),
                         legend.title = element_text(size = 15),
                         legend.text = element_text(size = 13))+xlab(pc1)+ylab(pc2)
  plot <- plot+geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), 
                         size = 4, vjust=1, color="black")
  plot <- plot+geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2),
                            arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black")
  
  return(plot)
  
}