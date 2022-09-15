library(rsatscan)
td = "E:\\zuomian\\DDW-RF\\cluster\\"
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
pm2.5$lgpm2.5 <- log10(pm2.5$pm2.5)
data1 <- split(pm2.5,pm2.5$time)
s <- data.frame()               
for (i in 1:365) {
  s <- data1[[i]][,c(2,29)]
  s$case <- 1
  s <- cbind(s[,1],s[,3],s[,2])
  colnames(s) <- c("field","case","lgpm2.5")
  s <- as.data.frame(s)
  f <- paste(i,sep = "")
  write.cas(s,"E:\\zuomian\\DDW-RF\\cluster\\",f)
  print(i)
}
library(readxl)
zd <- read_excel("E:\\zuomian\\zd.xlsx")
zd <- as.data.frame(zd) 
write.geo(zd, td, "zd") 
pm2.5$lgpm2.5 <- log10(pm2.5$pm2.5)
data1 <- split(pm2.5,pm2.5$time)
data <- list.files("E:\\zuomian\\cluster\\",pattern='.cas')
mean <- read.csv('E:\\zuomian\\mean.csv')
sd <- read.csv('E:\\zuomian\\sd.csv')
llrmcs <- read.csv('E:\\zuomian\\llrmcs.csv') 
N <- 102
{
  llrmcs <- data.frame()
  for (i in 1:365) {
    miu <- mean(data1[[i]]$lgpm2.5)
    theta <- sd(data1[[i]]$lgpm2.5)
    xz <- subset(data1[[i]]$lgpm2.5,data1[[i]]$lgpm2.5>miu)
    x_z <- subset(data1[[i]]$lgpm2.5,data1[[i]]$lgpm2.5<=miu)
    labmda <- mean(x_z)
    f1 <- 1/102*(sum(xz^2)-2*sum(xz)*mean(xz)+length(xz)*mean(xz)^2+
                  sum(x_z^2)-2*sum(x_z)*labmda+length(x_z)*labmda^2)
    llrmcs[i,1] <- 102*log(theta)+sum(((data1[[i]]$lgpm2.5-miu)^2)/(2*(theta^2)))-102*log(sqrt(f1))-102/2##llr
  }
}
site <- read_excel("E:\\zuomian\\zd.xlsx") #注意位置信息lat在前
region.id <- as.character(site$field)
sp::coordinates(site)<-~long+lat
site <- sp::SpatialPoints(site)
sp::proj4string(site)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
kn <- spdep::knearneigh(site, k = 3)
nblist <- spdep::knn2nb(kn,row.names = region.id) # 邻近单元list
mcs <- data.frame()
mcsp <- list()
clust <- list()
timestart<-Sys.time()
for (j in 1:length(data)){
  a <- paste(j,".cas",sep = "")
  filename <- a
  pmmean <- mean[j,3]
  pmsd <- sd[j,3]
  pmllrmcs <- llrmcs[j,3]
  f <- log(pmsd)
  lnL0 <- -N*log(sqrt(2*pi))-N*log(pmsd)-sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))
  invisible(ss.options(reset=TRUE))
  ss.options(list(CaseFile=filename,CoordinatesFile="zd.geo",PrecisionCaseTimes=0,
                  CoordinatesType=1,AnalysisType=1,ModelType=5,ScanAreas=1,
                  CriteriaForReportingSecondaryClusters=0,ReportGiniClusters="n",
                  MaxSpatialSizeInDistanceFromCenter=1,SpatialWindowShapeType=0,
                  NonCompactnessPenalty=1))
  mcs_p <- sapply(c(1:50),function(i){
    ss.options(list(MaxSpatialSizeInPopulationAtRisk=i))
    write.ss.prm(td,substring(filename,1,1))
    cluster = satscan(td, substring(filename,1,1), sslocation="D:\\sat\\") ##the location of rsatscan
    b <- cluster$gis[which(cluster$gis$P_VALUE<0.05),]
    a <- split(b,b$CLUSTER)
    if (length(a)<2){
      xz <- b$LOC_MEAN
      x_z <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% b$LOC_ID)]
      f2 <- 1/N*((sum(xz^2)-2*sum(xz)*mean(xz)+length(xz)*mean(xz)^2)+
                   (sum(x_z^2)-2*sum(x_z)*mean(x_z)+length(x_z)*mean(x_z)^2))
      llr <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f2))-N*1/2
      return(llr/pmllrmcs)
    }else if(length(a)==2){
      neigh1 <- match(b[which(b$CLUSTER==1),1],region.id)##ID of cluster1
      neigh2 <- match(b[which(b$CLUSTER==2),1],region.id)##ID of cluster2
      neighber <- unique(unlist(nblist[neigh1]))####neighber of cluster1  1的邻近与2
      if(length(unique(c(neigh2,neighber)))==(length(neigh2)+length(neighber))){##判断为不相邻，分别计算
        xz1 <- a[[1]]$LOC_MEAN
        xz2 <- a[[2]]$LOC_MEAN
        x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID)]
        x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[2]]$LOC_ID)]
        f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                      (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
        f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                      (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
        llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
        llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
        llr <- llr1+llr2+lnL0
        return(llr/pmllrmcs)
      }else{
        xz <- b$LOC_MEAN##相邻，全部计算
        x_z <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% b$LOC_ID)]
        f2 <- 1/N*((sum(xz^2)-2*sum(xz)*mean(xz)+length(xz)*mean(xz)^2)+
                     (sum(x_z^2)-2*sum(x_z)*mean(x_z)+length(x_z)*mean(x_z)^2))
        llr <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f2))-N*1/2
        return(llr/pmllrmcs)
      }
    }else if(length(a)==3){
      neigh1 <- match(b[which(b$CLUSTER==1),1],region.id)##ID of cluster1
      neigh2 <- match(b[which(b$CLUSTER==2),1],region.id)##ID of cluster2
      neigh3 <- match(b[which(b$CLUSTER==3),1],region.id)##ID of cluster3
      neighber1 <- unique(unlist(nblist[neigh1]))####neighber of cluster1
      neighber2 <- unique(unlist(nblist[neigh2]))####neighber of cluster11的邻近与2
      neighber3 <- unique(unlist(nblist[neigh3]))####neighber of cluster11的邻近与2
      nei5 <- c(neigh2,neigh3)
      neighber5 <- unique(unlist(nblist[nei5]))
      nei6 <- c(neigh1,neigh3)
      neighber6 <- unique(unlist(nblist[nei6]))
      nei7 <- c(neigh1,neigh2)
      neighber7 <- unique(unlist(nblist[nei7]))
      if(length(unique(c(neigh2,neighber1)))==(length(neigh2)+length(neighber1))
         &length(unique(c(neigh3,neighber1)))==(length(neigh3)+length(neighber1))
         &length(unique(c(neigh3,neighber2)))==(length(neigh3)+length(neighber2))){##判断为不相邻，分别计算
        xz1 <- a[[1]]$LOC_MEAN
        xz2 <- a[[2]]$LOC_MEAN
        xz3 <- a[[3]]$LOC_MEAN
        x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID)]
        x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[2]]$LOC_ID)]
        x_z3 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[3]]$LOC_ID)]
        f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                      (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
        f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                      (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
        f32 <- 1/N*((sum(xz3^2)-2*sum(xz3)*mean(xz3)+length(xz3)*mean(xz3)^2)+
                      (sum(x_z3^2)-2*sum(x_z3)*mean(x_z3)+length(x_z3)*mean(x_z3)^2))
        llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
        llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
        llr3 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f32))-N*1/2
        llr <- llr1+llr2+llr3+2*lnL0
        return(llr/pmllrmcs)}
      else if(length(unique(c(neigh1,neighber5)))==(length(neigh1)+length(neighber5))){ #T
        xz1 <- a[[1]]$LOC_MEAN
        xz2 <- c(a[[2]]$LOC_MEAN,a[[3]]$LOC_MEAN)
        x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID)]
        x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[2]]$LOC_ID &!data1[[j]]$field %in% a[[3]]$LOC_ID)]
        f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                      (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
        f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                      (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
        llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
        llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
        llr <- llr1+llr2+lnL0
        return(llr/pmllrmcs)
      } else if(length(unique(c(neigh2,neighber6)))==(length(neigh2)+length(neighber6))){ #T
        xz1 <- a[[2]]$LOC_MEAN
        xz2 <- c(a[[1]]$LOC_MEAN,a[[3]]$LOC_MEAN)
        x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[2]]$LOC_ID)]
        x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID &!data1[[j]]$field %in% a[[3]]$LOC_ID)]
        f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                      (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
        f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                      (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
        llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
        llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
        llr <- llr1+llr2+lnL0
        return(llr/pmllrmcs)
      }else if(length(unique(c(neigh3,neighber7)))==(length(neigh3)+length(neighber7))){ #T
        xz1 <- a[[3]]$LOC_MEAN
        xz2 <- c(a[[1]]$LOC_MEAN,a[[2]]$LOC_MEAN)
        x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[3]]$LOC_ID)]
        x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID &!data1[[j]]$field %in% a[[2]]$LOC_ID)]
        f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                      (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
        f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                      (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
        llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
        llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
        llr <- llr1+llr2+lnL0
        return(llr/pmllrmcs)
      }
      else{
        return('3neighber')
      }
    }else if(length(a)==4){ 
      xz1 <- c(a[[1]]$LOC_MEAN,a[[3]]$LOC_MEAN)
      xz2 <- c(a[[2]]$LOC_MEAN,a[[4]]$LOC_MEAN)
      x_z1 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[1]]$LOC_ID & !data1[[j]]$field %in% a[[3]]$LOC_ID)]
      x_z2 <- data1[[j]]$lgpm2.5[which(!data1[[j]]$field %in% a[[2]]$LOC_ID & !data1[[j]]$field %in% a[[4]]$LOC_ID)]
      f12 <- 1/N*((sum(xz1^2)-2*sum(xz1)*mean(xz1)+length(xz1)*mean(xz1)^2)+
                    (sum(x_z1^2)-2*sum(x_z1)*mean(x_z1)+length(x_z1)*mean(x_z1)^2))
      f22 <- 1/N*((sum(xz2^2)-2*sum(xz2)*mean(xz2)+length(xz2)*mean(xz2)^2)+
                    (sum(x_z2^2)-2*sum(x_z2)*mean(x_z2)+length(x_z2)*mean(x_z2)^2))
      llr1 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f12))-N*1/2
      llr2 <- N*f+sum(((data1[[j]]$lgpm2.5-pmmean)^2)/(2*(pmsd^2)))-N*log(sqrt(f22))-N*1/2
      llr <- llr1+llr2+lnL0
      return(llr/pmllrmcs)
    }
  })
  mcs_p <- as.data.frame(mcs_p)
  colnames(mcs_p) <- 'mcs_p'
  window <- 1:50
  mcs_p <- cbind(window,mcs_p)
  mcsp[[j]] <- mcs_p
  mcs[j,1] <- j
  mcs[j,2] <- min(mcs_p$window[which(mcs_p$mcs_p==max(mcs_p$mcs_p,na.rm = T))])
  print(j)
}
timeend<-Sys.time()
timeend-timestart



pm2.5 <- read_excel("E:\\pm2.5.csv")
read_excel("E:\\zuomian\\zd.xlsx")
data <- list.files("E:\\zuomian\\cluster\\",pattern='.cas')
cc <- list()
timestart<-Sys.time()
region.id <- as.character(zd$field)
sp::coordinates(zd)<-~long+lat
site <- sp::SpatialPoints(zd)
sp::proj4string(site)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
kn <- spdep::knearneigh(site, k = 3)
nblist <- spdep::knn2nb(kn,row.names = region.id) # ?
for (j in 1:365){
  a <- paste(j,".cas",sep = "")
  filename <- a
  i <- mcs[j,2]
  invisible(ss.options(reset=TRUE)) 
  ss.options(list(CaseFile=filename,CoordinatesFile="zd.geo",PrecisionCaseTimes=0,
                  CoordinatesType=1,AnalysisType=1,ModelType=5,ScanAreas=1,
                  CriteriaForReportingSecondaryClusters=0,ReportGiniClusters="n",
                  MaxSpatialSizeInDistanceFromCenter=1,SpatialWindowShapeType=0,
                  NonCompactnessPenalty=1))
  ss.options(list(MaxSpatialSizeInPopulationAtRisk=i))
  write.ss.prm(td,substring(filename,1,1))
  cluster = satscan(td, substring(filename,1,1), sslocation="D:\\sat\\")
  b <- cluster$gis[which(cluster$gis$P_VALUE<0.05),]
  a <- split(b,b$CLUSTER)
  if(length(a)==1){
    Cluster <- b[,1:2]
  }else{
    neigh1 <- match(b[which(b$CLUSTER==1),1],region.id)##ID of cluster1
    neigh2 <- match(b[which(b$CLUSTER==2),1],region.id)##ID of cluster2
    neighber <- unique(unlist(nblist[neigh1]))####neighber of cluster1
    if(length(unique(c(neigh2,neighber)))==(length(neigh2)+length(neighber))){
      Cluster <- b[,1:2]
    }else{
      Cluster <- cbind(b[,1],1)
    }
  }
  cc[[j]] <- Cluster
  print(j)
}
save(cc,file = "cc.Rdata")
