library("openxlsx");library("broom");library("car");library("sp") ;library("maptools");library("spdep")
library("dplyr");library("plyr");library("ggplot2");library(readxl); library(geosphere);library(ranger)
dat_site <- read_excel("E:\\zuomian\\cyzd.xlsx") #
dat_site <- dat_site[order(dat_site$field),] # dat_site order by field
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
###site-based cross-validation
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)  
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  
  return(cvlist)
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))

###rf###########
cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  test <- do.call('rbind',data.test)
  model.rf<- ranger(pm2.5 ~d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi,data=train,num.trees = 400)
  prediction <- predict(model.rf,test)
  temp <- data.frame(cbind(test$pm2.5,prediction$predictions,test$time,test$field))
  temp
}
for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\rf\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
###r2/rme of RF
r2 <- data.frame()
rmse <- data.frame()
for (i in 1:10) {
  f <- paste("E:\\zuomian\\rf\\",i,".csv",sep="")
  result <- read.csv(f)
  rmse[i,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[i,1] <- 1-sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2)
  print(i)
}
mean(rmse$V1)
mean(r2$V1)

###idw
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))
cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  train.t <- split(train,train$time)
  zdpoint <- train.t[[1]]
  coordinates(zdpoint) <-~long+lat
  kn1 <- knearneigh(zdpoint, k=23) 
  region.id <- train.t[[1]]$field
  w_kn1 <- knn2nb(kn1, row.names=region.id) 
  w_kn1 <- knn2nb(kn1)
  map_crd <- coordinates(zdpoint)
  dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
  dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
  idw <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) #the idw matrix 
  c <- list() #calculate the autocorrelation terms
  for (q in 1:365) {
    c[[q]] <- data.frame()
    for (w in 1:nrow(train.t[[q]])) {
      a <- idw$neighbours[w]
      b <- idw$weights[w]
      IDs <- region.id
      c[[q]][w,1] <- IDs[w]
      c[[q]][w,2] <- train.t[[q]][a[[1]],]$pm2.5%*%b[[1]]
      colnames(c[[q]])[1] <- "field"
    }
    train.t[[q]] <- merge(train.t[[q]],c[[q]],by='field')
  }
  train.tt <- do.call("rbind",train.t)
  train.tt$sar <- as.numeric(train.tt$V2[,1])
  model.rf<- ranger(pm2.5 ~wind+d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi+sar,data=train.tt,num.trees = 500)  
  result <- matrix(nrow=length(data.test)*365,ncol=4)
  ###the  autocorrelation terms in test
  for (m in 1:length(data.test)) {
    a <- length(data.test)
    data.train[[a+1]] <- data.test[[m]]
    newtrain <- do.call("rbind",data.train)
    newtrain.t <- split(newtrain,newtrain$time)
    dat_site2 <- newtrain.t[[1]]
    region.id2 <- dat_site2$field
    coordinates(dat_site2) <-~long+lat
    kn1 <- knearneigh(dat_site2, k=23) 
    w_kn1 <- knn2nb(kn1, row.names=region.id2) 
    w_kn1 <- knn2nb(kn1)
    map_crd <- coordinates(dat_site2)
    dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
    dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
    weight3 <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
    c <- list() 
    for (q in 1:365) {
      c[[q]] <- data.frame()
      for (w in 1:nrow(newtrain.t[[q]])) {
        a <- weight3$neighbours[w]
        b <- weight3$weights[w]
        IDs <- region.id2
        c[[q]][w,1] <- IDs[w]
        c[[q]][w,2] <- newtrain.t[[q]][a[[1]],]$pm2.5%*%b[[1]]
        colnames(c[[q]])[1] <- "field"
      }
      newtrain.t[[q]] <- merge(newtrain.t[[q]],c[[q]],by='field')
    }
    newtrain.tt <- do.call("rbind",newtrain.t)
    newtrain.tt$sar <- as.numeric(newtrain.tt$V2[,1])
    nt <- split(newtrain.tt,newtrain.tt$field)
    a <- nt[[data.test[[m]]$field[1]]]
    prediction <- predict(model.rf,a)
    c1 <- (m-1)*365+1
    c2 <- m*365
    result[c1:c2,1] <- a$pm2.5
    result[c1:c2,2] <- prediction$predictions
    result[c1:c2,3] <- a$time
    result[c1:c2,4] <- a$field
  }
  result
}
for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\idw\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
###r2/rme of IDW-RF
r2 <- data.frame()
rmse <- data.frame()
for (i in 1:10) {
  f <- paste("E:\\zuomian\\idw\\",i,".csv",sep="")
  result <- read.csv(f)
  rmse[i,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[i,1] <- 1-sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2)
  print(i)
}
mean(rmse$V1)
mean(r2$V1)

########################DDW-RF###############
#####part0: common setting ####################################################
########W00-RF###############################
library(readxl);library(spdep);library(dplyr);library(spatialreg)
`%+%` <- function(x,y) paste0(x,y) 
`%>%` <- magrittr::`%>%` 
path <- "E:\\zuomian\\DDW0.1\\"  #the Founction of DDW0.1
load(file = "E:\\zuomian\\cc.Rdata")  ##the cluster each day in 2019
funpath <- path%+%"Rfun\\"
invisible(lapply(funpath%+%dir(funpath),source))  # loading function
datpath <- path%+%"data\\"
dat_site <- read_excel("E:\\zuomian\\cyzd.xlsx")
dat_site <- dat_site[order(dat_site$field),] 
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x]) 
  return(cvlist)
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))
cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  train.t <- split(train,train$time)
  dat_site1 <- train.t[[1]]
  region.id <- dat_site1$field
  for (l in 1:365) {
    dat_cluster <- as.data.frame(cc[[l]])
    if(nrow(dat_cluster)!=0){
      dat_site1 <- train.t[[1]]
      dat_cluster1 <- dat_cluster[dat_cluster$field%in%dat_site1$field,]
      clusterID <- as.character(dat_cluster1$field)
      region.id <- as.character(dat_site1$field)
      cluster <- as.character(dat_cluster1$cluster)
      nbNN <- cluster2nb(clusterID = clusterID,cluster = cluster,region.id = region.id)
      ddwlist <- spdep::nb2listw(neighbours = nbNN,zero.policy = T)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        if(length(a[[1]])>1){
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }else{
          zdpoint <- train.t[[1]]
          zdpoint <- zdpoint[order(zdpoint$field),]
          coordinates(zdpoint) <-~long+lat 
          kn1 <- knearneigh(zdpoint, k=23) 
          region.id <- train.t[[1]]$field
          w_kn1 <- knn2nb(kn1, row.names=region.id) 
          w_kn1 <- knn2nb(kn1)
          map_crd <- coordinates(zdpoint)
          dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
          dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
          ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
          a <- ddwlist$neighbours[w]
          b <- ddwlist$weights[w]
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }}}
    else{
      zdpoint <-train.t[[1]]
      zdpoint <- zdpoint[order(zdpoint$field),]
      coordinates(zdpoint) <-~long+lat 
      kn1 <- knearneigh(zdpoint, k=23) 
      region.id <- train.t[[1]]$field
      w_kn1 <- knn2nb(kn1, row.names=region.id) 
      w_kn1 <- knn2nb(kn1)
      map_crd <- coordinates(zdpoint)
      dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
      dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
      ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        train.t[[l]]$sar[w] <-train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
      } 
    }}
  train.tt <- do.call("rbind",train.t)
  train.tt <- na.omit(train.tt)
  model.rf<- ranger(pm2.5 ~d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi+sar,data=train.tt,num.trees = 400)
  result <- matrix(nrow=length(data.test)*365,ncol=4)
  for (m in 1:length(data.test)) {
    a <- length(data.train)
    data.train[[a+1]] <- data.test[[m]]
    newtrain <- do.call("rbind",data.train)
    newtrain.t <- split(newtrain,newtrain$time)
    dat_site2 <- newtrain.t[[1]]
    region.id2 <- dat_site2$field
    for (l in 1:365) {
      dat_cluster <- as.data.frame(cc[[l]])
      if(nrow(dat_cluster)!=0){
        dat_cluster2 <- dat_cluster[dat_cluster$field%in%dat_site2$field,]
        clusterID2 <- as.character(dat_cluster2$field)
        dat_site2 <- newtrain.t[[1]]
        region.id2 <- as.character(dat_site2$field)
        cluster <- as.character(dat_cluster2$cluster)
        nbNN <- cluster2nb(clusterID = clusterID2,cluster = cluster,region.id = region.id2)
        ddwlist2 <- spdep::nb2listw(neighbours = nbNN,zero.policy = T)
        for (w in 1:nrow(newtrain.t[[l]])) {
          a <- ddwlist2$neighbours[w]
          b <- ddwlist2$weights[w]
          if(length(a[[1]])>1){
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }else{
            zdpoint <- newtrain.t[[1]]
            zdpoint <- zdpoint[order(zdpoint$field),]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T)
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) #最终反距离加权后的权重矩阵 
            a <- ddwlist$neighbours[w]
            b <- ddwlist$weights[w]
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }}}else{
            zdpoint <- newtrain.t[[1]]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist2 <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
            for (w in 1:nrow(newtrain.t[[l]])) {
              a <- ddwlist2$neighbours[w]
              b <- ddwlist2$weights[w]
              newtrain.t[[l]]$sar[w] <-  newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
            } 
          }
    }
    newtrain.tt <- do.call("rbind",newtrain.t)
    nt <- split(newtrain.tt,newtrain.tt$field)
    a <- nt[[data.test[[m]]$field[1]]]
    prediction <- predict(model.rf,a)
    c1 <- (m-1)*365+1
    c2 <- m*365
    result[c1:c2,1] <- a$pm2.5
    result[c1:c2,2] <- prediction$predictions
    result[c1:c2,3] <- a$time
    result[c1:c2,4] <- a$field
  }
  result
}
for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\w00\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
###r2/rmse of W00-RF
r2 <- data.frame()
rmse <- data.frame()
r21 <- data.frame()
for (j in 1:10) {
  prediction1 <- cvtest(j)
  result <- prediction1
  rmse[j,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[j,1] <- (cor(result[,2],result[,1]))^2
  r21[j,1] <- 1-(sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2))
  print(j)
}
#the RMSE of W00-RF inside and outside of clusters
rmse <- data.frame()
for (k in 1:10) {
  f <- paste("E:\\zuomian\\w00\\",k,".csv",sep="")
  id <- read.csv(f)
  data <- split(id,id$time)
  da <- list()
  for (i in 1:365) {
    da[[i]] <- data[[i]][data[[i]]$field%in%cc[[i]]$field,]
  }
  result <- do.call("rbind",da)
  rmse[k,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  print(k)
}
mean(rmse$V1)


####W0d-RF#################
dat_site <- read_excel("E:\\zuomian\\cyzd.xlsx")
dat_site <- dat_site[order(dat_site$field),] 
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x]) 
  return(cvlist)
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))
cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  train.t <- split(train,train$time)
  dat_site1 <- train.t[[1]]
  region.id <- dat_site1$field
  sp::coordinates(dat_site1)<-~long+lat
  wuranpoint<-sp::SpatialPoints(dat_site1)
  sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
  kn3 <- spdep::knearneigh(wuranpoint, k = 23) 
  nblist <- spdep::knn2nb(kn3,row.names = region.id) 
  for (l in 1:365) {
    dat_cluster <- as.data.frame(cc[[l]])
    dat_site1 <- train.t[[1]]
    dat_cluster1 <- dat_cluster[dat_cluster$field%in%dat_site1$field,]
    if(nrow(dat_cluster)!=0&nrow(dat_cluster1)!=0){
      clusterID <- as.character(dat_cluster1$field)
      clusterID <- as.character(dat_cluster1$field)
      region.id <- as.character(dat_site1$field)
      cluster <- as.character(dat_cluster1$cluster)
      nbNG <- NeighberNG(clusterID = clusterID,cluster = cluster,gal = nblist)
      sp::coordinates(dat_site1)<-~long+lat
      dlistNG <- spdep::nbdists(nbNG, dat_site1@coords) %>% lapply(function(x) 1/x^2)
      dlistbase <- lapply(dlistNG[!region.id%in%clusterID],function(x) rep(1,length(x)))
      dlistNG[!region.id%in%clusterID] <- dlistbase 
      class(nbNG) <- c("nb","nbNG")
      attr(nbNG, "region.id") <- region.id
      attr(nbNG, "gis") <- TRUE
      attr(nbNG, "call") <- TRUE
      nbNG <- sym.attr.nb(nbNG)
      ddwlist <- spdep::nb2listw(neighbours = nbNG,glist = dlistNG,zero.policy = T)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        if(length(a[[1]])>1){
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }else{
          zdpoint <- train.t[[1]]
          zdpoint <- zdpoint[order(zdpoint$field),]
          coordinates(zdpoint) <-~long+lat 
          kn1 <- knearneigh(zdpoint, k=23) 
          region.id <- train.t[[1]]$field
          w_kn1 <- knn2nb(kn1, row.names=region.id) 
          w_kn1 <- knn2nb(kn1)
          map_crd <- coordinates(zdpoint)
          dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
          dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
          ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
          a <- ddwlist$neighbours[w]
          b <- ddwlist$weights[w]
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }}}
    else{
      zdpoint <-train.t[[1]]
      zdpoint <- zdpoint[order(zdpoint$field),]
      coordinates(zdpoint) <-~long+lat 
      kn1 <- knearneigh(zdpoint, k=23) 
      region.id <- train.t[[1]]$field
      w_kn1 <- knn2nb(kn1, row.names=region.id) 
      w_kn1 <- knn2nb(kn1)
      map_crd <- coordinates(zdpoint)
      dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
      dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
      ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        train.t[[l]]$sar[w] <-train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
      } 
    }}
  train.tt <- do.call("rbind",train.t)
  train.tt <- na.omit(train.tt)
  model.rf<- ranger(pm2.5 ~d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi+sar,data=train.tt,num.trees = 400)
  result <- matrix(nrow=length(data.test)*365,ncol=4)
  for (m in 1:length(data.test)) {
    a <- length(data.train)
    data.train[[a+1]] <- data.test[[m]]
    newtrain <- do.call("rbind",data.train)
    newtrain.t <- split(newtrain,newtrain$time)
    dat_site2 <- newtrain.t[[1]]
    region.id2 <- dat_site2$field
    sp::coordinates(dat_site2)<-~long+lat
    wuranpoint<-sp::SpatialPoints(dat_site2)
    sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
    kn3 <- spdep::knearneigh(wuranpoint, k = 23) 
    nblist <- spdep::knn2nb(kn3,row.names = region.id2)
    for (l in 1:365) {
      dat_cluster <- as.data.frame(cc[[l]])
      dat_cluster2 <- dat_cluster[dat_cluster$field%in%dat_site2$field,]
      if(nrow(dat_cluster)!=0&nrow(dat_cluster2)!=0){
        clusterID2 <- as.character(dat_cluster2$field)
        dat_site2 <- newtrain.t[[1]]
        region.id2 <- as.character(dat_site2$field)
        cluster <- as.character(dat_cluster2$cluster)
        sp::coordinates(dat_site2)<-~long+lat
        nbNG <- NeighberNG(clusterID = clusterID2,cluster = cluster,gal = nblist)
        dlistNG <- spdep::nbdists(nbNG, dat_site2@coords) %>% lapply(function(x) 1/x^2)
        dlistbase <- lapply(dlistNG[!region.id2%in%clusterID2],function(x) rep(1,length(x)))
        dlistNG[!region.id2%in%clusterID2] <- dlistbase 
        class(nbNG) <- c("nb","nbNG")
        attr(nbNG, "region.id") <- region.id2
        attr(nbNG, "gis") <- TRUE
        attr(nbNG, "call") <- TRUE
        nbNG <- sym.attr.nb(nbNG)
        ddwlist2 <- spdep::nb2listw(neighbours = nbNG,glist = dlistNG,zero.policy = T)
        for (w in 1:nrow(newtrain.t[[l]])) {
          a <- ddwlist2$neighbours[w]
          b <- ddwlist2$weights[w]
          if(length(a[[1]])>1){
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }else{
            zdpoint <- newtrain.t[[1]]
            zdpoint <- zdpoint[order(zdpoint$field),]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
            a <- ddwlist$neighbours[w]
            b <- ddwlist$weights[w]
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }}}else{
            zdpoint <- newtrain.t[[1]]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T)
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist2 <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
            for (w in 1:nrow(newtrain.t[[l]])) {
              a <- ddwlist2$neighbours[w]
              b <- ddwlist2$weights[w]
              newtrain.t[[l]]$sar[w] <-  newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
            } 
          }
    }
    newtrain.tt <- do.call("rbind",newtrain.t)
    nt <- split(newtrain.tt,newtrain.tt$field)
    a <- nt[[data.test[[m]]$field[1]]]
    prediction <- predict(model.rf,a)
    c1 <- (m-1)*365+1
    c2 <- m*365
    result[c1:c2,1] <- a$pm2.5
    result[c1:c2,2] <- prediction$predictions
    result[c1:c2,3] <- a$time
    result[c1:c2,4] <- a$field
  }
  result
}

for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\w0d\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
r2 <- data.frame()
rmse <- data.frame()
r2 <- data.frame()
for (j in 1:10) {
  prediction1 <- cvtest(j)
  result <- prediction1
  rmse[j,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[j,1] <- 1-(sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2))
  print(j)
}
mean(r2$V1)
mean(rmse$V1)

rmse <- data.frame()
for (k in 1:10) {
  f <- paste("E:\\zuomian\\w0d\\",k,".csv",sep="")
  id <- read.csv(f)
  data <- split(id,id$time)
  da <- list()
  for (i in 1:365) {
    da[[i]] <- data[[i]][data[[i]]$field%in%cc[[i]]$field,]
  }
  result <- do.call("rbind",da)
  rmse[k,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  print(k)
}
mean(rmse$V1)


#####Wd0-RF####################
dat_site <- read_excel("E:\\zuomian\\cyzd.xlsx")
dat_site <- dat_site[order(dat_site$field),] 
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x]) 
  return(cvlist)
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))
cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  train.t <- split(train,train$time)
  dat_site1 <- train.t[[1]]
  region.id <- dat_site1$field
  sp::coordinates(dat_site1)<-~long+lat
  wuranpoint<-sp::SpatialPoints(dat_site1)
  sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
  kn3 <- spdep::knearneigh(wuranpoint, k = 23) 
  nblist <- spdep::knn2nb(kn3,row.names = region.id) 
  for (l in 1:365) {
    dat_cluster <- as.data.frame(cc[[l]])
    dat_site1 <- train.t[[1]]
    dat_cluster1 <- dat_cluster[dat_cluster$field%in%dat_site1$field,]
    if(nrow(dat_cluster)!=0&nrow(dat_cluster1)!=0){
      clusterID <- as.character(dat_cluster1$field)
      region.id <- as.character(dat_site1$field)
      cluster <- as.character(dat_cluster1$cluster)
      nbGN <- NeighberGN(clusterID = clusterID,cluster = cluster,gal = nblist)
      sp::coordinates(dat_site1)<-~long+lat
      dlistGN <- spdep::nbdists(nbGN, dat_site1@coords) %>% lapply(function(x) 1/x^2)
      dlistcluster <- lapply(dlistGN[region.id%in%clusterID],function(x) rep(1,length(x)))
      dlistGN[region.id%in%clusterID] <- dlistcluster 
      ddwlist <- spdep::nb2listw(neighbours = nbGN,glist = dlistGN,zero.policy = T)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        if(length(a[[1]])>1){
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }else{
          zdpoint <- train.t[[1]]
          zdpoint <- zdpoint[order(zdpoint$field),]
          coordinates(zdpoint) <-~long+lat 
          kn1 <- knearneigh(zdpoint, k=23) 
          region.id <- train.t[[1]]$field
          w_kn1 <- knn2nb(kn1, row.names=region.id) 
          w_kn1 <- knn2nb(kn1)
          map_crd <- coordinates(zdpoint)
          dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
          dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
          ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
          a <- ddwlist$neighbours[w]
          b <- ddwlist$weights[w]
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }}}
    else{
      zdpoint <-train.t[[1]]
      zdpoint <- zdpoint[order(zdpoint$field),]
      coordinates(zdpoint) <-~long+lat 
      kn1 <- knearneigh(zdpoint, k=23) 
      region.id <- train.t[[1]]$field
      w_kn1 <- knn2nb(kn1, row.names=region.id) 
      w_kn1 <- knn2nb(kn1)
      map_crd <- coordinates(zdpoint)
      dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
      dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
      ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        train.t[[l]]$sar[w] <-train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
      } 
    }}
  train.tt <- do.call("rbind",train.t)
  train.tt <- na.omit(train.tt)
  model.rf<- ranger(pm2.5 ~d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi+sar,data=train.tt,num.trees = 400)
  result <- matrix(nrow=length(data.test)*365,ncol=4)
  for (m in 1:length(data.test)) {
    a <- length(data.train)
    data.train[[a+1]] <- data.test[[m]]
    newtrain <- do.call("rbind",data.train)
    newtrain.t <- split(newtrain,newtrain$time)
    dat_site2 <- newtrain.t[[1]]
    region.id2 <- dat_site2$field
    sp::coordinates(dat_site2)<-~long+lat
    wuranpoint<-sp::SpatialPoints(dat_site2)
    sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
    kn3 <- spdep::knearneigh(wuranpoint, k = 23) 
    nblist <- spdep::knn2nb(kn3,row.names = region.id2)
    for (l in 1:365) {
      dat_cluster <- as.data.frame(cc[[l]])
      dat_cluster2 <- dat_cluster[dat_cluster$field%in%dat_site2$field,]
      if(nrow(dat_cluster)!=0&nrow(dat_cluster2)!=0){
        clusterID2 <- as.character(dat_cluster2$field)
        dat_site2 <- newtrain.t[[1]]
        region.id2 <- as.character(dat_site2$field)
        cluster <- as.character(dat_cluster2$cluster)
        sp::coordinates(dat_site2)<-~long+lat
        nbGN <- NeighberGN(clusterID = clusterID2,cluster = cluster,gal = nblist)
        dlistGN <- spdep::nbdists(nbGN, dat_site2@coords) %>% lapply(function(x) 1/x^2)
        dlistcluster <- lapply(dlistGN[region.id2%in%clusterID2],function(x) rep(1,length(x)))
        dlistGN[region.id2%in%clusterID2] <- dlistcluster 
        ddwlist2 <- spdep::nb2listw(neighbours = nbGN,glist = dlistGN,zero.policy = T)
        for (w in 1:nrow(newtrain.t[[l]])) {
          a <- ddwlist2$neighbours[w]
          b <- ddwlist2$weights[w]
          if(length(a[[1]])>1){
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }else{
            zdpoint <- newtrain.t[[1]]
            zdpoint <- zdpoint[order(zdpoint$field),]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
            a <- ddwlist$neighbours[w]
            b <- ddwlist$weights[w]
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }}}else{
            zdpoint <- newtrain.t[[1]]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T)
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist2 <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
            for (w in 1:nrow(newtrain.t[[l]])) {
              a <- ddwlist2$neighbours[w]
              b <- ddwlist2$weights[w]
              newtrain.t[[l]]$sar[w] <-  newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
            } 
          }
    }
    newtrain.tt <- do.call("rbind",newtrain.t)
    nt <- split(newtrain.tt,newtrain.tt$field)
    a <- nt[[data.test[[m]]$field[1]]]
    prediction <- predict(model.rf,a)
    c1 <- (m-1)*365+1
    c2 <- m*365
    result[c1:c2,1] <- a$pm2.5
    result[c1:c2,2] <- prediction$predictions
    result[c1:c2,3] <- a$time
    result[c1:c2,4] <- a$field
  }
  result
}

for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\wdd\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
r2 <- data.frame()
rmse <- data.frame()
r2 <- data.frame()
for (j in 1:10) {
  prediction1 <- cvtest(j)
  result <- prediction1
  rmse[j,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[j,1] <- 1-(sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2))
  print(j)
}
mean(r2$V1)
mean(rmse$V1)

rmse <- data.frame()
for (k in 1:10) {
  f <- paste("E:\\zuomian\\wdd\\",k,".csv",sep="")
  id <- read.csv(f)
  data <- split(id,id$time)
  da <- list()
  for (i in 1:365) {
    da[[i]] <- data[[i]][data[[i]]$field%in%cc[[i]]$field,]
  }
  result <- do.call("rbind",da)
  rmse[k,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  print(k)
}
mean(rmse$V1)



####Wdd-RF###########
dat_site <- read_excel("E:\\zuomian\\cyzd.xlsx")
dat_site <- dat_site[order(dat_site$field),] 
pm2.5 <- read.csv("E:\\zuomian\\pm.csv")
data <- split(pm2.5,pm2.5$field)
CVgroup <- function(k,datasize){
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize] 
  temp <- sample(n,datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
set.seed(12365)
data.cv <- replicate(10,CVgroup(k=10,datasize = length(data)))

cvtest <- function(i){
  data.train <- data[-data.cv[[i]]]
  data.test <- data[data.cv[[i]]]
  train <- do.call('rbind',data.train)
  train.t <- split(train,train$time)
  dat_site1 <- train.t[[1]]
  region.id <- dat_site1$field
  sp::coordinates(dat_site1)<-~long+lat
  wuranpoint<-sp::SpatialPoints(dat_site1)
  sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
  kn3 <- spdep::knearneigh(wuranpoint, k = 23) 
  nblist <- spdep::knn2nb(kn3,row.names = region.id) 
  for (l in 1:365) {
    dat_cluster <- as.data.frame(cc[[l]])
    dat_site1 <- train.t[[1]]
    dat_cluster1 <- dat_cluster[dat_cluster$field%in%dat_site1$field,]
    if(nrow(dat_cluster)!=0&nrow(dat_cluster1)!=0){
      clusterID <- as.character(dat_cluster1$field)
      region.id <- as.character(dat_site1$field)
      cluster <- as.character(dat_cluster1$cluster)
      nbGG <- NeighberGG(clusterID = clusterID,cluster = cluster,gal = nblist) # neighbor list
      class(nbGG) <- c("nb","nbGG")
      attr(nbGG, "region.id") <- region.id
      attr(nbGG, "gis") <- TRUE
      attr(nbGG, "call") <- TRUE
      nbGG <- sym.attr.nb(nbGG)
      sp::coordinates(dat_site1)<-~long+lat
      dlist <- spdep::nbdists(nbGG, dat_site1@coords) %>% lapply(function(x) 1/x^2) 
      ddwlist <- spdep::nb2listw(neighbours = nbGG,glist = dlist,zero.policy = T)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        if(length(a[[1]])>1){
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }else{
          zdpoint <- train.t[[1]]
          zdpoint <- zdpoint[order(zdpoint$field),]
          coordinates(zdpoint) <-~long+lat 
          kn1 <- knearneigh(zdpoint, k=23) 
          region.id <- train.t[[1]]$field
          w_kn1 <- knn2nb(kn1, row.names=region.id) 
          w_kn1 <- knn2nb(kn1)
          map_crd <- coordinates(zdpoint)
          dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) #求距离
          dlist_p1 <- lapply(dlist, function(x) 1/(x*x))#幂函数反距离
          #反距离加权
          ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) #最终反距离加权后的权重矩阵 
          a <- ddwlist$neighbours[w]
          b <- ddwlist$weights[w]
          train.t[[l]]$sar[w] <- train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
        }}}
    else{
      zdpoint <-train.t[[1]]
      zdpoint <- zdpoint[order(zdpoint$field),]
      coordinates(zdpoint) <-~long+lat 
      kn1 <- knearneigh(zdpoint, k=23) 
      region.id <- train.t[[1]]$field
      w_kn1 <- knn2nb(kn1, row.names=region.id) 
      w_kn1 <- knn2nb(kn1)
      map_crd <- coordinates(zdpoint)
      dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
      dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
      ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
      for (w in 1:nrow(train.t[[l]])) {
        a <- ddwlist$neighbours[w]
        b <- ddwlist$weights[w]
        train.t[[l]]$sar[w] <-train.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
      } 
    }}
  train.tt <- do.call("rbind",train.t)
  train.tt <- na.omit(train.tt)
  model.rf<- ranger(pm2.5 ~wind+d2m+t2m+fal+blh+vis+tcc+sp+rh+tp+v10+u10+ndvi+sar,data=train.tt,num.trees = 500)
  result <- matrix(nrow=length(data.test)*365,ncol=4)
  for (m in 1:length(data.test)) {
    a <- length(data.train)
    data.train[[a+1]] <- data.test[[m]]
    newtrain <- do.call("rbind",data.train)
    newtrain.t <- split(newtrain,newtrain$time)
    dat_site2 <- newtrain.t[[1]]
    region.id2 <- dat_site2$field
    sp::coordinates(dat_site2)<-~long+lat
    wuranpoint<-sp::SpatialPoints(dat_site2)
    sp::proj4string(wuranpoint)<-"+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +datum=WGS84"
    kn3 <- spdep::knearneigh(wuranpoint, k = 23) # the three nearest sites
    nblist <- spdep::knn2nb(kn3,row.names = region.id2)
    for (l in 1:365) {
      dat_cluster <- as.data.frame(cc[[l]])
      dat_cluster2 <- dat_cluster[dat_cluster$field%in%dat_site2$field,]
      if(nrow(dat_cluster)!=0&nrow(dat_cluster2)!=0){
        clusterID2 <- as.character(dat_cluster2$field)
        dat_site2 <- newtrain.t[[1]]
        region.id2 <- as.character(dat_site2$field)
        cluster <- as.character(dat_cluster2$cluster)
        sp::coordinates(dat_site2)<-~long+lat
        nbGG <- NeighberGG(clusterID = clusterID2,cluster = cluster,gal = nblist) # neighbor list
        class(nbGG) <- c("nb","nbGG")
        attr(nbGG, "region.id") <- region.id2 
        attr(nbGG, "gis") <- TRUE
        attr(nbGG, "call") <- TRUE
        nbGG <- sym.attr.nb(nbGG)
        dlist <- spdep::nbdists(nbGG, dat_site2@coords) %>% lapply(function(x) 1/x^2) # 邻近单元反距离
        ddwlist2 <- spdep::nb2listw(neighbours = nbGG,glist = dlist,zero.policy = T)
        for (w in 1:nrow(newtrain.t[[l]])) {
          a <- ddwlist2$neighbours[w]
          b <- ddwlist2$weights[w]
          if(length(a[[1]])>1){
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }else{
            zdpoint <- newtrain.t[[1]]
            zdpoint <- zdpoint[order(zdpoint$field),]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE)
            a <- ddwlist$neighbours[w]
            b <- ddwlist$weights[w]
            newtrain.t[[l]]$sar[w] <- newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
          }}}else{
            zdpoint <- newtrain.t[[1]]
            coordinates(zdpoint) <-~long+lat 
            kn1 <- knearneigh(zdpoint, k=23) 
            region.id2 <- newtrain.t[[1]]$field
            w_kn1 <- knn2nb(kn1, row.names=region.id2) 
            w_kn1 <- knn2nb(kn1)
            map_crd <- coordinates(zdpoint)
            dlist <- nbdists(w_kn1, map_crd[,c(1,2)], longlat = T) 
            dlist_p1 <- lapply(dlist, function(x) 1/(x*x))
            ddwlist2 <- nb2listw(w_kn1, glist=dlist_p1, zero.policy=TRUE) 
            for (w in 1:nrow(newtrain.t[[l]])) {
              a <- ddwlist2$neighbours[w]
              b <- ddwlist2$weights[w]
              newtrain.t[[l]]$sar[w] <-  newtrain.t[[l]][a[[1]],]$pm2.5%*%b[[1]]
            } 
          }
    }
    newtrain.tt <- do.call("rbind",newtrain.t)
    nt <- split(newtrain.tt,newtrain.tt$field)
    a <- nt[[data.test[[m]]$field[1]]]
    prediction <- predict(model.rf,a)
    c1 <- (m-1)*365+1
    c2 <- m*365
    result[c1:c2,1] <- a$pm2.5
    result[c1:c2,2] <- prediction$predictions
    result[c1:c2,3] <- a$time
    result[c1:c2,4] <- a$field
  }
  result
}

for (j in 0:9) {
  prediction1 <- cvtest(j*10+1)
  for (i in (j*10+2):(j*10+10)) {
    prediction0 <- cvtest(i)
    prediction1 <- rbind(prediction1,prediction0)
  }
  result <- na.omit(prediction1)
  colnames(prediction1) <- c("pm2.5","pre","time","field")
  setwd("E:\\zuomian\\wdd\\")
  l <- j+1
  a <- paste(l,".csv",sep = "")
  write.csv(prediction1,a,row.names = F)
  print(j)
}
r2 <- data.frame()
rmse <- data.frame()
for (i in 1:10) {
  f <- paste("E:\\zuomian\\wdd\\",i,".csv",sep="")
  result <- read.csv(f)
  rmse[i,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  r2[i,1] <- 1-sum((result[,2]-result[,1])^2)/sum((mean(result[,1])-result[,1])^2)
  print(i)
}
mean(rmse$V1)
mean(r2$V1)

rmse <- data.frame()
for (k in 1:10) {
  f <- paste("E:\\zuomian\\w0d\\",k,".csv",sep="")
  id <- read.csv(f)
  data <- split(id,id$time)
  da <- list()
  for (i in 1:365) {
    da[[i]] <- data[[i]][data[[i]]$field%in%cc[[i]]$field,]
  }
  result <- do.call("rbind",da)
  rmse[k,1] <- sqrt(sum((result[,1]-result[,2])^2)/length(result[,1]))
  print(k)
}
mean(rmse$V1)

#############calculated the average prediction of the 10 cross-validation:RF\IDW-RF\DDW-RF###
data <- list()
for (i in 1:10) {
  f <- paste("E:\\zuomian\\rf\\",i,".csv",sep="")
  data[[i]] <- read.csv(f)
}
da <- data.frame()
pm <- as.data.frame(cbind(data[[1]]$pm2.5,data[[2]]$pm2.5,data[[3]]$pm2.5,data[[4]]$pm2.5,data[[5]]$pm2.5,data[[6]]$pm2.5data[[7]]$pm2.5data[[8]]$pm2.5data[[9]]$pm2.5data[[10]]$pm2.5))
pre <- as.data.frame(cbind(data[[1]]$pre,data[[2]]$pre,data[[3]]$pre,data[[4]]$pre,data[[5]]$pre,data[[6]]$predata[[7]]$predata[[8]]$predata[[9]]$predata[[10]]$pre))
mean <- data.frame()
for (i in 1:33945) {
  mean[i,1] <- (sum(pm[i,1:5]))/5
  mean[i,2] <- (sum(pre[i,1:5]))/5
}
colnames(mean) <- c("pm","pre") 
mean$time <- data[[1]]$time
mean$field <- data[[1]]$field
write.csv(mean,"ssidw.csv")


####plot
rfdata <- read.csv("E:\\zuomian\\sswrf.csv")
w00data <- read.csv("E:\\zuomian\\ssw00.csv")
wd0data <- read.csv("E:\\zuomian\\sswd0.csv")
w0ddata <- read.csv("E:\\zuomian\\ssw0d.csv")
wdddata <- read.csv("E:\\zuomian\\sswdd.csv")
idwdata <- read.csv("E:\\zuomian\\ssidw.csv")
summary(lm(idwdata$pre~idwdata$pm2.5))
p1 <- ggplot(data = rfdata, mapping = aes(x = rfdata[,3], y = rfdata[,2]),ylim=c(0,300)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep=""), size = 3),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep=""),, size = 3))+
  annotate("text", x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.877)", size = 4) +
  annotate("text", x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==9.024)", size = 4)	+ theme_bw() +
  annotate("text", x = 150, y = 15, parse = TRUE, 
           label = "Y == -2.139+1.052*X", size = 4) +
  geom_abline(intercept = -2.139, slope = 1.052,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  labs(tag = "(a)") +
  scale_fill_viridis() +
  ggtitle("RF") +
  scale_x_continuous(limits = c(0, 200))+
  #theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()
p1+theme(plot.title = element_text(hjust = 0.5))  
p2 <- ggplot(data = w00data, mapping = aes(x = w00data[,3], y = w00data[,2]),ylim=c(0,300)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep="")),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep="")))+
  annotate("text", x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.910)", size = 4) +
  annotate("text", x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==7.726)", size = 4)	+ theme_bw() +
  annotate("text", x = 150, y = 15, parse = TRUE, 
           label = "Y == -1.473+1.012*X", size = 4) +
  geom_abline(intercept = -0.595, slope = 1.036,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  scale_fill_viridis() +
  labs(tag = "(c)") +
  ggtitle(expression(paste("w"[0][0],"-RF",sep=""))) +
  theme_classic()
p2+theme(plot.title = element_text(hjust = 0.5)) 
p3 <- ggplot(data = w0ddata, mapping = aes(x = w0ddata[,3], y = w0ddata[,2]),ylim=c(0,300)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep="")),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep="")))+
  annotate("text", x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.910)", size = 4) +
  annotate("text", x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==7.706)", size = 4)	+ theme_bw() +
  annotate("text", x = 150, y = 15, parse = TRUE, 
           label = "Y == -1.404+1.034*X", size = 4) +
  geom_abline(intercept = -1.404, slope = 1.034,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  scale_fill_viridis() +
  labs(tag = "(d)") +
  ggtitle(expression(paste("w"[0][d],"-RF",sep=""))) +
  theme_classic()
p3+theme(plot.title = element_text(hjust = 0.5)) 
p4 <- ggplot(data = wd0data, mapping = aes(x = wd0data[,3], y = wd0data[,2]),ylim=c(0,300)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep="")),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep="")))+
  annotate("text",x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.917)", size = 4) +
  annotate("text", x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==7.420)", size = 4)	+ theme_bw() +
  annotate("text", x = 150, y = 15, parse = TRUE, 
           label = "Y == -0.653+1.014*X", size = 4) +
  geom_abline(intercept = -0.653, slope = 1.014,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  scale_fill_viridis() +
  labs(tag = "(e)") +
  ggtitle(expression(paste("w"[d][0],"-RF",sep=""))) +
  theme_classic()
p4
p5 <- ggplot(data = wdddata, mapping = aes(x = wdddata[,3], y = wdddata[,2]),ylim=c(0,300)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep="")),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep="")))+
  annotate("text", x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.918)", size = 4) +
  annotate("text",x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==7.372)", size = 4)	+ theme_bw() +
  annotate("text",x = 150, y = 15, parse = TRUE, 
           label = "Y == -0.595+1.012*X", size = 4) +
  geom_abline(intercept = -0.595, slope = 1.012,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  scale_fill_viridis() +
  labs(tag = "(f)") +
  ggtitle(expression(paste("w"[d][d],"-RF",sep=""))) +
  theme_classic()
p5
p6 <- ggplot(data = idwdata, mapping = aes(x = idwdata[,3], y = idwdata[,2]),ylim=c(0,300),xlim=c(0,150)) +
  geom_bin2d(bins = 200) + # bins控制着图中每个箱子的大小
  labs(x=expression(paste("PM"[2.5]," ","Predictions(ug/","m"^3,")",sep="")),
       y=expression(paste("PM"[2.5]," ","Measurements(ug/","m"^3,")",sep="")))+
  annotate("text", x = 150, y = 35, parse = TRUE, 
           label = "atop(R^2==0.917)", size = 4) +
  annotate("text", x = 150, y = 20, parse = TRUE, 
           label = "atop(RMSE==7.400)", size = 4)	+ theme_bw() +
  annotate("text", x = 150, y = 15, parse = TRUE, 
           label = "Y == -0.554+1.011*X", size = 4) +
  geom_abline(intercept = -0.554, slope = 1.011,
              colour = "73", size = 0.5) +
  geom_abline(intercept = 0, slope = 1,lty=2,
              colour = "red", size = 0.5)+
  scale_fill_viridis() +
  expand_limits(x = 0, y = 0)+
  labs(tag = "(b)") +
  ggtitle("idw-RF") +
  theme_classic()
p6
library(cowplot)
plot_grid(p1+theme(plot.title = element_text(hjust = 0.5)) , p6+theme(plot.title = element_text(hjust = 0.5)) , 
          p2+theme(plot.title = element_text(hjust = 0.5)) , p3+theme(plot.title = element_text(hjust = 0.5)) ,
          p4+theme(plot.title = element_text(hjust = 0.5)) ,p5+theme(plot.title = element_text(hjust = 0.5)) ,nrow = 2)
ggsave("cv.tiff", units="in", width=13.5, height=7.5, dpi=200, compression = 'lzw')






