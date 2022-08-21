
#Continous distances
  
x=c(7,4,2,0,9,6)
y=c(3,5,4,1,7,8)
E=cbind(x,y);E

##Minkowski distance----

d_minkowski<-dist(E,method = "minkowski");round(d_minkowski,4)

##Manhattan or city block------

d_manhattan<-dist(E,method = "manhattan");round(d_manhattan,4)

##Euclidean distance----

d_euclidean<-dist(E,method = "euclidean");round(d_euclidean,4)

##Mahalanobis distance -------

Se<-cov(E);Se
meane<-colMeans(E);meane
d_mahalanobis<-mahalanobis(E,meane,Se);d_mahalanobis
#2.8624661 0.1727642 0.7994580 2.3678862 1.6971545 2.1002710

#or

X_cov <- var(scale(E)) # standardize first
X_cov <- var(E)
X_mean <- apply(E,2,mean)
X_mah <- mahalanobis(E, X_mean, X_cov);X_mah

#2.8624661 0.1727642 0.7994580 2.3678862 1.6971545 2.1002710


##Chebyshev(maximum) distance----------


d_maximum<-dist(E,method = "maximum");round(d_maximum,4)



##Canberra distance---------


d_canberra<-dist(E,method = "canberra");round(d_canberra,4)



#Binary variables

x1=c(0,1,0,0,1)
x2=c(0,1,1,0,1)
x3=c(0,1,0,1,1)
x4=c(1,0,1,0,0)
x5=c(1,1,1,1,0)
x6=c(1,0,0,1,0)
X=cbind(x1,x2,x3,x4,x5,x6);X


#Russel index-------

russel<-function(obj1,obj2) {
  tab=as.numeric((table(obj1,obj2)))
  d=tab[1]
  c=tab[2]
  b=tab[3]
  a=tab[4]
  russel=a/(a+b+c+d)
  return(russel)
  
}

russel(X[,1],X[,2])




#Jaccard Index-------------

jac<-function(obj1,obj2) {
  tab=as.numeric((table(obj1,obj2)))
  d=tab[1]
  c=tab[2]
  b=tab[3]
  a=tab[4]
  jac=a/(a+b+c)
  return(jac)
  
}

jac(X[,1],X[,2])


#Other function

library(clusteval)
jaccard <- function(df, margin) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}

df2 <- data.frame(
  IDS = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
  CESD = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

jaccard(df2, 1) # JSim     JDist 
                 #0.1860465 0.8139535 



df3<-data.frame(X[,1],X[,2])

jaccard(df3, 1) #gives the same first value as jac function
# JSim     JDist 
#0.6666667 0.3333333


#Sorenson index--------------

sorenson<-function(obj1,obj2) {
  tab=as.numeric((table(obj1,obj2)))
  d=tab[1]
  c=tab[2]
  b=tab[3]
  a=tab[4]
  sorenson=2*a/(2*a+b+c)
  return(sorenson)
  
}

sorenson(X[,1],X[,2])


#Sokal and Sneath index------------


sokal<-function(obj1,obj2) {
  tab=as.numeric((table(obj1,obj2)))
  d=tab[1]
  c=tab[2]
  b=tab[3]
  a=tab[4]
  sokal=2*(a+d)/(2*(a+d)+b+c)
  return(sokal)
  
}

sokal(X[,1],X[,2])


#Similarity between units using any binary method----

sim=matrix(NA,5,5)  #similarity between 5 unities 
for (i in 1:5) {
  for (j in 1:5) {
    sim[i,j]=jac(X[i,],X[j,])
  }
}

sim


1-dist(X,method="bin")






# Methods to aggregate clusters----
  
d1<-matrix(c(0,9,3,6,11,9,0,7,5,10,3,7,0,9,2,6,5,9,0,8,11,10,2,8,0), 5, 5);d1

library(readxl)
data=read_xlsx("A11_Ex03.xlsx", col_names = TRUE)
head(data)
d=dist(data[,-1],method="euclidean") #euclidean is default method
d

  
#nearest neighbor ou single-linkage method-------
#minimum distance
fit1=hclust(d,method = "single")
plot(fit1,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit1,k=4,border="red")


#furthest neighbor ou complete-linkage method---------
fit2=hclust(d,method = "complete")
plot(fit2,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit2,k=3,border="red")


#average-linkage method-------------------


fit3=hclust(d,method = "average")
plot(fit3,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit3,k=3,border="red")


#centroid clustering-----------
?hclust
fit5=hclust(d,method = "centroid")
plot(fit5,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit5,k=3,border="red")

#median clustering---------
fit6=hclust(d,method = "median")
plot(fit6,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit6,k=3,border="red")


#Ward's method-----------
fit4=hclust(d,method = "ward.D2")
plot(fit4,hang=-1,ylab="Distance euclidean",labels = data$City)
rect.hclust(fit4,k=3,border="red")





#kmeans clustering-----------
library(readxl)
dados<-read_xlsx("A12_Ex04.xlsx", col_names = TRUE)
dados
variance<-diag(cov(dados[,-1]));variance #variance checking if to use cov or corr
z=scale(dados[,-1], center=TRUE, scale=TRUE) #scale data since variance is huge


#Kmeans with 5 centers-------

set.seed(1)
fitkm=kmeans(z,centers = 5, nstart = 5)
str(fitkm)
fitkm


#Centroids--------
fitkm$centers
aggregate(z,list(fitkm$cluster), mean) #group based on which cluster it belongs to


#variation-------------
fitkm$withinss #clusters 2, 3 and 5 are more heterogeneos

#5.900318 16.994661 22.110431  8.012133 18.925874



# method WWSj within sum of squares------------


temp=cbind(z, fitkm$cluster)
temp
cent=fitkm$centers
zclus1=temp[temp[,10]==1,-10]
D1=zclus1-c(rep(1,4))%*%t(cent[1,])
sum(diag(t(D1)%*%D1)) #variation for this cluster #I got 8.012 with nstart=20 and teacher 5.900318 nstart=5


#Visualisation-----------
library(cluster)
library(fastcluster)
library(ggplot2)
library(factoextra)
library(gridExtra)

#1st method
par(pty="s")
clusplot(z,fitkm$cluster, lines = 0, color=TRUE, shade=TRUE, labels=4) #labels 2 shows



#2nd method
set.seed(1)
km2=kmeans(z,centers=2,nstart = 20)#kmeans will run 5 times simulation and choose best with less variability
km3=kmeans(z,centers=3,nstart = 20)
km4=kmeans(z,centers=4,nstart = 20)
km5=kmeans(z,centers=5,nstart = 20)



p1=fviz_cluster(km2, data=z) +ggtitle("k=2") #show data
p1


#How many clusters to choose? Elbow method----------

#Numero otimo de clusters - elbow method 

fviz_nbclust(z,kmeans,method = "wss",nstart=20) #elbow would tell minimum number of clusters for analysis


#Silhouette method
fviz_nbclust(z,kmeans,method = "sil",nstart=20) 


#how well each cluster is defined? Silhouette method-----------

plot(silhouette(km3$cluster, dist(z))) #the closer to 1, better defined cluster





























