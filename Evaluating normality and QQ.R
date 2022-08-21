
#Study normality by plotting-----------------------
library(readxl)
data<-read_xlsx("A4_Ex01.xlsx", col_names =FALSE)
colnames(data)<-paste0("x",1:6)
data$x6<-factor(data$x6)
table(data$x6)
str(data)
s<-seq(min(data$x1), max(data$x1), length.out = nrow(data)) #building min and max for variable X1
dx1<-dnorm(s, mean=mean(data$x1), sd=sd(data$x1))
hist(data$x1, probability = TRUE) #probability =TRUE to put on top after density function

lines(s,dx1, col="red") # (Kernel desity plot), plotting distribution of my variables assuming normal distribution
abline(v=mean(data$x1), col=3, lwd=2)
lines(density(data$x1), col="blue") #ploting density of my variables - Kernel plot?


#Plot QQp for variable x1 from data---------
install.packages("car")
library("car")
qqp(data$x1, distribution="norm", mean=mean(data$x1), sd=sd(data$x1)) #graphic library for distributions


#Relate observations relate to normal distibution--------
data2<-c(-1,-0.1, 0.16,0.41,0.62,0.8,1.26,1.54,1.71,2.3)
length(data2) #checking that I typed all 10 values in data2
pi<-(rank(data2)-0.5)/10;pi  #ranking my real data
qi<-qnorm(pi,mean=0,sd=1);qi #quantiles for my ranked data

plot(qi,data2, xlab="quantis teorical",ylab="quantis empiricos", #emprical means my own data before ranking
     ylim=c(-1,2.5),xlim = c(-1,2.5)) #teoretical quantiles are ranked and estimated, 
abline(lm(data2~qi),col=2) #addind lm between my data and teoretical

#Calculate c2 and relate to qchisq distribution----------

x1<-c(126974, 96933, 86656, 63438,55264,50976, 39069,36156,35209,32416);x1
x2<-c(4224, 3835, 3510, 3758,3939,1809, 2946,359,2480,2413);x2
length(x1)
length(x2)
x<-matrix(c(x1,x2),10,2);x #make matrix with vectors
xbar<-colMeans(x);xbar
S<-var(x);S


#C2 --------------------
c2<-numeric()
for(i in 1:nrow(x)){
  
  c2[i]<-t(matrix(x[i,]-xbar,2,1))%*%solve(var(x))%*%matrix(x[i,]-xbar,2,1) #this is function to do for all matrix
  
  
}
c2  # C1_x value is the same as for c1_x estimated above, function works
sortedc2<-sort(c2);sortedc2

#Plot result 
car::qqPlot(sortedc2, distribution="chisq", df=2)

#Plot other way c2 (xi) against qchisq (teoretical quantiles fo distribution Xp2)

p<-nrow(S);p
n<-length(sortedc2);n
pi<-(rank(sortedc2)-0.5)/n;pi
qi_qcisq<-qchisq(pi,df=2);qi_qcisq #df=p dont forget

plot(sortedc2,qi_qcisq, xlab="quantis teorical",ylab="quantis empiricos", #emprical means my own data before ranking
     ylim=c(0,6),xlim = c(0,6)) #teoretical quantiles are ranked and estimated, 
abline(lm(sortedc2~qi_qcisq),col=2) #adding lm between my data and teoretical

#Plot result by calculating pi and qi - correct!

c2example<-c(1.265,2.044,2.649,2.865,3.177);c2example #this is variables xi
pi<-(rank(c2example)-0.5)/5;pi  #empirical probabilities
qi_qcisq<-qchisq(pi,df=3);qi_qcisq #df=p based on mahanobis or qchisq distribution 

plot(c2example,qi_qcisq, xlab="quantis teorical",ylab="quantis empiricos", #emprical means my own data before ranking
     ylim=c(0.2,7),xlim = c(0.2,7)) #teoretical quantiles are ranked and estimated, 
abline(lm(c2example~qi_qcisq),col=2) #adding lm between my data and teoretical


#Calculate u and relate to beta distribution--------------------------------
n<-nrow(x)
p<-ncol(x);p
alpha <-p/2;alpha
beta<-(n-p-1)/2;beta
u<-(n*c2)/(n-1)^2;u

pi<-(rank(u)-0.5)/n;pi
qi_beta<-qbeta(pi, shape1=p/2, shape2=(n-p-1)/2);qi_beta #this is based on u distance calculations

#Plot results

plot(u,qi_beta, xlab="quantis teorical",ylab="quantis empiricos", #emprical means my own data before ranking
     ylim=c(0,0.6),xlim = c(0,0.6)) #teoretical quantiles are ranked and estimated, 
abline(lm(u~qi_beta),col=2)

car::qqp(u,distribution="beta",shape1=p/2,shape2=(n-p-1)/2) 





