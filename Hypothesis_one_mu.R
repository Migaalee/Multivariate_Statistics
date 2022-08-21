#Hypothesis T2obs with one mean; mu=60 or not equal------------------
x<-c(62,62,68,48,51,60,51,57,57,41,62,50,53,34,62,61)
m0<-60
mean<-mean(x);mean
sd<-sd(x);sd
n<-length(x);n
# Calculate t2obs and p-value----------------------------
t2obs<-(mean-m0)/(sd/sqrt(length(x)));t2obs #this is hypothesis test
pcum<-pt(abs(t2obs),df=n-1);pcum
pvalue<-2*(1-pcum);pvalue #with 2 tails

#Test hypothesis with given alpha----
alpha<-0.01
n=16 #or nrow in matrix
p=1  #or ncol in matrix

T<-qt(p=(1-alpha)/2, df=n-1);T #this is equal to -0.070799, not sure if qt or pt

#If t2obs is more that T, we reject H0


#Test hypothesis T2obs with vector; mu0=mu=[11,3]--------------------
x5<-matrix(c(6,10,8,9,6,3),3,2);x5
mu<-matrix(c(11,3),2,1);mu
x5med<-colMeans(x5);x5med
S5<-var(x5);S5
Ss<-solve(S5);Ss
n<-nrow(x5);n
p<-ncol(x5);p
alpha<-0.05

#T2obs of vector mu0=mu=[11,3]---------------
T2obs<-n*t((x5med-mu))%*%solve(S5)%*%(x5med-mu);T2obs

#Significance and if we accept or reject----------
F<-qf(p=1-alpha, df1=p,df2=n-p);F #199.5
coef<-(p*(n-1)/(n-p));coef
tt<-(p*(n-1)/(n-p))*F;tt  #798
#we reject if T2obs is more that tt
#T2=7 < 4F(2,1);0.95= 798, we do not reject
#Rejection region ---------
#RR[4F(2,1);0.95;+infinity[798;+infinity]]

#Confidence region of mu-----

alpha<-0.05
F6<-qf(p=1-alpha, df1=p,df2=n-p);F6
conf_r<-(p*(n-1))/(n*(n-p))*F6;conf_r #confidence interval/region for mu!
#Confidence region T2obs<=tt6, so 7<=266

#Ellipsoide confidence region for mu with centre x(mean) and axis defined ------------------

F<-qf(p=1-alpha, df1=p,df2=n-p);F
c<-sqrt((p*(n-1)/n*(n-p))*F);c
c<-sqrt((p*(n-1))/(n*(n-p))*F);c
c<-sqrt(conf_r);c

alpha = 0.05


#Plot axis of ellipsoide--------------------

center<-x5med;center

length<-eigen(S5)$val;length

direction<-eigen(S5)$vec;direction

axis1<-c*sqrt(length[1])*matrix(c(direction[,1]),2,1);axis1

#vectors of principal axis-----------------
axis11<-x5med-axis1;axis11
axis12<-x5med+axis1;axis12

axis2<-c*sqrt(length[2])*matrix(c(direction[,2]),2,1);axis2

#vectors of secondary axis----------------
axis21<-x5med-axis2;axis21
axis22<-x5med+axis2;axis22

#Make a plot
#plot
# merge into one matrix
elipse_vectors = cbind(axis11,axis12,axis21,axis22);elipse_vectors
# draw
plot(elipse_vectors[1,],elipse_vectors[2,],col='blue')
points(center[1],center[2],col='green')

plot(mu[1],mu[2],ylim=c(2,12),xlim=c(2,12))
# axis1_1
arrows(center[1],center[2], axis11[1], axis11[2])
#axis1_2
arrows(center[1],center[2], axis12[1], axis12[2])
# axis2_1
arrows(center[1],center[2], axis21[1], axis21[2])
#axis2_2
arrows(center[1],center[2], axis22[1], axis22[2])


#if c is given find length and direction of axis?----



#Full Hypothesis test from matrix X----------
library(readxl)
data6<-read_xlsx("A4_Ex06.xlsx", col_names =TRUE)
x1<-data6$v1^0.25
x2<-data6$v2^0.25
x6<-matrix(c(x1,x2), nrow = nrow(data6),ncol=2) #matrix any
x6med<-colMeans(x6);x6med
mu6<-matrix(c(0.562,0.589),2,1);mu6
S6<-var(x6);S6
n<-nrow(x6);n
p<-ncol(x6);p

#Test statistic T2obs----
Ss<-solve(S6);Ss
Ts6<-n*t((x6med-mu6))%*%solve(S6)%*%(x6med-mu6);Ts6 #T2obs is 1.2573
alpha<-0.05
F6<-qf(p=1-alpha, df1=p,df2=n-p);F6 #p is 1-alpha, df1=p and df2=n-p=3-2=1
tt<-(p*(n-1)/(n-p))*F6;tt #6.62 #Ts2(41)=(2*41/41-2+1)*F(2;40);0.95

#Confidence region for mu-------
alpha<-0.05
F6<-qf(p=1-alpha, df1=p,df2=n-p);F6 #3.23
tt6<-(p*(n-1))/(n*(n-p))*F6;tt6  #confidence interval/region for mu!


#Data for confidence intervals
n = nrow(x6);n
p = ncol(x6);p
xbar = 1/n * t(x6) %*% rep(1,n) ; 
xbar=colMeans(x6);xbar
S = cov(x6);S[2,2]

alpha = 0.05
coef = p*(n-1)/(n*(n-p))
df1 = p
df2 = n-p

# intervals simultaneos--------------
interval_matrix_simultaneo = matrix(0,ncol=2,nrow=p)
for (i in 1:p){
  interval_up = xbar[i] + sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(S[i,i]/n)
  interval_down = xbar[i] - sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(S[i,i]/n)
  interval_matrix_simultaneo[i,] = c(interval_down,interval_up)
}
interval_matrix_simultaneo

pIC1P <- xbar[1] + sqrt(coef*qf(1-alpha,df1,df2) * (S[1,1]/n));pIC1P 
pIC1M <- xbar[1] - sqrt(coef*qf(1-alpha,df1,df2) * (S[1,1]/n));pIC1M
pIC2P <- xbar[2] + sqrt(coef*qf(1-alpha,df1,df2) * (S[2,2]/n));pIC2P 
pIC1M <- xbar[2] - sqrt(coef*qf(1-alpha,df1,df2) * (S[2,2]/n));pIC1M


# one at a time----------------

interval_matrix_one_at_time = matrix(0,ncol=2,nrow=p)
df = n-1
for (i in 1:p){
  interval_up = xbar[i] + qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_down = xbar[i] - qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_matrix_one_at_time[i,] = c(interval_down,interval_up)
}

interval_matrix_one_at_time

# Bonferroni------------------
alpha = alpha/p

interval_matrix_bonf = matrix(0,ncol=2,nrow=p)
df = n-1
for (i in 1:p){
  interval_up = xbar[i] + qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_down = xbar[i] - qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_matrix_bonf[i,] = c(interval_down,interval_up)
}

interval_matrix_bonf

interval_matrix_simultaneo[2,]
interval_matrix_one_at_time[2,]
interval_matrix_bonf[2,]











