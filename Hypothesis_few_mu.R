#1) PAIRED POPULATIONS ------------

#Data given
A<-matrix(c(6,6,18,8,11,34,28,71,43,33,20,27,23,64,44,30,75,26,124,54,30,14),11,2) ;A #with x1 and x2
B<-matrix(c(25,28,36,35,15,44,42,54,34,29,39,15,13,22,29,31,64,30,64,56,20,21),11,2) ;B #with x1 and x2


#T2obs test-----------
d1<-A[,1]-B[,1];d1
d2<-A[,2]-B[,2];d2
d <- matrix(c(d1,d2),11,2);d
d_mean <- colMeans(d);d_mean
S <- cov(d);S
p=2
n=11
Tobs <- n %*% t(d_mean) %*% solve(S) %*% d_mean; Tobs #13.6391

#T2 distribution-----------

#T2~T2p(n-1)=(p(n-1)/n-p)*F(p;n-p)


#T2 rejection test-----------
#Tobs>=(p(n-1)/n-p)*F(p,n-p);1-alpha
n=11
p=2
alpha=0.05
F<-qf(p=1-alpha, df1=p,df2=n-p);F #4.2565
tt<-(p*(n-1)/(n-p))*F;tt  #9.458877
#we reject if T2obs is more that tt
#Because T2obs(13.63) is more than tt (9.45), we reject hypothesis
(p*(n-1)/(n-p))
#p(n-1)/n-p=2.22
#Rejection region [2.22 F(2,9);0.95;+infinity[9.45;+infinity]]

#T2 confidence region for mu-----------
alpha<-0.05
F<-qf(p=1-alpha, df1=p,df2=n-p);F  #4.25
conf_r<-(p*(n-1))/(n*(n-p))*F;conf_r #conf_region 0.8598

#T2obs(13.6391)<=0.8598

#T2 elipsoide and axis-----------

#Ellipsoide confidence region for mu with centre x(mean)
F<-qf(p=1-alpha, df1=p,df2=n-p);F #4.25
c<-sqrt((p*(n-1))/(n*(n-p))*F);c #0.92
c<-sqrt(conf_r);c  #0.92
alpha = 0.05

center<-d_mean;center

length<-eigen(S)$val;length

direction<-eigen(S)$vec;direction

axis1<-c*sqrt(length[1])*matrix(c(direction[,1]),2,1);axis1

#vectors of principal axis-----------------
axis11<-center-axis1;axis11
axis12<-center+axis1;axis12

axis2<-c*sqrt(length[2])*matrix(c(direction[,2]),2,1);axis2

#vectors of secondary axis----------------
axis21<-center-axis2;axis21
axis22<-center+axis2;axis22


#T2 Simultaneous confidence intervals-----------
#data
n = nrow(A);n
p = ncol(A);p
xbar=d_mean;xbar
S = cov(d);S
alpha = 0.05
coef = (p*(n-1)/(n-p));coef
df1 = p
df2 = n-p


interval_matrix_simultaneo = matrix(0,ncol=2,nrow=p)
for (i in 1:p){
  interval_up = xbar[i] + sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(S[i,i]/n)
  interval_down = xbar[i] - sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(S[i,i]/n)
  interval_matrix_simultaneo[i,] = c(interval_down,interval_up)
}

interval_matrix_simultaneo

pIC1P <- xbar[1] + sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(S[1,1]/n);pIC1P 
pIC1M <- xbar[1] - sqrt(tt * (S[1,1]/n));pIC1M
pIC2P <- xbar[2] + sqrt(tt * (S[2,2]/n));pIC2P 
pIC1M <- xbar[2] - sqrt(tt * (S[2,2]/n));pIC1M


#T2 Bonferroni intervals-----------

alpha = alpha/p;alpha

interval_matrix_bonf = matrix(0,ncol=2,nrow=p)
df = n-1;df
for (i in 1:p){
  interval_up = xbar[i] + qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_down = xbar[i] - qt(1-alpha/2, df) * sqrt(S[i,i]/n)
  interval_matrix_bonf[i,] = c(interval_down,interval_up)
}

interval_matrix_bonf


#2) INDEPENDENT POPULATIONS ------------

#Data given
library(readxl)
data<- readxl::read_excel("A5_Ex02.xlsx")
data1<-data[data$grupo==1,-1];data1 #has x1 and x2 columns
data2<-data[data$grupo==2,-1];data2 #has x1 and x2 
mu0<-c(170,250)

#Calculate statistics
n1=nrow(data1);n1
p1=ncol(data1);p1
xbar1<-colMeans(data1);xbar1
S1<-cov(data1);S1
n2=nrow(data2);n2
p2=ncol(data2);p2
xbar2<-colMeans(data2);xbar2
S2<-cov(data2);S2


#T2obs test if first group has mu(170,250)-----------
alpha=0.01
t2obs<-n1*t(xbar1-mu0)%*%solve(S1)%*%(xbar1-mu0);t2obs  #2.0406
F<-qf(p=1-alpha, df1=p1,df2=n1-p1);F #5.098579
tt<-(p1*(n1-1)/(n1-p1))*F;tt  #10.41884

#we reject if T2obs is more that tt
#T2obs (2.04)<tt, we do not reject, we accept H0 hypothesis

#T2obs test if second group has mu(170,250)-----------

alpha=0.01
t2obs2<-n2*t(xbar2-mu0)%*%solve(S2)%*%(xbar2-mu0);t2obs2  #26.26257
F<-qf(p=1-alpha, df1=p2,df2=n2-p2);F #5.05661
tt2<-(p2*(n2-1)/(n2-p2))*F;tt2  #10.315
#we reject if T2obs is more that tt
#T2obs (26.26)>tt2(10.315), we reject H0 hypothesis

#T2obs test for mean is same between two groups-Spool-----------

n1=145
n2=244
p=2
S1<-matrix(c(0.8593,0.9178,0.9178,1.7951),2,2);S1
S2<-matrix(c(0.9996,1.0110,1.0110,1.9927),2,2);S2
xbar1<-c(5.9827,16.0784);xbar1
xbar2<-c(5.1167,15.1926);xbar2

Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled
t2obs_S<-t(xbar1-xbar2)%*%solve((1/n1+1/n2)*Spooled)%*%(xbar1-xbar2);t2obs_S  #9.24
tt_S<-((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);tt_S #9.76 confidence region
ifelse(t2obs_S>=tt_S, "reject H0", "Do not reject H0")

#T2 distribution between two mean-----------

#T2~T2p(n1+n2-2))=(((n1+n2-2)*p)/(n1+n2-p-1))*F(p;n1+n2-p-1)

#T2 rejection region between two mean-----------
#[98*2/97 F(2,97);0.95;+infinity[9.76;+infinity]]

#T2obs test for mean difference between two groups-Spool-----------
n1=145
n2=244
p=2
S1<-matrix(c(0.8593,0.9178,0.9178,1.7951),2,2);S1
S2<-matrix(c(0.9996,1.0110,1.0110,1.9927),2,2);S2
xbar1<-c(5.9827,16.0784);xbar1
xbar2<-c(5.1167,15.1926);xbar2
mu1minusmu2<-c(0,1);mu1minusmu2

Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled
t2obs_S<-t((xbar1-xbar2)-(mu1minusmu2))%*%solve((1/n1+1/n2)*Spooled)%*%((xbar1-xbar2)-(mu1minusmu2));t2obs_S  
tt_S<-((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);tt_S 
ifelse(t2obs_S>=tt_S, "reject H0", "Do not reject H0")

# ???T2Obs for a=(0,1) amu1=amu2 for mean the same between groups----
n1=145
n2=244
p=2
S1<-matrix(c(0.8593,0.9178,0.9178,1.7951),2,2);S1
S2<-matrix(c(0.9996,1.0110,1.0110,1.9927),2,2);S2
xbar1<-c(5.9827,16.0784);xbar1
xbar2<-c(5.1167,15.1926);xbar2
a<-c(1,0)

Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled
t2obs_S<-t(xbar1-xbar2)%*%solve((1/n1+1/n2)*Spooled)%*%(xbar1-xbar2);t2obs_S  #9.24
tt_S<-((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);tt_S #9.76
ifelse(t2obs_S>=tt_S, "reject H0", "Do not reject H0")


#T2 confidence region between mu1 and mu2-----------
n1=145
n2=244
p=2
alpha<-0.05
c2=(1/n1+1/n2)*((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);c2 #0.06656
#T2obs(calculated)<=0.06656

#T2 elipsoide and axis-----------
p<-2
n1<-50
n2<-50
m1<-matrix(c(8.3, 4.1),2,1) ; m1
m2<-matrix(c(10.2, 3.9),2,1) ; m2
S1<-matrix(c(2, 1, 1, 6),2,2) ; S1
S2<-matrix(c(2,1,1,4),2,2) ; S2
alpha=0.05
num<-(n-1)*S1+(n2-1)*S2
den<-n1+n2-2
Spool<-num/den

c2=(1/n1+1/n2)*((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);c2
center=m1-m2;center #Ellipse is centered at -1.9 and 0.2

eval <- eigen(Spooled)$val;eval #length of semi-axis

first_axis<-sqrt(c2*eval[1]);first_axis #1.1508 length of axis
second_axis<-sqrt(c2*eval[2]);second_axis #0.65107 length
evec <- eigen(Spooled)$vec;evec #direction of ellipse

#First axis
center-sqrt(c2*eval[1])*evec[1]
center+sqrt(c2*eval[1])*evec[1]
#Second axis
center-sqrt(c2*eval[2])*evec[2]
center+sqrt(c2*eval[2])*evec[2]


#T2 Simultaneous confidence intervals-----------

n1=145
n2=244
p=2
S1<-matrix(c(0.8593,0.9178,0.9178,1.7951),2,2);S1
S2<-matrix(c(0.9996,1.0110,1.0110,1.9927),2,2);S2
xbar1<-c(5.9827,16.0784);xbar1
xbar2<-c(5.1167,15.1926);xbar2
xbar<-xbar1-xbar2;xbar
Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled


alpha = 0.04 #DO NOT FORGET TO CHANGE ALPHA


df1 = p;p
df2 = n1+n2-p-1;df2

coef = ((n1+n2-2)*p)/(n1+n2-p-1);coef

interval_matrix_simultaneo = matrix(0,ncol=2,nrow=p)
for (i in 1:p){
  interval_up = xbar[i] + sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[i,i])
  interval_down = xbar[i] - sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[i,i])
  interval_matrix_simultaneo[i,] = c(interval_down,interval_up)
}

interval_matrix_simultaneo



pIC1P <- xbar[1]+ sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[1,1]);pIC1P 
pIC1M <- xbar[1] -sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[1,1]);pIC1M
pIC2P <- xbar[2] + sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[2,2]);pIC2P 
pIC1M <- xbar[2] - sqrt(coef*qf(1-alpha,df1,df2)) * sqrt(((1/n1)+(1/n2))*Spooled[2,2]);pIC1M






#T2 Bonferroni intervals-----------
n1=145
n2=244
p=2
S1<-matrix(c(0.8593,0.9178,0.9178,1.7951),2,2);S1
S2<-matrix(c(0.9996,1.0110,1.0110,1.9927),2,2);S2
xbar1<-c(5.9827,16.0784);xbar1
xbar2<-c(5.1167,15.1926);xbar2
xbar<-xbar1-xbar2;xbar
Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled


alpha = 0.04 #DO NOT FORGET TO CHANGE ALPHA
alpha = alpha/p;alpha

interval_matrix_bonf = matrix(0,ncol=2,nrow=p)
df = n1+n2-2;df
for (i in 1:p){
  interval_up = xbar1[i]-xbar2[i] + qt(1-alpha/2, df) *  sqrt(((1/n1)+(1/n2))*Spooled[i,i])
  interval_down = xbar1[i]-xbar2[i] - qt(1-alpha/2, df) * sqrt(((1/n1)+(1/n2))*Spooled[i,i])
  interval_matrix_bonf[i,] = c(interval_down,interval_up)
}

interval_matrix_bonf

pIC1P <- xbar[1]+  qt(1-alpha/2, df) * sqrt(((1/n1)+(1/n2))*Spooled[1,1]);pIC1P 
pIC1M <- xbar[1] - qt(1-alpha/2, df) * sqrt(((1/n1)+(1/n2))*Spooled[1,1]);pIC1M
pIC2P <- xbar[2] + qt(1-alpha/2, df) * sqrt(((1/n1)+(1/n2))*Spooled[2,2]);pIC2P 
pIC1M <- xbar[2] -  qt(1-alpha/2, df)* sqrt(((1/n1)+(1/n2))*Spooled[2,2]);pIC1M









