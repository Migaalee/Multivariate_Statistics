#Pvalue for mu=60, one population, one mu-------
t2obs<-(mean-m0)/(sd/sqrt(length(x)));t2obs #this is hypothesis test
pcum<-pt(abs(t2obs),df=n-1);pcum
pvalue<-2*(1-pcum);pvalue #with 2 tails

#Pvalue for mu=mu0=mu=[11,3], one population, few observations x1, x2(need to check)-----
#T2obs<-n*t((x5med-mu))%*%solve(S5)%*%(x5med-mu);T2obs
T2obs<-7;T2obs
n=3;n
p=2;p

df1=p;p
df2=n-p;df2
virsus<-n-p;virsus
apacia<-p*(n-p);apacia
dalyba<-virsus/apacia;dalyba
pvalue_M<-1-pf(q=T2obs*dalyba,df1 = df1,df2=df2);pvalue_M


#Pvalue should for two mean test-------

T2OBS=9.2484;T2OBS
n1=48
n2=52
p=2
df1=p;p
df2=n1+n2-p-1;df2
virsus<-(n1+n2-p-1);virsus
apacia<-(n1+n2-2)*p;apacia
dalyba<-virsus/apacia;dalyba

pvalue_M<-1-pf(q=T2OBS*dalyba,df1 = df1,df2=df2);pvalue_M

#reject H0 if pvalue<=alpha, but our pvalue(0.12604)>alpha(0.01), we do not reject H0

#Another example
n1=32
n2=32
p=4
S1<-matrix(c(5.192,4.545,6.522,5.250,4.545,13.18,6.760,6.266,6.522,6.760,28.67,14.47,5.250,6.266,14.47,16.65),4,4);S1
S2<-matrix(c(9.136,7.549,4.864,4.151,7.549,18.60,10.22,5.446,4.864,10.22,30.04,13.49,4.151,5.446,13.49,28),4,4);S2
xbar1<-c(15.97,15.91,27.19,22.75);xbar1
xbar2<-c(12.34,13.91,16.66,21.94);xbar2
Spooled <- (n1-1)*S1/(n1+n2-2) + (n2-1)*S2/(n1+n2-2); Spooled #ok

alpha<-0.05
F<-qf(1-alpha,p,n1+n2-p-1);F  #4.25
conf_r<-(p*(n-1))/(n*(n-p))*F;conf_r #conf_region 0.8598
t2obs_S<-t((xbar1-xbar2))%*%solve((1/n1+1/n2)*Spooled)%*%((xbar1-xbar2));t2obs_S  #ok
tt_S<-((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-alpha,p,n1+n2-p-1);tt_S 
ifelse(t2obs_S>=tt_S, "reject H0", "Do not reject H0")


df1=p;p
df2=n1+n2-p-1;df2
virsus<-(n1+n2-p-1);virsus
apacia<-(n1+n2-2)*p;apacia
dalyba<-virsus/apacia;dalyba

pvalue_M<-1-pf(q=t2obs_S*dalyba,df1 = df1,df2=df2);pvalue_M


#Pvalue for one mean test--------------

T2OBS=2.041;T2OBS
n=48
p=2
df1=p;p
df2=n-p;df2
virsus<-(n1+n2-p-1);virsus
apacia<-(n1+n2-2)*p;apacia
dalyba<-virsus/apacia;dalyba
pvalue_M<-1-pf(q=T2OBS*dalyba,df1 = df1,df2=df2);pvalue_M

#Because pvalue(0.374)>alpha(0.05), we do not reject H0



#Equation P(T(p=3,n=40))<2.77933332--------------

Test<-2.779332
n=41
p=3
df1=p;p
df2=n-p;df2
pvalue_M<-1-pf(q=Test,df1 = df1,df2=df2);pvalue_M  #0.054



#Probability P(ax>19.2291) with given a not sure if correct, close
mu<-c(3,-4,-2)
S<-matrix(c(17,18,19,18,22,21,19,21,22),3,3);S
a<-matrix(c(9,3,3),3,1);a
n=41
p=3
value=19.2291


mu_a<-t(a)%*%mu;mu_a
S_a<-t(a)%*%S%*%a;S_a
p_value<-1-pnorm(q=value,mean=mu_a,sd=S_a/n);p_value


#P(t(x-mu)*W-1*(x-mu))=c2) find P, when p=3-----------------------

pchisq(6.2514,3) #if distance c2 is given

qchisq(0.9,3) #if probability is given

#calculate probability more that 0.3518463 mahalanobis

pchisq(0.3518463,3) #if distance c2 is given 

1-pchisq(0.3518463,3)


#Ellipses bivariade pops (p=2), we have c2 mahalanobis?, get probability

pchisq(7.824046,2) #if distance c2 is given #mine was 98%

pchisq(7.01311,2) #if distance c2 is given #mine was 97% Pedro had

qchisq(0.9,3) #if probability is given



#???Ellipses 89% probability for mu, N7 distribution,variance unknown, find axis length

n=167
p=7
F6<-qf(p=0.89, df1=7,df2=160);F6 #1.71
conf_r<-(p*(n-1))/(n*(n-p))*F6;conf_r #for mu
#0.074385
c<-sqrt(conf_r);c #0.2727361 

#answer has to be 0.18???






