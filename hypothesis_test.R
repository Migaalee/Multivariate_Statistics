#Number Hypothesis testing - question 2

S2<-matrix(c(0.86,0.78,0.78,1.71),2,2);S2
mu<-c(5,10);mu
x2med<-c(6.02,16.1);x2med
n<-167
p<-2;p
alpha<-0.05

test<-solve(S2);test

#T2obs of vector mu0=mu=[5,10]---------------
T2obs<-n*t((x2med-mu))%*%solve(S2)%*%(x2med-mu);T2obs #4662.783

#Significance and if we accept or reject----------
F<-qf(p=1-alpha, df1=p,df2=n-p);F #3.0507
coef<-(p*(n-1)/(n-p));coef  #2.012121
tt<-(p*(n-1)/(n-p))*F;tt  #6.138553

#Rejection region ---------
#RR[[6.13855;+infinity]]

#p-value

df1=p;p
df2=n-p;df2
virsus<-n-p;virsus
apacia<-p*(n-p);apacia
dalyba<-virsus/apacia;dalyba
pvalue_M<-1-pf(q=4662.783*0.5,2,165,lower.tail = TRUE, log.p = FALSE);pvalue_M
format.pval(pvalue_M)



