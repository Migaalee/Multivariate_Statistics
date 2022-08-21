
library(HSAUR2)
DATA=heptathlon
head(DATA)

dados<-scale(DATA,center = TRUE,scale=TRUE);dados
head(dados)
S=cov(dados);S
pca_results<-prcomp(dados)

# eigenvalues------------

eigenvals <- pca_results$sdev^2;eigenvals

#eigenvectors(rotation, loading)-----------

eigenvectors<-pca_results$rotation;eigenvectors




#Kaiser criteria--------------------------
Sval<-round(eigen(S)$val,4);Sval
eigenvals <- pca_results$sdev^2;eigenvals

#plot eigenvalues (which are variance)


plot(
  x = seq(1:length(eigenvals)), y = eigenvals,
  type = "o",
  xlab = "Principle Component", ylab = "Eigenvalues(variance)")


plot(pca_results,type="l")


fviz_eig(pca_results, choice = "eigenvalue", 
         addlabels=TRUE)


##Scree plot-----------

fviz_eig(pca_results, 
         addlabels=TRUE)


#proportion of variability from retained PC1 and PC2------
num=eigen(S)$val;num
den=sum(eigen(S)$val);den
percentage_variab<-round(num/den,3);percentage_variab

#0.681 0.150 0.065 0.057 0.031 0.009 0.006 0.001


#Variation for each individual(or genes)--------------------
y=dados%*%pca_results$rotation;y
first_individual<-pca_results$x[1,];first_individual #variation for first individual in a list

head(y[,1])
signif(sort(y[,1],decreasing=TRUE)[1:5],2)


#Adequality of PCA test-----

#a)Bivariate correlations-----

R=cov(dados);R #because dados are standardised we can use cov
round(R,2)
(sum(R>=0.5)-6)/(nrow(S)^2-6) #proportion of correlations that are larger than 0.5
#0.31 1/3rd of all correlations are well correlated


#b)Parcial correlations----

library(corpcor)

round(cor2pcor(R),2)
#smallest values of partial correlations will be best???

pR=(sum(R<0.5)-6)/(nrow(S)^2-6);pR
#0.5 our threshold we decide
#pR=0.5862 



#Mauchly's test of sphericity-------------
#H0 : S = sigma^2xI vs. H1 :S neg= sigma^2xI
#when pvalue is small we reject null hypothesis, we can use data for PCA analysis.
#that also means at least some variables have covariance that is not equal to 0.

dados1<-as.matrix(dados)
mauchly.test(lm(dados1~1))

#pvalue p-value < 2.2e-16
#W = 5.5768e-06


#Mauchlys hypothesis test-------------

n=nrow(dados);n  #25

p=ncol(dados);p  #note, check that name is not counted as column #8


#W(U)=p^p|S|/(tr(S))^p
S<-cov(dados);S
det_S<-det(S);det_S
trS<- sum(diag(S));trS
W<-(p^p)%*%det(S)/trS^p;W #same value as in Mauchlys test 5.576783e-06

#lambda    lambda=|S|^(n/2)/(tr(S)/p)^(n*p/2)

lam<-(det_S^(n/2))/(trS/p)^(n*p/2);lam

#U=labmbda^(2/n)

W1<-lam^(2/n);W1  #same result 5.576783e-06


#test

#U_ast=-(n-1-((2p^2+p+2)/(6*p))*ln(U)

U_ast=-(n-1-(2*p^2+p+2)/(6*p))*log(U);U_ast

#pvalue--------

pval<-pchisq(q=U_ast,df=35);pval  
format.pval(pval)  # < 2.22e-16

#distribution---------

n=nrow(dados);n  #25

p=ncol(dados);p  #note, check that name is not counted as column #8

df=0.5*p*(p+1)-1;df #35
#dstribution X^2 with df=35









#Correlation between variables-----------------
#non-standardized data

correlation <-cor(DATA,DATA);correlation

#scaled data

dados<-scale(DATA,center = TRUE,scale=TRUE);dados
correlation1 <-cor(dados,dados);correlation1




#correlation between variables xj and PCk------

head(DATA)

#non-standardized variables

pca_non_st<-prcomp(DATA)

new_var_n<-pca_non_st$x;new_var_n
old_var1<-DATA;old_var1

correlation<-round(cor(new_var_n,old_var1),2);correlation


#standardized variables

dados<-scale(DATA,center = TRUE,scale=TRUE);dados
pca_st<-prcomp(dados)

new_var_st<-pca_st$x;new_var_st
old_var2<-dados;old_var2


correlation2<-round(cor(new_var_st,old_var2),2);correlation2

#|rkj|>=0.5 or rkj^2>=0.25 are elevated correlations



#Plot all correlations----------------------
plot(heptathlon[,-8])



#linear combinations (score)------------
pca_results$rotation 

#The elements of an eigenvector are the weights aij, and are known as loadings

#rotation = the matrix of variable loadings (columns are eigenvectors)

#from up to down

#hurdles   0.40734940 -0.17794781 -0.043796357 -0.025041525  0.112963713 -0.794814407
#highjump -0.33892709  0.26354914  0.367070576  0.678862033  0.005560801 -0.080264774
#shot     -0.33224949 -0.27292253 -0.677276565  0.123452280  0.493920625  0.077124708
#run200m   0.37039666  0.24019801  0.085500579  0.363014027  0.664656345  0.009835088
#longjump -0.41221153  0.07586727 -0.140635234  0.109886198 -0.197835144 -0.572166386
#javelin  -0.07767939 -0.83595531  0.470652388  0.119717205  0.123954423  0.039480775
#run800m   0.33775722 -0.23949840 -0.394971096  0.604197195 -0.491178981  0.132793809
#score    -0.42612100 -0.06897840

#y1=0.40734940x1 -0.33892709

#y1 is called scores



#Confidence intervals for PCs (eigenvalues)---------------

#THIS IS CI for PC1 (eigenvalues) only from page 7/28 aula 10

eval=pca_results$sdev^2;eval   #variance <- (pca$sdev)^2 
eval1=eval[1]  #variance
n=nrow(dados);n
alpha=0.05

L1=eval1/(1+qnorm(1-alpha/2)*sqrt(2/(n-1)))
U1=eval1/(1-qnorm(1-alpha/2)*sqrt(2/(n-1)))
IC1=c(L1,U1);IC1
#3.478281 12.543017

L2=eval1-qnorm(1-alpha/2)*eval1*sqrt(2/(n-1))
U2=eval1+qnorm(1-alpha/2)*eval1*sqrt(2/(n-1))
IC2=c(L2,U2);IC2
#2.364808 8.527727

L3=log(eval1)-qnorm(1-alpha/2)*sqrt(2/(n-1))
U3=log(eval1)+qnorm(1-alpha/2)*sqrt(2/(n-1))
IC3=c(exp(L3),exp(U3));IC3
#3.092988 9.590024


#plot results
plot(IC1,c(1,1),type="b",ylim=c(1,3),xlim=c(1,5))
lines(IC2,c(2,2),type="b",col=2)
lines(IC3,c(3,3),type="b",col=4)
abline(v=eval1,lty=2)


#Distribution of eigenvalues N(0,2) 5/28 aula10-----

#Distribution of eigenvectors X2(df=p-1) 6/28 aula10-----