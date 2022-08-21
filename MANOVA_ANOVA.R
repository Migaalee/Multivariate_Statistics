#ANOVA one way------
#Hypothesis mu1=mu2=...----
#By hand--------
pop1<-c(9,6,9);pop1
pop2<-c(0,2);pop2
pop3<-c(3,1,2);pop3
allpop<-c(9,6,9,0,2,3,1,2);allpop
alpha<-0.01
g<-c(rep(1,3),rep(2,2),rep(3,3));g
x<-as.data.frame(cbind(g,allpop));x
n<-nrow(x);n
mean1<-mean(pop1);mean1
mean2<-mean(pop2);mean2
mean3<-mean(pop3);mean3
mean<-mean(x$allpop);mean
n1<-nrow(x[x$g==1,]);n1
n2<-nrow(x[x$g==2,]);n2
n3<-nrow(x[x$g==3,]);n3
g_size=3
#SQR residual---------

SQR<-sum((pop1-mean1)^2)+sum((pop2-mean2)^2)+sum((pop3-mean3)^2);SQR

#SQF=SST in english---------

SQF<-sum(3*(mean1-mean)^2)+sum(2*(mean2-mean)^2) +sum(3*(mean3-mean)^2);SQF

#SQT total-------

SQT<-SQR+SQF;SQT

#QMT---------

QMT<-SQT/(n-1);QMT

#QMR--------

QMR<-SQR/(n-g_size);QMR


#QMF-------------

QMF<-SQF/(g_size-1);QMF

#Hypothesis test---------
alpha=0.01
df1=g_size-1
df2=n-g_size
F<-QMF/QMR;F
F_critical<-qf(1-alpha,df1,df2);F_critical

pvalue<-1-pf(F,df1,df2);pvalue

#F is more than F_critical, we reject our null hypothesis

#Rejection region

#Fobs(19.5) belongs to RR[F_critical(13.27);+begalybe]







#Using built-in function----
x$g <- as.factor(x$g) # change as factors!!!
levels(x$g)
x$g <- ordered(x$g,levels = c("1", "2", "3")) #ordered data
#Calculate basic statistics like mean, sd
library(dplyr)
group_by(x, g) %>%
  summarise(
    count = n(),
    mean = mean(allpop, na.rm = TRUE),
    sd = sd(allpop, na.rm = TRUE)
  )
# Compute the analysis of variance
res.aov <- aov(allpop ~ g, data = x)
# Summary of the analysis
summary(res.aov)
#As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the groups highlighted with "*" in the model summary.

#Multiple pairwise-comparison between the means of groups--------
#We test it as we dont know which of the group is behaving differently

#Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
#diff: difference between means of the two groups
#lwr, upr: the lower and the upper end point of the confidence interval at 95% (default)
#p adj: p-value after adjustment for the multiple comparisons.
#seems like group 3 and 2 is not significant, while other are

#Another test using generalised linear model hypothesis
install.packages("multcomp")
library(multcomp)
summary(glht(res.aov, linfct = mcp(g = "Tukey")))

#Pairewise t-test

pairwise.t.test(x$allpop, x$g,p.adjust.method = "BH")


#Check the homogeneity of variance assumption-----
plot(res.aov, 1) #we can plot residuals vs fitted values, we can visualise outliers

#Levene test for homogeneity
library(car)
leveneTest(allpop ~ g, data = x)
#Because pvalue is larger than 0.05, This means that there is no evidence to suggest that the variance across groups is statistically significantly different. Therefore, we can assume the homogeneity of variances in the different treatment groups.

#Check the normality assumption-----
plot(res.aov, 2)
#As all the points fall approximately along this reference line, we can assume normality.

plot(res.aov)

# Shapiro-Wilk test for normality--------

# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)
#pvalue is less than 0.05, normality might have been violated

#MANOVA One-way --------
#By hand---------
p1<-matrix(c(9,6,9,3,2,7),3,2);p1
p2<-matrix(c(0,2,4,0),2,2);p2
p3<-matrix(c(3,1,2,8,9,7),3,2);p3
mean1<-colMeans(p1);mean1
mean2<-colMeans(p2);mean2
mean3<-colMeans(p3);mean3

x1<-c(9,6,9,0,2,3,1,2)
x2<-c(3,2,7,4,0,8,9,7)
alpha<-0.01
g<-c(rep(1,3),rep(2,2),rep(3,3));g
x<-as.data.frame(cbind(g,x1,x2));x
n<-nrow(x);n
n1<-nrow(x[x$g==1,]);n1
n2<-nrow(x[x$g==2,]);n2
n3<-nrow(x[x$g==3,]);n3
p<-ncol(x)-1;p #removing one as its a name
s<-cov(x[,-1]);s  #dont forget to remove col name
s1<-cov(x[x$g==1,-1]);s1
s2<-cov(x[x$g==2,-1]);s2
s3<-cov(x[x$g==3,-1]);s3
g<-unique(x$g);g
x$g <- as.factor(x$g) # change as factors!!!
g_size=3

#W is residual
W<-(n1-1)*s1+(n2-1)*s2+(n3-1)*s3;W 
W<-(n-g_size)*Spooled;W  #if we have Spooled

#Total
total<-(n-1)*s;total
total2<-B+W;total2

#B is fatorial
B<-total-W;B


#Wilks-----
#Wilks if we have W and total--------
wilks<-det(W)/det(total);wilks

#Wilks if we have W and B--------
BinvW<-B%*%solve(W);BinvW
eval<-eigen(BinvW)$val;eval
eval[eval!=0] # noce function to see only eval that has no 0 values
wilks2<-prod(1/(1+eval));wilks2


#Hypothesis test p=1,g>=2--------
#how to get wilks?Whats the difference from one-way anova?

df1=(g_size-1);df1
df2=(n-g_size-1);df2

Fobs<-((n-g_size-1)/(g_size-1))*((1-(wilks))/(wilks));Fobs #50.00837
F_critical<-qf(1-alpha,df1,df2);F_critical
pvalue<-1-pf(Fobs,df1,df2);pvalue


#Hypothesis test p=2,g>=2--------
df1=(g_size-1)*2;df1
df2=(n-g_size-1)*2;df2
Fobs<-((n-g_size-1)/(g_size-1))*((1-sqrt(wilks))/sqrt(wilks));Fobs
F_critical<-qf(1-alpha,df1,df2);F_critical
#If Fobs(8.19886) is more than F_critical(7.006077), we reject our null hypothesis
pvalue<-1-pf(Fobs,df1=4,df2=8);pvalue #last number is our pvalue that came back from Wilks built-in function


#MANOVA model p=2,g>=2---------
x[,2] #this is x1
x[,3] #this is x2
x[,1] #this is our groups
modelo<-manova(cbind(x[,2],x[,3])~as.factor(x[,1]));modelo #binding columns and making sure column 1 is factor
summary(modelo,test="Wilks") #Wilks built in function
#num Df = df1 
#den Df = df2
#Df-first is g_size-1=3-1=2, second N=g_size=8(nrow)-g_size=5


#Hypothesis test p>=1,g=2--------

df1=p;df1
df2=n-p-1;df2
Fobs<-((n-p-1)/(p))*((1-(wilks))/(wilks));Fobs
F_critical<-qf(1-alpha,df1,df2);F_critical
#If Fobs(62.51046) is more than F_critical(13.27393), we reject our null hypothesis
pvalue<-1-pf(Fobs,df1=df1,df2=df2);pvalue

#Hypothesis test p>=1,g=3--------

g_s=3
df1=2*p;df1
df2=2*(n-p-2);df2
Fobs<-((n-p-2)/(p))*((1-sqrt(wilks))/sqrt(wilks));Fobs
F_critical<-qf(1-alpha,df1,df2);F_critical
#If Fobs(8.19886) is more than F_critical(7.006077), we reject our null hypothesis
pvalue<-1-pf(Fobs,df1=df1,df2=df2);pvalue


#Hypothesis test with large g and p--------

#p^2 + (n-g)^2 -5 has to be more than 5 to use this formula, else use t=1
g_s=3 #imaginary numbers
n=5
p=2
wilks=wilks
alpha=0.05
v1=p*(g_s-1);v1
w=(n-g_s)-0.5*(p-g_s+2);w
t=sqrt(((p^2)*((g_s-1)^2)-4))/sqrt(((p^2)+((n-g_s)^2)-5));t
v2=w*t-0.5*(p*(g_s-1)-2);v2

Fobs<-(((1-(wilks)^(1/t))*v2)/(wilks)^(1/t))*v1;Fobs #32.79
F_critical<-qf(1-alpha,df1=v1,df2=v2);F_critical #19.246

#Alternatively hypothesis for large n can be calculated:

Fobs<-(-(n-1)-(p+g_s)/2)*log(wilks);Fobs #21.17
Q<-qchisq(1-alpha,df=g_s-1);Q   #5.99


#Roy's largest root------
BinvW<-B%*%solve(W);BinvW
eval<-eigen(BinvW)$val;eval

# eval[eval!=0] # noce function to see only eval that has no 0 values
wilks2<-prod(1/(1+eval));wilks2

Roy<-prod(eval/(1+eval));Roy

#Hotelling´s trace-------
library(matlib)
hotelling_trace<-tr(BinvW);hotelling_trace


#Pillai's trace----------

pilau<-tr(B%*%solve(B+W));pilau


#Variance (meanXij-meanxkj)----
x1<-c(9,6,9,0,2,3,1,2)
x2<-c(3,2,7,4,0,8,9,7)
g<-c(rep(1,3),rep(2,2),rep(3,3));g
x<-as.data.frame(cbind(g,x1,x2));x
n<-nrow(x);n
n1<-nrow(x[x$g==1,]);n1
n2<-nrow(x[x$g==2,]);n2
n3<-nrow(x[x$g==3,]);n3
p<-ncol(x)-1;p #removing one as its a name
s<-cov(x[,-1]);s  #dont forget to remove col name
s1<-cov(x[x$g==1,-1]);s1
s2<-cov(x[x$g==2,-1]);s2
s3<-cov(x[x$g==3,-1]);s3

#Var(meanxij-meanxkj)=var(meanxij)+var(meanXkj)
#This is two sample t-based confidence interval
n=8
Wjj<-diag(W)/(n-3);Wjj
var<-((1/n1)+(1/n2))*(Wjj/(n-3));var #this is var(meanxij-meanxkj)

#This is simultaneous confidence statements
g=3
m<-p*(g*(g-1))/2;m #number of comparisons we will make
alpha=0.05

#Bonferroni simultaneous confidence intervals-----------

quantilt<-qt(1-alpha/(2*m),df=n-3);quantilt #4.21 this -s plusminus
mgrupo<-aggregate(x[,-1],list(x[,1]),mean);mgrupo # groups 1 to 3, mean for each group for variables x1 and x2


#group  1 vs. 2 
#var x1
aIC1plus<-(mgrupo[1,"x1"]-mgrupo[2,"x1"])+quantilt*sqrt(Wjj[1]*(1/3+1/2));aIC1plus
bIC1minus<-(mgrupo[1,"x1"]-mgrupo[2,"x1"])-quantilt*sqrt(Wjj[1]*(1/3+1/2));bIC1minus
#answer (1.5529, 12.4471)
#var x2
cIC2plus<-(mgrupo[1,"x2"]-mgrupo[2,"x2"])+quantilt*sqrt(Wjj[2]*(1/3+1/2));cIC2plus
dIC2minus<-(mgrupo[1,"x2"]-mgrupo[2,"x2"])-quantilt*sqrt(Wjj[2]*(1/3+1/2));dIC2minus
# answer x2(-6.43,10.4386)

#grupo 1 vs. 3 

#var x1
eIC1plus<-(mgrupo[1,"x1"]-mgrupo[3,"x1"])+quantilt*sqrt(Wjj[1]*(1/3+1/3));eIC1plus
fIC1minus<-(mgrupo[1,"x1"]-mgrupo[3,"x1"])-quantilt*sqrt(Wjj[1]*(1/3+1/3));fIC1minus
#answer(1.1280,10.8720)
#var x3
gIC2plus<-(mgrupo[1,"x2"]-mgrupo[3,"x2"])+quantilt*sqrt(Wjj[2]*(1/3+1/3));gIC2plus
hIC2minus<-(mgrupo[1,"x2"]-mgrupo[3,"x2"])-quantilt*sqrt(Wjj[2]*(1/3+1/3));hIC2minus
#answer(3.54773,-11.54773)

#group 2 vs. 3 
#var x2
iIC1plus<-(mgrupo[2,"x1"]-mgrupo[3,"x1"])+quantilt*sqrt(Wjj[1]*(1/3+1/2));iIC1plus
jIC1minus<-(mgrupo[2,"x1"]-mgrupo[3,"x1"])-quantilt*sqrt(Wjj[1]*(1/3+1/2));jIC1minus
#answer(-6.4471,4.447)
#var x3
kIC2plus<-(mgrupo[2,"x2"]-mgrupo[3,"x2"])+quantilt*sqrt(Wjj[2]*(1/3+1/2));kIC2plus
lIC2minus<-(mgrupo[2,"x2"]-mgrupo[3,"x2"])-quantilt*sqrt(Wjj[2]*(1/3+1/2));lIC2minus
#answer (-14.4886,2.4386)


#Two-way MANOVA---------





#Distributions----
