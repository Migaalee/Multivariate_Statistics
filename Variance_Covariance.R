#Matrix of observations

x<-matrix(c(4,6,5,6,5,4,3,4,3,3,4,5),4,3);x
n<-nrow(x);n
p<-ncol(x);p


#Euclidean distance between 2 points--------------
x<-matrix(c(2,4,3,2,4,3,4,2),4,2);x
x1<-matrix(c(2,4,3,2),4,1);x1
x2<-matrix(c(4,3,4,2),4,1);x2
euclidean_x1X2<-sqrt(sum((x1-x2)^2));euclidean_x1X2

#Euclidean distance for coordinates p--------------
d<-sqrt(sum(x1^2+x2^2));d


#Calculate variance/covariance matrix S--------
S<-cov(x);S

#Calculate mean---------

x_mean<-colMeans(x);x_mean


#Calculate mean using matrices formula-------------


matrixOne_x<-matrix(c(1), nrow =n);matrixOne_x
x_mean2<-(1/n)*t(x)%*%matrixOne_x;x_mean2 


#Using identity matrix to calculate S------------

x1<-matrix(c(42,52,48,58,4,5,4,3),4,2);x1
n1<-nrow(x1);n1
I<-diag(nrow = n1);I
matrixOne<-matrix(c(1), nrow =n1,ncol=n1);matrixOne
S1<-1/(n1-1)*t(x1)%*%(I-(1/n1)*matrixOne)%*%x1;S1


#Calculate generalized variance-----------
genS<-det(S);genS
genS2<-prod(eigen(S)$val);genS2


#calculate total variance--------------

trS<-sum(eigen(S)$val);trS

#Calculate correlation matrix R if we have X-------------

R<-cor(x); R  # if we have data X


#Calculate correlation matrix If we only have S---------------
SD<-diag(x=diag(S)^(1/2),nrow = p, ncol=p);SD #Standard dev matrix
R1<-solve(SD)%*%S%*%solve(SD);R1

#Calculate R determinant------------

R_<-prod(eigen(R)$val);R_



#Calculate S If we only have correlation matrix R and standard deviation matrix S-----------

S1<-SD%*%R1%*%SD;S1



#Linear combinations----------------

#matrix 3 rows, 5 columns
x_3<-matrix(c(51,27,37,36,20,22,50,26,41,35,17,37,42,27,30),3,5);x_3
x_4<-matrix(c(2,4,3,2,4,3,4,2),4,2);x_4

#vectors z = 3y1 - 2y2 + 4y3 - y4 + y5 
z<-matrix(c(3,-2,4,-1,1),5,1);z
u<-matrix(c(3,2),2,1);u

#vector w = y1 - 3y2 - y3 + y4 - 2y5.

w<-matrix(c(1,-3,-1,1,-2),5,1);w
v<-matrix(c(2,4),2,1);v


#Find Sz, Sw, Suv,Ruv, correlation between z and w, u and v CORRECT ONE----

med_y<-colMeans(x_4);med_y

S_y<-var(x_4);S_y

Su<-t(u)%*%S_y%*%u; Su
Sv<-t(v)%*%S_y%*%v;Sv
Suv<-t(v)%*%S_y%*%u; Suv
ruv<-Suv/sqrt(Su*Sv); ruv


#Another exercise
x<-matrix(c(4,6,5,6,5,4,3,4,3,3,4,5),4,3);x

u<-matrix(c(1,2,-1),3,1);u
v<-matrix(c(-2,1,-1),3,1);v



#Calculate mean of new vector u = x1 +2x2 -x3

u_mean<-t(u)%*%x_mean;u_mean


#Calculate variance S of new vector u = x1 +2x2 -x3  


var_u<-t(u)%*%S%*%u;var_u

#Calculate variance/covariance S of new vector v = -2x1 +x -x3 

var_v<-t(v)%*%S%*%v;var_v


#Correlation coeficient Pearsons

data1<-c(2,3,4,4);data1
data2<-c(2,3,4,5);data1

pearson<-cor(data1,data2);pearson #0.94 highly correlated


#Spectral Decomposition Of A Matrix

#only symmetric matrix


#Is matrix positive or negative?

x1<-matrix(c(-1,2,5,3,4,2,-2,2,3),3,3) #square matrix
Xposneg<-eigen(x1)$val;Xposneg  #negative as it has negative values


#Determinant only for square matrix 

x1_det<-det(x1);x1_det
x1_det<-det(x1);x1_det


#are vectors linearly independent?

u<-matrix(c(1,2,-1),3,1);u
v<-matrix(c(-2,1,-1),3,1);v

#they are if covariance between them is 0 (not diagonal, but other)

uv<-matrix(c(u,v),2,3,byrow=TRUE);uv
S_uv<-cov(uv);S_uv











  

