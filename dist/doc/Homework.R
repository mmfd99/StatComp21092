## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- test-plot---------------------------------------------------------------
plot(1)         # high-level plot
abline(0, 1)    # low-level change
plot(rnorm(10)) # high-level plot
# many low-level changes in a loop (a single R expression)
for(i in 1:10) {
    abline(v = i, lty = 2)
}

## -----------------------------------------------------------------------------
library(knitr)
res<-data.frame(sep=1:12,name=LETTERS[1:12],month=month.abb[1:12])
knitr::kable(res)

## -----------------------------------------------------------------------------
n<-1e6
U<-runif(n)
par(mfrow=c(1,2))
#sigma=1
x<-sqrt(-2*log(1-U))
hist(x,probability = TRUE,main = "sigma=1")
y<-seq(0,5,0.001)
lines(y,y*exp(-y^2/(2*1*1))/(1*1))
#sigma=10
x<-sqrt(-2*100*log(1-U))
hist(x,probability = TRUE,main = "sigma=10")
y<-seq(0,50,0.001)
lines(y,y*exp(-y^2/(2*10*10))/(10*10))

## -----------------------------------------------------------------------------
n<-1e6
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
par(mfrow=c(1,2))
#p=0.75
p<-0.75
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z, breaks=10000,main = "p=0.75")
#p=0.5
p<-0.5
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z, breaks=10000,main = "p=0.5")

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
#p=0.25
p<-0.25
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z, breaks=10000,main = "p=0.25")
#p=0.375
p<-0.375
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z,breaks=10000, main = "p=0.375")

## -----------------------------------------------------------------------------
x<-seq(0,3,0.001)
y1<-x*exp(-x*x/2)
y2<-(x-3)*exp(-(x-3)*(x-3)/2)
y<-y2/(y2-y1)
plot(x,y)

## -----------------------------------------------------------------------------
n<-1e6
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
par(mfrow=c(1,2))
#p=0.2
p<-0.2
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z, breaks=10000,main = "p=0.2")
#p=0.375
p<-0.8
k<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
Z<-k*X1+(1-k)*X2
hist(Z,breaks=10000, main = "p=0.8")

## -----------------------------------------------------------------------------
# lambda=1,alpha=1,beta=1, we have E[Y]=1, E[Y^2]=2, the results should be E[X]=10,Var[X]=20
lambda <- 1
t0 <- 10
n <-1e5
X <- 0
for(k in 1:n)
{
  N<-rpois(1,lambda*t0)
  Yk<-rgamma(N,1)
  X[k]<-sum(Yk)
}
mean(X)
var(X)

# lambda=10,alpha=2,beta=0.5,we have E[Y]=4,E[Y^2]=24,the results should be E[X]=400,Var[X]=2400
lambda <- 10
t0 <- 10
n <-1e5
X <- 0
for(k in 1:n)
{
  N<-rpois(1,lambda*t0)
  Yk<-rgamma(N,2,scale = 2)
  X[k]<-sum(Yk)
}
mean(X)
var(X)


## -----------------------------------------------------------------------------
MC_beta<-function(x)#生成一个函数来模拟beta(3,3)
{
  n<-1e6
  y<-runif(n,0,x)#生成n个(0,x)上均匀分布的随机数
  theta=x*mean(y*y*(1-y)*(1-y)/beta(3,3))
  theta
}
for (i in 1:9) {
  x=i/10
  print(c(MC_beta(x),pbeta(x,3,3)))
}

## -----------------------------------------------------------------------------
#MC模拟cdf
MC_Rayleigh<-function(x, m=1e5, antithetic=TRUE)
{
  u<-runif(m/2)
  if(!antithetic)v<-runif(m/2)else
    v<-1-u
  u<-c(u,v)
  cdf<-numeric(length(x))
  for(i in 1:length(x))
  {
    g<-x[i]*u*x[i]*exp(-(u*x[i])^2/2)
    cdf[i]<-mean(g)
  }
  cdf
}
par(mfrow=c(1,2))
x<-seq(0,3,length=100)
MC1<-MC_Rayleigh(x)
MC2<-MC_Rayleigh(x,antithetic = FALSE)
plot(x,MC1)
plot(x,MC2)

## -----------------------------------------------------------------------------
#下面计算方差的减少
m<-1e2
MC1<-MC2<-numeric(m)
for (j in 1:5) {
  x=3/j
  for(i in 1:m)
  {
    MC1[i]=MC_Rayleigh(x)
    MC2[i]=MC_Rayleigh(x,antithetic = FALSE)
  }
  print(c(x,(var(MC2)-var(MC1))/var(MC2)))
}

## -----------------------------------------------------------------------------
m<-1e3
g<-function(x)
{
  x*x*exp(-x*x/2)/sqrt(2*pi)
}

drayleigh = function (x, sigma) {
  stopifnot(sigma > 0)
  y = x / (sigma^2) * exp(- x^2 / (2 * sigma^2))
  y[x < 0] = 0
  y
}

xs = seq(0,10,0.1)

ys.g = g(xs)
#f1用rayleigh分布
ys.rayleigh = drayleigh(xs, sigma = 1)
#f2用正态分布
ys.norm = dnorm(xs, mean = 1)
lim = max(c(ys.g, ys.rayleigh, ys.norm))

plot(xs, ys.g, type = "l", ylim = c(0, lim))
lines(xs, ys.rayleigh, col="red", ylim = c(0, lim))
lines(xs, ys.norm, col="blue", ylim = c(0, lim))


## -----------------------------------------------------------------------------
mean = 1
n = 1e5
g = function (x) {
  x ^ 2 / sqrt(2*pi) * exp(-x^2/2) * (x > 1)
}
f2 = function (x) {
  dnorm(x, mean = mean) * (x > 1)
}
rf2 = function () {
  rnorm(n, mean = mean)
}
is.norm = function () {
  xs = rf2()
  return(mean(g(xs)/f2(xs), na.rm = TRUE)/2)  
}
theta2 = is.norm()
theta2

## -----------------------------------------------------------------------------
m=1e4
n=20
y=numeric(m)
for (j in 1:m) 
{
  x=rchisq(n,2)  
  lci=t.test(x)[4]$conf.int[1]
  rci=t.test(x)[4]$conf.int[2]
  if (lci<2 & 2<rci)
    y[j]=1
}
sum(y)/m

## -----------------------------------------------------------------------------
m=1e4
n=20
I=numeric(m)
L=numeric(m)
# 自由度为1的卡方分布
for (j in 1:m) 
{
  x=rchisq(n,1)-1 
  if(t.test(x)[3]<=0.05)
    I[j]=1#reject mu=mu0
}
sum(I)/m

# (0,2)上的均匀分布
I=numeric(m)
for (j in 1:m) 
{
  x=runif(n,0,2)-1
 if(t.test(x)[3]<=0.05)
    I[j]=1#reject mu=mu0
}
sum(I)/m
# 均值为1的指数分布
I=numeric(m)
for (j in 1:m) 
{
  x=rexp(n,1)-1
  if(t.test(x)[3]<=0.05)
    I[j]=1#reject mu=mu0
}
sum(I)/m

## -----------------------------------------------------------------------------
#计算b_1d
t1=proc.time()
b<-function(n)
{
  library(MASS)
  mean=c(0,0)
  sigma=matrix(c(1,0,0,1),nrow=2,ncol=2)
  X=mvrnorm(n,mean,sigma)
  barX=c(mean(X[,1]),mean(X[,2]))
  #计算方差阵的极大似然
  DX=matrix(0,nrow=n,ncol=2)
  DX[,1]=X[,1]-barX[1]
  DX[,2]=X[,2]-barX[2]
  Sigma=matrix(0,nrow=2,ncol=2)
  Sigma=t(DX)%*%DX/n
  #计算偏度
  b1d=mean((DX%*%solve(Sigma)%*%t(DX))^3)
    
  b1d
}

# 生成渐进分布
n<-c(100,400,900,1600,2500,3600,4900,6400,8100,10000)
cv<-6*qchisq(0.975,2)/n

#比较

# p.reject=numeric(length(n)) 
# m=100
# for (i in 1:length(n))
# {
#   sktests <- numeric(m)
#   for (j in 1:m)
#   {
# 
#     sktests[j] <- as.integer(abs(b(n[i])) >= cv[i] )
#   }
#   p.reject[i] <- mean(sktests) #proportion rejected
# }
# t2=proc.time()
# p.reject
# t2-t1

## -----------------------------------------------------------------------------
alpha=0.1
n=100
m=1000
epsilon=c(seq(0, .15, .01), seq(.15, 1, .05))
N=length(epsilon)
pwr=numeric(N)

cv=6*qchisq(0.95,2)/n

b1<-function(n,e)
{
  library(MASS)
  X=matrix(0,nrow=n,ncol=2)
  mean=c(0,0)
  sigma <- sample(c(1, 10), replace = TRUE, size = n, prob = c(1-e, e))
  for(i in 1:n)
    X[i,]=mvrnorm(n=1,mu=mean,Sigma=matrix(c(sigma[i],0,0,sigma[i]),ncol=2,nrow=2 ))
  
  barX=c(mean(X[,1]),mean(X[,2]))
  #计算方差阵的极大似然
  DX=matrix(0,nrow=n,ncol=2)
  DX[,1]=X[,1]-barX[1]
  DX[,2]=X[,2]-barX[2]
  Sigma=matrix(0,nrow=2,ncol=2)
  Sigma=t(DX)%*%DX/n
  #计算偏度
  b1d=mean((DX%*%solve(Sigma)%*%t(DX))^3)
    
  b1d
}

# for (j in 1:N)
# {
#   e <- epsilon[j]
#   sktests <- numeric(m)
#   for (i in 1:m)
#   {
#     sktests[i] <- as.integer(abs(b1(n,e)) >= cv)
#   }
#   pwr[j] <- mean(sktests)
# }
# #plot power vs epsilon
# plot(epsilon, pwr, type = "b",
# xlab = bquote(epsilon), ylim = c(0,1))
# abline(h = .1, lty = 3)
# se <- sqrt(pwr * (1-pwr) / m) #add standard errors
# lines(epsilon, pwr+se, lty = 3)
# lines(epsilon, pwr-se, lty = 3)


## -----------------------------------------------------------------------------
library(bootstrap)
hat_lambda=eigen(cov(scor))$values
hat_theta=hat_lambda[1]/sum(hat_lambda)
hat_theta

## -----------------------------------------------------------------------------
#估计bias
m=1e3
n=nrow(scor)
theta=numeric(m)

for(i in 1:m)
{
  index=sample(1:n,size = n,replace = TRUE)
  dat=scor[index,]
  lambda=eigen(cov(dat))$values
  theta[i]=lambda[1]/sum(lambda)
}

mean(theta-hat_theta)

## -----------------------------------------------------------------------------
#估计se
m=1e3
n=nrow(scor)
theta=numeric(m)

for (i in 1:m)
{
  index=sample(1:n,size = n,replace = TRUE)
  dat=scor[index,]
  lambda=eigen(cov(dat))$values
  theta[i]=lambda[1]/sum(lambda)
}

sd(theta)

## -----------------------------------------------------------------------------
#估计bias
n=nrow(scor)
theta=numeric(n)
for(i in 1:n)
{
  dat=scor[-i,]
  lambda=eigen(cov(dat))$values
  theta[i]=lambda[1]/sum(lambda)
}
(n-1)*(mean(theta)-hat_theta)

## -----------------------------------------------------------------------------
#估计se
n=nrow(scor)
theta=numeric(n)
for(i in 1:n)
{
  dat=scor[-i,]
  lambda=eigen(cov(dat))$values
  theta[i]=lambda[1]/sum(lambda)
}
sqrt((n-1)*mean((theta-mean(theta))^2))

## -----------------------------------------------------------------------------
#计算percentile CI
library(boot)
data(scor, package = "bootstrap")

theta_boot=function(da,index)
{
  dat=da[index,]
  lambda=eigen(cov(dat))$values
  lambda[1]/sum(lambda)
}

boot.ci(boot(scor,statistic = theta_boot,R=1e3),type="perc")

## -----------------------------------------------------------------------------
#计算BCa CI
boot.ci(boot(scor,statistic = theta_boot,R=1e3),type="bca")

## -----------------------------------------------------------------------------
sk=function(x,index)
{
  xbar=mean(x[index])
  m3=mean((x[index] - xbar)^3)
  m2=mean((x[index] - xbar)^2)
  m3/m2^1.5
}

## -----------------------------------------------------------------------------
#norm，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)

for(j in 1:m)
{
  x=rnorm(n)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "norm")[4]$normal[2]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "norm")[4]$normal[3]
  if(0<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<0)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#basic，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)

for(j in 1:m)
{
  x=rnorm(n)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "basic")[4]$basic[4]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "basic")[4]$basic[5]
  if(0<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<0)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#perc，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)

for(j in 1:m)
{
  x=rnorm(n)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "perc")[4]$percent[4]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "perc")[4]$percent[5]
  if(0<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<0)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#norm，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)
skew=2*sqrt(10)/5
for(j in 1:m)
{
  x=rchisq(n,df=5)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "norm")[4]$normal[2]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "norm")[4]$normal[3]
  if(skew<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<skew)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#basic，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)

for(j in 1:m)
{
  x=rchisq(n,df=5)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "basic")[4]$basic[4]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "basic")[4]$basic[5]
  if(skew<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<skew)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#perc，输出结果依次为cp，left_pr，right_pr
library(boot)
m=1e2
n=1e1
y=right=left=numeric(m)

for(j in 1:m)
{
  x=rchisq(n,df=5)
  theta1=boot.ci(boot(x,statistic =sk,R=1e3),type = "perc")[4]$percent[4]
  theta2=boot.ci(boot(x,statistic =sk,R=1e3),type = "perc")[4]$percent[5]
  if(skew<theta1)
  {
    left[j]=1
  }
  else 
    if(theta2<skew)
    {
      right[j]=1
    }
    else
      y[j]=1
  
}
print(mean(y))
print(mean(left))
print(mean(right))

## -----------------------------------------------------------------------------
#生成两个随机样本
m=1e2
x=rnorm(m,mean=0,sd=1)
y=rnorm(m,mean=100,sd=1)

#写出Spearman rank 相关系数统计量
spearman=function(x,y)
{
  n=length(x)
  d=rank(x)-rank(y)
  rho=1-(6*sum(d^2))/(n^3-n)
  rho
}

#利用Spearman rank做置换检验

library(boot)
boot.obj=boot(data=rbind(x,y),statistic=spearman,R=1000,sim="permutation")

tb=c(boot.obj$t0,boot.obj$t)

mean(tb>=boot.obj$t0)

#利用cor.test
cor.test(x,y,method = "spearman")



## -----------------------------------------------------------------------------
n1=10
n2=8
x=rnorm(n1,0,1)
y=runif(n2,-10,10)
z=c(x,y)
# NN method
library(RANN)
Tn=function(z,index,sizes,k)
{
  n1=sizes[1]
  n2=sizes[2]
  n=n1+n2
  if(is.vector((z)))
    z=data.frame(z,0)
  z=z[index,]
  NN=nn2(data=z,k=k+1)
  block1=NN$nn.idx[1:n1,-1]
  block2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(block1<n1+0.5)
  i2=sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}
set.seed(1)
N=c(n1,n2)
library(boot)
eqdist.nn=function(z,sizes,k)
{
  boot.obj=boot(data=z,statistic=Tn,R=2000,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  mean(ts>=ts[1])
}
eqdist.nn(z,N,k)
# Energy method
library(energy)
boot.obs=eqdist.etest(z,sizes=c(n1,n2),R=2000)
boot.obs$p.value

# Ball method 
library(Ball)
bd.test(x,y,num.permutations = 2000)$p.value

#Power Comparison
m=1e2
k=3
p=2
set.seed(1)
R=2e3
p.values <- matrix(NA,m,3)

# for(i in 1:m)
# {
#   x=rnorm(n1,0,1)
#   y=runif(n2,-10,10)
#   z=c(x,y)
#   p.values[i,1]=eqdist.nn(z,N,k)
#   p.values[i,2]=eqdist.etest(z,N,R=R)$p.value
#   p.values[i,3]=bd.test(x,y,num.permutations = R)$p.value
# }
# alpha=0.1
# pow=colMeans(p.values<alpha)
# names(pow)=c("NN","Energy","Ball")
# pow

## -----------------------------------------------------------------------------
n1=10
n2=8
x=rnorm(n1,0,1)
y=runif(n2,10,20)
z=c(x,y)
# NN method
library(RANN)
Tn=function(z,index,sizes,k)
{
  n1=sizes[1]
  n2=sizes[2]
  n=n1+n2
  if(is.vector((z)))
    z=data.frame(z,0)
  z=z[index,]
  NN=nn2(data=z,k=k+1)
  block1=NN$nn.idx[1:n1,-1]
  block2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(block1<n1+0.5)
  i2=sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}
set.seed(1)
N=c(n1,n2)
library(boot)
eqdist.nn=function(z,sizes,k)
{
  boot.obj=boot(data=z,statistic=Tn,R=2000,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  mean(ts>=ts[1])
}
eqdist.nn(z,N,k)
# Energy method
library(energy)
boot.obs=eqdist.etest(z,sizes=c(n1,n2),R=2000)
boot.obs$p.value

# Ball method 
library(Ball)
bd.test(x,y,num.permutations = 2000)$p.value

#Power Comparison
m=1e2
k=3
p=2
set.seed(1)
R=2e3
p.values <- matrix(NA,m,3)

# for(i in 1:m)
# {
#   x=rnorm(n1,0,1)
#   y=runif(n2,10,20)
#   z=c(x,y)
#   p.values[i,1]=eqdist.nn(z,N,k)
#   p.values[i,2]=eqdist.etest(z,N,R=R)$p.value
#   p.values[i,3]=bd.test(x,y,num.permutations = R)$p.value
# }
# alpha=0.1
# pow=colMeans(p.values<alpha)
# names(pow)=c("NN","Energy","Ball")
# pow

## -----------------------------------------------------------------------------
n1=10
n2=8
x=rt(n1,1)
X1<-rnorm(n2,0,1)
X2<-rnorm(n2,10,1)
k<-sample(c(0,1),n2,replace = TRUE,prob = c(0.5,0.5))
y<-k*X1+(1-k)*X2
z=c(x,y)
# NN method
library(RANN)
Tn=function(z,index,sizes,k)
{
  n1=sizes[1]
  n2=sizes[2]
  n=n1+n2
  if(is.vector((z)))
    z=data.frame(z,0)
  z=z[index,]
  NN=nn2(data=z,k=k+1)
  block1=NN$nn.idx[1:n1,-1]
  block2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(block1<n1+0.5)
  i2=sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}
set.seed(1)
N=c(n1,n2)
library(boot)
eqdist.nn=function(z,sizes,k)
{
  boot.obj=boot(data=z,statistic=Tn,R=2000,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  mean(ts>=ts[1])
}
eqdist.nn(z,N,k)
# Energy method
library(energy)
boot.obs=eqdist.etest(z,sizes=c(n1,n2),R=2000)
boot.obs$p.value

# Ball method 
library(Ball)
bd.test(x,y,num.permutations = 2000)$p.value

#Power Comparison
m=1e2
k=3
p=2
set.seed(1)
R=2e3
p.values <- matrix(NA,m,3)

# for(i in 1:m)
# {
#   x=rt(n1,1)
#   X1<-rnorm(n2,0,1)
#   X2<-rnorm(n2,10,1)
#   k<-sample(c(0,1),n2,replace = TRUE,prob = c(0.5,0.5))
#   y<-k*X1+(1-k)*X2
#   z=c(x,y)
#   p.values[i,1]=eqdist.nn(z,N,k)
#   p.values[i,2]=eqdist.etest(z,N,R=R)$p.value
#   p.values[i,3]=bd.test(x,y,num.permutations = R)$p.value
# }
# alpha=0.1
# pow=colMeans(p.values<alpha)
# names(pow)=c("NN","Energy","Ball")
# pow

## -----------------------------------------------------------------------------
n1=20
n2=4
x=rt(n1,1)
X1<-rnorm(n2,0,1)
X2<-rnorm(n2,10,1)
k<-sample(c(0,1),n2,replace = TRUE,prob = c(0.5,0.5))
y<-k*X1+(1-k)*X2
z=c(x,y)
# NN method
library(RANN)
Tn=function(z,index,sizes,k)
{
  n1=sizes[1]
  n2=sizes[2]
  n=n1+n2
  if(is.vector((z)))
    z=data.frame(z,0)
  z=z[index,]
  NN=nn2(data=z,k=k+1)
  block1=NN$nn.idx[1:n1,-1]
  block2=NN$nn.idx[(n1+1):n,-1]
  i1=sum(block1<n1+0.5)
  i2=sum(block2>n1+0.5)
  (i1+i2)/(k*n)
}
set.seed(1)
N=c(n1,n2)
library(boot)
eqdist.nn=function(z,sizes,k)
{
  boot.obj=boot(data=z,statistic=Tn,R=2000,sim="permutation",sizes=N,k=3)
  ts=c(boot.obj$t0,boot.obj$t)
  mean(ts>=ts[1])
}
eqdist.nn(z,N,k)
# Energy method
library(energy)
boot.obs=eqdist.etest(z,sizes=c(n1,n2),R=2000)
boot.obs$p.value

# Ball method 
library(Ball)
bd.test(x,y,num.permutations = 2000)$p.value

#Power Comparison
m=1e2
k=3
p=2
set.seed(1)
R=2e3
p.values <- matrix(NA,m,3)

# for(i in 1:m)
# {
#   x=rt(n1,1)
#   X1<-rnorm(n2,0,1)
#   X2<-rnorm(n2,10,1)
#   k<-sample(c(0,1),n2,replace = TRUE,prob = c(0.5,0.5))
#   y<-k*X1+(1-k)*X2
#   z=c(x,y)
#   p.values[i,1]=eqdist.nn(z,N,k)
#   p.values[i,2]=eqdist.etest(z,N,R=R)$p.value
#   p.values[i,3]=bd.test(x,y,num.permutations = R)$p.value
# }
# alpha=0.1
# pow=colMeans(p.values<alpha)
# names(pow)=c("NN","Energy","Ball")
# pow

## -----------------------------------------------------------------------------
set.seed(1234)
N=2e4
M_H=function(N)
{
  x=numeric(N)
  x[1]=rnorm(1)
  k=0
  for (i in 2:N) 
  {
    xt=x[i-1]
    y=rnorm(1,xt)
    u=runif(1)
    r=dcauchy(y)*dnorm(xt,y)/(dcauchy(xt)*dnorm(y,xt))
    if(u<=r)
      x[i]=y
    else
    {
      k=k+1
      x[i]=xt
    }
  }
  x
}


## -----------------------------------------------------------------------------
index=1001:N
plot(index,M_H(N)[index])

## -----------------------------------------------------------------------------
index2=5e3:N
p=seq(0.1,0.9,0.1)
qcauchy(p)
quantile(M_H(N)[index2],probs = p)

## -----------------------------------------------------------------------------
N=1e4           #length of chain
burn=2e3        #burn


G=function(N)
{
  X=matrix(0,N,2) #chain
  a=2
  b=3
  n=100
  X[1,]=c(0.5,0.5)    #initial

  for(i in 2:N)
  {
    x2=X[i-1,2]
    X[i,1]=rbinom(1,n,x2)
    x1=X[i,1]
    X[i,2]=rbeta(1,x1+a,n-x1+b)
  }
  X
}


x=G(N)[(burn+1):N,]

colMeans(x)
cov(x)
cor(x)
par(mfrow=c(1,2))
hist(x[,1],breaks = 20)
hist(x[,2],breaks = 20)

## -----------------------------------------------------------------------------
G.R=function(psi)
{
  psi=as.matrix(psi)
  n=ncol(psi)
  k=nrow(psi)
  
  psi.means=rowMeans(psi)
  B=n*var(psi.means)
  psi.w=apply(psi, 1, "var")
  W=mean(psi.w)
  v.hat=W*(n-1)/n+(B/n)
  r.hat=v.hat/W
  r.hat
}
set.seed(10)
n=3e4
k=4
X=matrix(0,nrow = k,ncol = n)
for (i in 1:k) 
{
  X[i,]=M_H(n)
}

psi=t(apply(X, 1, cumsum))
for(i in 1:nrow(psi))
  psi[i,]=psi[i,]/(1:ncol(psi))
print(G.R(psi))

## -----------------------------------------------------------------------------
set.seed(19)
N=6e5
k=4
X=matrix(0,nrow = k,ncol = N)
Y=matrix(0,nrow = k,ncol = N)
for (i in 1:k) 
{
  X[i,]=G(N)[,1]
  Y[i,]=G(N)[,2]
}

psi.x=t(apply(X, 1, cumsum))
psi.y=t(apply(Y, 1, cumsum))
for(i in 1:nrow(psi.x))
  psi.x[i,]=psi.x[i,]/(1:ncol(psi.x))
print(G.R(psi.x))
for(i in 1:nrow(psi.y))
  psi.y[i,]=psi.y[i,]/(1:ncol(psi.y))
print(G.R(psi.y))

## -----------------------------------------------------------------------------
k_term=function(a,k)
{
  d=length(a)
  (-1)^k*norm(a,"2")^(2*k+2)*gamma((d+1)/2)*gamma(k+1.5)/(factorial(k)*2^k*(2*k+1)*(2*k+2)*gamma(k+1+d/2))
}

## -----------------------------------------------------------------------------
xk=function(a,x)
{
  d=length(a)
  k=floor(x)
  p1=(2*k+2)*log(norm(a,"2"))+log(gamma((d+1)/2))+log(gamma(k+1.5))
  p2=log(factorial(k))+k*log(2)+log(2*k+1)+log(2*k+2)+log(gamma(k+1+d/2))
  return((x>=k)*(x<k+1)*(-1)^k*exp(p1-p2))
}

s=function(a)
{
  integrate(xk,lower=0,upper=140,a=a)
}

## -----------------------------------------------------------------------------
s(matrix(c(1,2)))

## -----------------------------------------------------------------------------
intersection = function (k) 
{
  sn = function (a,n) 
  {
    1-pt(sqrt(a^2 * n / (n+1 - a^2)), df = n)
  }
  
  equa = function (a) {
    sn(a,k) - sn(a,k-1)
  }
  
  e = 1e-6
  uniroot(equa, interval = c(e, sqrt(k)-e))$root
}
k=c(4,25,100,500,1000)
for (i in k) {
  print(intersection(i))
}

## -----------------------------------------------------------------------------
solve.equa = function (k) 
{

  f = function(u, n)
  {
    (1 + u^2/(n-1))^(-n/2)
  }
  
  cn = function (n, a) 
  {
    sqrt(a^2 * n / (n + 1 - a^2))
  }
  
  xn = function (n, a) 
  {
    
    ff = function (u) 
    {
      f(u, n)
    }
    
    c = cn(n - 1, a)
    
    2/sqrt(pi*(n-1))*exp(lgamma(n/2)-lgamma((n-1)/2)) * integrate(ff, lower = 0, upper = c)$value
  }
  
  equa = function (a) 
  {
    xn(k,a)-xn(k+1,a)
  }
  
  e = 1e-6
  s=seq(0,sqrt(k),e)
  w=length(s)
  r=s[w]
  while(equa(e)*equa(r)>0)
  {
    w=w-1
    r=s[w]
  }
  uniroot(equa,lower=e,upper=r)$root
}
k=c(4,25,100,500,1000)
for (i in k) {
  print(solve.equa(i))
}

## -----------------------------------------------------------------------------
set.seed(123)
y=c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
n=length(y)
N=1e4
lambda.t=1 
e=1e-3
lambda.tplus=1


for(t in 1:N)
{
  #EM
  lambda.tplus=1/(mean(y)+mean(pexp(y,1)/lambda.t))
  
  
  #更新&判断
  if(abs(lambda.tplus-lambda.t)<e) break
  
  lambda.t=lambda.tplus
}
print(lambda.tplus)

## -----------------------------------------------------------------------------
y=c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
length(y)/sum(y)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
unlist(lapply(trims, function(trim) mean(x, trim = trim)))
unlist(lapply(trims, mean, x = x))

## -----------------------------------------------------------------------------
data("mtcars")

formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

rsq <- function(mod) summary(mod)$r.squared

unlist(lapply(formulas, function(formulae) rsq(lm(formulae,mtcars))))


## -----------------------------------------------------------------------------
data("mtcars")

formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

rsq <- function(mod) summary(mod)$r.squared

matrix(unlist(lapply(1:10, function(i) 
{
rows <- sample(1:nrow(mtcars), rep = TRUE)
lapply(formulas, function(formulae) rsq(lm(formulae,mtcars[rows,])))
})),nrow=10)

## -----------------------------------------------------------------------------
data("mtcars")

vapply(mtcars, sd, numeric(1))

## -----------------------------------------------------------------------------
sd1=function(dat)
{
  vapply(dat[vapply(dat, is.numeric, logical(1))], sd, numeric(1))#先用一个vapply判断是不是数值列，再用一个vapply进行方差计算
}
sd1(mtcars)

## -----------------------------------------------------------------------------
mcsapply=function(X,FUN)
{
  result=parallel::mclapply(X,FUN)
  unlist(result)
}

mcsapply(mtcars,sd)
parallel::mclapply(mtcars, sd)

## -----------------------------------------------------------------------------
library(Rcpp)
library(installr)
N=1e3;thin=1e2
cppFunction('#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x=0, y=0;
  int n=20;
  double a=2,b=1;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0)=x;
    mat(i, 1)=y;
  }
  return(mat);
}')
dataC=gibbsC(N,thin)
print(head(dataC))

## -----------------------------------------------------------------------------
gibbsR=function(N, thin) {
  mat=matrix(nrow = N, ncol = 2)
  x=y=0
  n=20
  a=2;b=1
  for (i in 1:N) {
    for (j in 1:thin) {
     x=rbinom(1, n, y)
     y=rbeta(1, x+a, n-x+b)
    }
    mat[i, ]=c(x, y)
  }
  mat
}
dataR=gibbsR(N,thin)

print(head(dataR))

## -----------------------------------------------------------------------------
data=dataR-dataC

par(mfrow=c(1,2))
x=data[,1]
q=qnorm(rank(x)/length(x))
plot(q,x)

y=data[,2]
q=qnorm(rank(y)/length(y))
plot(q,y)

## -----------------------------------------------------------------------------
library(microbenchmark)
ts=microbenchmark(gibbsR=gibbsR(N,thin),gibbsC=gibbsC(N,thin))
summary(ts)[,c(1,3,5,6)]

