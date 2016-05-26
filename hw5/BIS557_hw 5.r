####2.18
##problem 1
##define y and M*g
f=function(x){
  exp(-x^2/2)*((sin(6*x))^2+3*(cos(x))^2*(sin(4*x))^2+1)
}
g=function(x){
  exp(-x^2/2)/sqrt(2*pi)
}
x=runif(10000,-5,5)
plot(x,f(x),type="h",main="Density Plot",ylim=c(0,5))
points(x,g(x),col="red",cex=0.1)
maxf=optimize(f,c(-2,2),maximum=T)
maxg=optimize(g,c(-2,2),maximum=T)
M=round(maxf$objective/maxg$objective)+2
points(x,M*g(x),col="green",cex=0.1)
##apply accept-reject algorithm
rv=NULL
j=0
for(i in 1:2500){
repeat{
u=runif(1)
y=rnorm(1)
j=j+1
if((u*M*g(y))<=f(y))
 break
}
  rv=c(rv,y)}
##Deduce normalizing constant from acceptance rate
pr=2500/j
c=1/(pr*M)
hist(rv,breaks=40)
plot(x,c*f(x),main="Normalized f(x)",cex=0.2)

####problem 2
##Generate data u
seed = 1.0
myrnd = function() {
  seed <<- ((2^16 + 3) * seed) %% (2^31)
  return(seed/(2^31))
}
u = sapply(1:3000, function(i)myrnd())
##plot u by order
ru=sort(u)
plot(ru,main="Sorted U",cex=0.1)
##histgram of u
hist(u)
##perform a ks test
ks.test(u,"punif")
##make a 3D scatterplot
library("rgl")
y = matrix(u, nrow=1000, ncol=3, byrow=TRUE)
plot3d(y[,1], y[,2], y[,3], axes=TRUE, xlab="", ylab="", zlab="")

####Problem 3
##poisson
expo=function(lam){
  -1/lam*log(runif(1))
}
a=4
rv=NULL
for (j in 1:1000){
i=0
sum=0
while(sum<=1){
  sum=sum+expo(a)
  i=i+1
}
rv=c(rv,i-1)
}
hist(rv,freq=F,breaks=seq(-1,15,1))

ppois=function(x){
  a^x*exp(-a)/factorial(x)
}
x=0:15
lines(x-1,ppois(x))

##weibull
w=function(x,alpha,beta){
  beta*(-log(1-x))^(1/alpha)
}
u=runif(1000)
x=w(u,1,1)
hist(x,freq=F,breaks=seq(0,8,1))
alpha=1
beta=1
wei=function(x){
  alpha/beta*(x/beta)^(alpha-1)*exp(-((x/beta)^alpha))
}
x=seq(0,8,by=0.2)
lines(x,wei(x))

##triangle
t=function(x){
  (sqrt(x/2))^(x<=0.5)*(1-0.5*sqrt(2-2*x))^(1-(x<=0.5))
  }
u=runif(1000)
x=t(u)
hist(x,freq=F,ylim=range(0,2))

tri=function(x){
  4*x^(x<=0.5)*(1-x)^(x>0.5)
}
x=seq(0,1,0.5)
lines(x,tri(x))
