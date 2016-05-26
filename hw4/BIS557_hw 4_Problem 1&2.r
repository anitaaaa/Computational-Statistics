########Probem 1: Weibull regression for a clinical trial (HW4)
####install package "survival"
install.packages("survival")
library(survival)
head(ovarian)
####setting covariate and observed values from dataset "ovarian"
attach(ovarian)
X=cbind(1,age,resid.ds,rx,ecog.ps)
t=futime
w=1-fustat ##Note that the definition of Wi in dataset is different from our definition
detach(ovarian)
####Newton Raphson Method for optimizing log likelihood and find alpha, beta
beta=1
alpha=c(1,1,1,1,1)
p=c(alpha,beta)
esp=10^-8
repeat{
  m=0
  n=0
  for (i in 1:length(t)){
    m1=-(c(X[i,],log(t[i]))%*%t(c(X[i,],log(t[i])))*(t[i]^beta)*as.numeric(exp(t(X[i,])%*%alpha))
         +(1-w[i])*matrix(c(rep(0,6*6-1),beta^(-2)),nrow=6,ncol=6))
    m=m+m1 ##hessian
    n1=c((1-w[i])*X[i,]-exp(t(X[i,])%*%alpha)*X[i,]*(t[i]^beta),
         (1-w[i])*(log(t[i])+beta^(-1))-exp(t(X[i,])%*%alpha)*(t[i]^beta)*log(t[i]))
    n=n+n1 ##gradient
  }
  inm=solve(m)
  pnew=p-inm%*%n
  if (as.numeric((t(pnew-p)%*%(pnew-p)))<esp){
    break
  }
  p=pnew
  alpha=p[1:5]
  beta=p[6]
}
####calculate the standard error of parameters using hessian matrix
variance=solve(-m) ##covariance matrix is the inverse of negative poisson
se=sqrt(diag(variance)) ##SE is square root of diagnal variance
####print estimated parameters and their SE
estimate=cbind(pnew,se)
####fit the model using survreg in R
fit = survreg(Surv(futime, fustat) ~ age + resid.ds + rx + ecog.ps,
              data=ovarian, dist='weibull')
summary(fit)
?survreg


########Problem 2: Poisson regression
####generate dataset x and z
n = 100
d = 13
z = array(NA, dim=c(n,d))
for(i in 1:n) z[i,] = c(1, runif(d-1, -3,1))
beta = rnorm(d,mean=0,sd=1)
x = rpois(n, exp(z %*% beta))
####Newton Raphson Method for optimizing log likelihood, define function poisreg
poisreg=function(x,z){
  p=rep(0.1,d)
  esp=10^-8
  repeat{
    h=0
    g=0
    for (i in 1:n){
      h1=-as.numeric(exp(t(z[i,])%*%p))*(z[i,]%*%t(z[i,]))
      h=h+h1 ##hessian
      g1=(x[i]-as.numeric(exp(t(z[i,])%*%p)))*z[i,]
      g=g+g1 ##gradient
    }
    inh=solve(h)
    pnew=p-inh%*%g
    if (as.numeric((t(pnew-p)%*%(pnew-p)))<esp){
      break
    }
    p=pnew
  }
  variance=solve(-h)
  se=sqrt(diag(variance)) 
  estimate=cbind(pnew,se)
  return(list(estimate,variance))
}
poisreg(x,z)
####glm 
z1=z[,-1]
result=glm(x ~ z1, family="poisson")
summary(result)

