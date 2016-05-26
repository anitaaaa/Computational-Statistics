######Problem 2 
##generate lambda form exp(alpha)
alpha=0.7
lam=rexp(1,alpha)
lam
##generate 13 x from poisson(lambda)
n=13
x=rpois(n,lam)
x
##map estimate for lambda
lam_map=sum(x)/(n+alpha)
##calculate f*g
m=10000
lam_fg_=rnorm(5*m,lam_map,sqrt((sum(x))/(n+alpha)))
lam_fg=(lam_fg_[lam_fg_>0])[1:m]
##calculate g
lamg=function(l){
  (n+alpha)/sqrt(2*pi*sum(x))*exp(-(n+alpha)^2*(l-sum(x)/(n+alpha))^2/(2*sum(x)))
}
lam_g=lamg(lam_fg)
##calculate p
lamp=function(l){
  product=1
  product1=1
  for(i in 1:n){
    product=l^x[i]/factorial(x[i])*exp(-l)
    product1=product*product1
  }
  product1*alpha*exp(-alpha*l)
}
lam_p=lamp(lam_fg)

##posterior mean and sd from importance sampling
e_imp=sum(lam_fg*lam_p/lam_g)/sum(lam_p/lam_g)
sd_imp=sqrt(sum(lam_fg^2*lam_p/lam_g)/sum(lam_p/lam_g)-e_imp^2)

##mle and asymptotic se from MLE
e_mle=mean(x)
se_mle=sqrt(sum(x))/n

##expectation from map
lam_map


######problem 3
Sigma=seq(1,10,0.5)
sumrho=rep(0,length(Sigma))
for (s in 1:length(Sigma)){
  sigma=Sigma[s]
  N=10000
  lam_1=1
  lam_mh=rep(1,N)
  sum_rho=0
  for (i in 2:N){
    repeat{
      lam_star=rnorm(1,lam_1,sigma)
      if(lam_star>0){
        break
      }
    }
    rho=rbinom(1,1,min(lamp(lam_star)/lamp(lam_1)*pnorm(lam_1/sigma)/pnorm(lam_star/sigma),1))
    if(rho==1){
      lam_mh[i]=lam_star
    }
    if(rho==0){
      lam_mh[i]=lam_1
    }
    sum_rho=sum_rho+rho
    lam_1=lam_mh[i]
  }
  sumrho[s]=sum_rho
}
plot(Sigma,sumrho,type="l")

#####
sigma=Sigma[1]
N=10000
lam_1=1
lam_mh=rep(1,N)
for (i in 2:N){
  repeat{
    lam_star=rnorm(1,lam_1,sigma)
    if(lam_star>0){
      break
    }
  }
  rho=rbinom(1,1,min(lamp(lam_star)/lamp(lam_1)*pnorm(lam_1/sigma)/pnorm(lam_star/sigma),1))
  if(rho==1){
    lam_mh[i]=lam_star
  }
  if(rho==0){
    lam_mh[i]=lam_1
  }
  lam_1=lam_mh[i]
}

e_mh=sum(lam_mh*lamp(lam_mh))/sum(lamp(lam_mh))
sd_mh=sum(lam_mh^2*lamp(lam_mh))/sum(lamp(lam_mh))-e_mh^2

hist(lam_mh,col=rgb(0,0,0,0.5),nclass=20,main="Overlapping Histogram",xlab="lambda")
hist(lam_fg,col=rgb(1,0,1,0.5),nclass=60,add=T)
box()


hist(lam_mh,col=rgb(0,0,0,0.5),nclass=20,main="Overlapping Histogram",xlab="lambda")
box()
abline(v=lam_map,col="green")
abline(v=e_mle,col="blue")
abline(v=e_imp,col="red")
