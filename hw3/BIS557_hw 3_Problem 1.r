#######Problem 1: Seasonal Influenza
###plot lamt
t=seq(1,365)
lamt=103*(1+0.8*cos(2*pi*t/365-3*pi/10))
plot(t,lamt,type="l",col="red",lwd=3)
###simulate X1-X365 and plot it over the curve of lamt
x=rep(0,365)
for(i in 1:365)
{
  x[i]=rpois(1,lamt[i]) 
}
points(t,x,pch=20,cex=0.5) 

###calculate loglikelihood on(a,phi) in [-1,1]*[0,4pi] 
a=seq(-1,1,length.out=100)
phi=seq(0,4*pi,length.out=100)
lglk=matrix(0,nrow=length(a),ncol=length(phi))
for(j in 1:length(a))
{
  for(k in 1:length(phi))
  {
    lamt=103*(1+a[j]*cos(2*pi*t/365-phi[k]))
    lglk[j,k]=sum(x*log(lamt)-lfactorial(x)-103*lamt)
  }
}

###plot contour
contour(x=a,y=phi,lglk,xlab = "a", ylab = "phi",main="Contour plot for loglikelihood at lam=103")
###iteration method finally works!!
###since sum(log(factorial(x))) is a constant,we don't need it in calculation
dlam=function(a,lam,phi){sum(x/lam-1-a*cos(2*pi*t/365-phi))}
ddlam=function(a,lam,phi){sum(-x/(lam^2))}
da=function(a,lam,phi){sum(x*cos(2*pi*t/365-phi)/(1+a*cos(2*pi*t/365-phi))-lam*cos(2*pi*t/365-phi))}
dda=function(a,lam,phi){sum(-x*(cos(2*pi*t/365-phi))^2/(1+a*cos(2*pi*t/365-phi))^2)}
dphi=function(a,lam,phi){sum(a*x*sin(2*pi*t/365-phi)/(1+a*cos(2*pi*t/365-phi))-a*lam*sin(2*pi*t/365-phi))}
ddphi=function(a,lam,phi){sum(-x*a*(cos(2*pi*t/365-phi)+a)/((1+a*cos(2*pi*t/365-phi))^2)+a*lam*cos(2*pi*t/365-phi))}
lam=10
a=0.5
phi=1
esp=10^-6
ii=1
rlam=NULL
ra=NULL
rphi=NULL
repeat{
  lamnew=lam-dlam(a,lam,phi)/ddlam(a,lam,phi)
  anew=a-da(a,lam,phi)/dda(a,lam,phi)
  phinew=phi-dphi(a,lam,phi)/ddphi(a,lam,phi)
  if(((lam-lamnew)^2+(a-anew)^2+(phi-phinew)^2)<esp){
    break
  }
  rlam[ii]=lamnew
  ra[ii]=anew
  rphi[ii]=phinew
  ii=ii+1
  lam=lamnew
  a=anew
  phi=phinew
}
###Plot your iterates on the heatmap/contour plot
points(ra,rphi,col="red",cex=0.5)
points(ra,rphi,col="blue",type="l")
##another method: gradient super unstable! Still doesn't work!
dlam=function(a,lam,phi){sum(x/lam-1-a*cos(2*pi*t/365-phi))}
da=function(a,lam,phi){sum(x*cos(2*pi*t/365-phi)/(1+a*cos(2*pi*t/365-phi))-lam*cos(2*pi*t/365-phi))}
dphi=function(a,lam,phi){sum(a*x*sin(2*pi*t/365-phi)/(1+a*cos(2*pi*t/365-phi))-a*lam*sin(2*pi*t/365-phi))}
lam=10
a=0.5
phi=1
esp=10^-6
ii=1
rlam=NULL
ra=NULL
rphi=NULL
repeat{
  lamnew=lam+0.1*dlam(a,lam,phi)
  anew=a+0.1*da(a,lam,phi)
  phinew=phi+0.1*dphi(a,lam,phi)
  if(((lam-lamnew)^2+(a-anew)^2+(phi-phinew)^2)<esp){
    break
  }
  rlam[ii]=lamnew
  ra[ii]=anew
  rphi[ii]=phinew
  ii=ii+1
  lam=lamnew
  a=anew
  phi=phinew
}

