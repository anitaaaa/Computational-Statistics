####Estimator 1-3
c=1.3
rv1=(rnorm(10000)>c)+0
p1=mean(rv1)
rv2=runif(10000,0,c)
p2=0.5-sum(exp(-rv2^2/2))*c/sqrt(2*pi)/10000
rv3=rexp(10000,1)
p3=sum(exp(-0.5*(rv3+c)^2+rv3))/sqrt(2*pi)/10000
###Calculate and plot Variance on the interval 0-5
c=seq(0,5,by=0.05)
var1=(1-pnorm(c))*(1-(1-pnorm(c)))/10000
var2=((c/(sqrt(pi)*2))*(pnorm(sqrt(2)*c)-0.5)-(0.5-(1-pnorm(c)))^2)/10000
var3=pnorm(1/sqrt(2)-sqrt(2)*c)*exp(1/4-c)/(2*sqrt(pi)*10000)-(1-pnorm(c))^2/10000
plot(c,var1,type="l",col="red",main="Variance of three estimators on the range of c from 0 to 5",xlab="c",ylab="Variance",ylim=range(c(var1,var2,var3)))
legend(4,4e-5,c("Estimator1","Estimator2","Estimator3"),lty=c(1,1,1),lwd=c(1,1,1),col=c("red","green","blue"),cex=0.5)
lines(c,var2,col="green")
lines(c,var3,col="blue")
###comparison of estimator variances
f12=function(c){
  (1-pnorm(c))*(1-(1-pnorm(c)))/10000-(((c/(sqrt(pi)*2))*(pnorm(sqrt(2)*c)-0.5)-(0.5-(1-pnorm(c)))^2)/10000)
}
f13=function(c){
  (1-pnorm(c))*(1-(1-pnorm(c)))/10000-(pnorm(1/sqrt(2)-sqrt(2)*c)*exp(1/4-c)/(2*sqrt(pi)*10000)-(1-pnorm(c))^2/10000)
}
f23=function(c){
  ((c/(sqrt(pi)*2))*(pnorm(sqrt(2)*c)-0.5)-(0.5-(1-pnorm(c)))^2)/10000-(pnorm(1/sqrt(2)-sqrt(2)*c)*exp(1/4-c)/(2*sqrt(pi)*10000)-(1-pnorm(c))^2/10000)
}
x1=uniroot(f23,c(0,5))$root
x2=uniroot(f12,c(0,5))$root
x3=uniroot(f13,c(0,5))$root

####alternative neater
p=1-pnorm(c)
var1=p*(1-p)/10000
var2=((c/(sqrt(pi)*2))*(pnorm(sqrt(2)*c)-0.5)-(0.5-p)^2)/10000
var3=pnorm(1/sqrt(2)-sqrt(2)*c)*exp(1/4-c)/(2*sqrt(pi)*10000)-p^2/10000