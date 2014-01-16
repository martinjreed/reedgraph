ddpln <- function(x,a=2.5,b=.01,t=0.45,v=6.5){
  # Density of Double Pareto LogNormal distribution
  c <- (a * b /(a+b))
  


  # z is dominant for larger x values
  norm1<-pnorm((log(x)-v-(a*t^2))/t)
  expo1<- a*v+(a^2*t^2/2)
  z<- (x^(-a-1)) * exp(expo1)*(norm1)
  # y is dominant for small x values
  norm2<-pnorm((log(x)-v+(b*t^2))/t)
  expo2<- -b*v+(b^2*t^2/2) 
  y<- (x^(b-1)) * exp(expo2)*(1-norm2) # 1-norm is the complementary CDF of N(0,1)
  
  c*(z+y)
}
#x<-10^seq(0,5,by=0.1)# from 1 to 10^5
#plot(x,ddpln(x,a=2.5,b=.01,t=0.45,v=6.5),log='xy',type='l')
#plot(x,ddpln(x,a=2.8,b=.01,v=3.8,t=0.35),log='xy',type='l')

# here we plot DPLN generated random numbers using equ 10 from Reed vs Exponentials numbers
rdpln2<-function(n,a=2.5,b=0.01,t=0.45,v=6.5)
  { 
    library(VGAM)
    y<-rlnorm(n,meanlog=v,sdlog=t)*ppareto(n,location=1,shape=b)/rpareto(n,location=1,shape=a)
    e <- rexp(1e6,1/120)
    plot(density(y)); lines(density(e),col="red")
  }


rdpln <- function(n,a=2.5,b=0.01,t=0.45,v=6.5) {
    return (exp(v + t * rnorm(n) + rexp(n)/a - rexp(n)/b) )
}

plotDPLN <- function(N=10000)  {
    x <- rdpln(N)
    breaks <- c(0,10^seq(-1,4,0.1),max(x))
    
    binned <- hist(x,plot=FALSE,breaks=breaks)
    #binned <- hist(x,plot=FALSE)
    plot(binned$mids,binned$density,log="xy",xlim=c(0.1,1e4))
    return(binned)
}
