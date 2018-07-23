log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}




zipf.mandelbrot <- function(N=1000,alpha=1.0,var=0){
  k = 1:N
  dist <- length(k)/((k+var)^alpha)
  nordist <- dist/sum(dist)
  return(nordist)
}

rzipfm <- function(n=1,N=1000,alpha=1.0,var=0.0) {

    tmp <- zipf.mandelbrot(N,alpha=alpha,var=0.0)
 
    sample(x = c(1:length(tmp)),size = n,replace = TRUE,prob = tmp)
}


kl <- function(p,q) {
    sum(p*log(p/q,2))
}

## Total variation distance
## for finite countable alphabet tvd(p,q) <= sqrt(0.5*kl(p,q))
tvd <- function(p,q) {
    0.5*sum(abs(p-q))
}


geomean <- function(x, na.rm = FALSE, trim = 0, ...) {
    exp(mean(log(x, ...), na.rm = na.rm, trim = trim, ...))
}

geosd <- function(x, na.rm = FALSE, ...) {
    exp(sd(log(x, ...), na.rm = na.rm, ...))
}

conf.int <- function(data,interval=0.95) {
    n <- length(data)
    m <- mean(data)
    se <- sd(data)/sqrt(n)
    
    ci <- qt(interval/2+0.5,(n-1)) * se
    return(list(mean=m,cimin=(m-ci),cimax=(m+ci)))
}

log.conf.int.cox <- function(data,interval=0.95) {
    ## see Confidence Intervals for the Mean of a Log-Normal Distribution
    ## Ulf Olsson
    ## Swedish University of Agricultural Sciences
    ## Journal of Statistics Education Volume 13, Number 1 (2005)
    ## url http://www.amstat.org/publications/jse/v13n1/olsson.html
    n <- length(data)
    m <- mean(log(data))
    s <- var(log(data))
    
    cil <- qt(interval/2+0.5,(n-1)) * sqrt(s/n + s^2/(2*n -2))
    cimax <- m + s/2 + cil
    cimin <- m + s/2 - cil
    est <- exp(m+s/2)
    cimax <- exp(cimax)
    cimin <- exp(cimin)
    return(list(mean=mean(data),est=est,cimin=cimin,cimax=cimax))
}

log.conf.int.bayes <- function(data) {
    ## see http://stats.stackexchange.com/questions/33382/how-do-i-calculate-a-confidence-interval-for-the-mean-of-a-log-normal-data-set
    library(bayesm)

    ## simulated data
    ##mu <- 0
    ##sdv <- 1
    ##y <- exp(rnorm(1000, mean=mu, sd=sdv))
    
    ## model matrix
    X <- model.matrix(log(data)~1)
    ## prior parameters
    Theta0 <- c(0)
    A0 <- 0.0001*diag(1)
    ## Jeffreys prior for the normal model; set nu0 to 1 for the lognormal model
    nu0 <- 0 
    sigam0sq <- 0
    ## number of simulations
    n.sims <- 5000

    ## run posterior simulations
    Data <- list(y=log(data),X=X)
    Prior <- list(betabar=Theta0, A=A0, nu=nu0, ssq=sigam0sq)
    Mcmc <- list(R=n.sims)
    ## supress ouptput of runireg
    sink("/dev/null")
    bayesian.reg <- runireg(Data, Prior, Mcmc)
    sink()
    mu.sims <- t(bayesian.reg$betadraw) # transpose of bayesian.reg$betadraw
    sigmasq.sims <- bayesian.reg$sigmasqdraw

    ## posterior simulations of the mean of y: exp(mu+sigma²/2)
    lmean.sims <- exp(mu.sims+sigmasq.sims/2)

    ## credibility interval about lmean:
    return(quantile(lmean.sims, probs = c(0.025, 0.975)))
}


summaryNECI <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                       c( Num    = length2(xx[,col], na.rm=na.rm),
                         mean = mean   (xx[,col], na.rm=na.rm),
                         sd   = sd     (xx[,col], na.rm=na.rm)
                         )
                   },
                   measurevar,
                   na.rm
                   )
    
                                        # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$Num)  # Calculate standard error of the mean
    
                                        # Confidence interval multiplier for standard error
                                        # Calculate t-statistic for confidence interval: 
                                        # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$Num-1)
    datac$ciPos <- datac[[measurevar]] + datac$se * ciMult
    datac$ciNeg <- datac[[measurevar]] - datac$se * ciMult
   
    return(datac)
}


summaryNEQuantile <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      quant.interval=.50, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                       c( Num    = length2(xx[,col], na.rm=na.rm),
                         mean = mean   (xx[,col], na.rm=na.rm),
                         sd   = sd     (xx[,col], na.rm=na.rm),
                         ciPos = quantile(xx[,col],
                             probs=c(0.5 + quant.interval/2),names=FALSE),
                         ciNeg = quantile(xx[,col],
                             probs=c(0.5 - quant.interval/2),names=FALSE)
                         )
                   },
                   measurevar,
                   na.rm
                   )
    
                                        # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))
    
    return(datac)
}

rg.lineplot <- function(results,summaryAlg=summaryNECI,eb=TRUE,
                        xscale="X",
                        yscale="Y",
                        keys=NULL) {
    reslong <- melt(results,variable.name="Type",
                    value.name="Entries",
                    measure.vars=colnames(results)[2:length(colnames(results))])
    names(reslong)[1] <- "x"
    ressum <- summaryAlg(reslong,"Entries",groupvars=c("Type","x"))
    if(is.null(keys))
        keys <- levels(reslong$Type)
    p <- NULL
    pd <- position_dodge(0)
    #pl <<- position_dodge(.1)
    pl <- ggplot(ressum,aes(x=x,y=Entries)) +
        geom_line(aes(fill=factor(Type)),position=pd) +
            geom_point(size=4,position=pd,aes(shape=factor(Type))) +
                scale_shape_discrete(labels=keys) +
                    scale_x_continuous(name=xscale) +
                        scale_y_continuous(name=yscale)+
                            theme(legend.position=c(.25, .8),
                                  legend.title=element_blank())
    if(eb)
        pl <- pl + geom_errorbar(aes(ymin=ciNeg, ymax=ciPos),
                                      width=1.5, position=pd)

    pl
}

## results - wide format results, first column must be x axis
rg.boxplot <- function(results,
                       yscale="Y",
                       xscale="X",
                       keytitle="",
                       keys=NULL) {
    #theme_set( theme_bw(base_size=24,base_family= "serif"))
    ## to plot results you can use something like
    reslong <- melt(results,variable.name="Type",
                    value.name="Entries",
                    measure.vars=colnames(results)[2:length(colnames(results))])
    ressum <- summarySE(reslong,"Entries",groupvars=c("Type",colnames(results)[1]))
    ###reslong$Type <- factor(reslong$Type,
###                       levels=c("BFwithTables","BFoneTable","L2switch"))
    update_geom_defaults("point", list(colour = NULL))
    if(is.null(keys))
        keys <- levels(reslong$Type)
                                        #print(factor(reslong[[1]]))

                                        #return(TRUE)
    names(reslong)[1] <- "x"
    p <- ggplot(reslong,aes(x=factor(x),y=Entries)) +
        geom_boxplot(aes(colour=factor(Type))) +
            ## see http://colorbrewer2.org
            scale_colour_manual(labels=keys,
                                values=c("#e66101","#5e3c99", "#fdb863", "#b2abd2")) +
            ##scale_colour_brewer(palette="PuOr") +
            ##scale_colour_grey(start=0,end=0.5,
            ##                  labels=c("A","B")) +
            ##scale_colour_discrete(name="Implementation") +
                theme(legend.position=c(.2, .85),legend.title=element_blank()) +
                    scale_x_discrete(name=xscale) +
                        ##scale_y_continuous(limits=c(0,100),name=
                        scale_y_continuous(name=yscale)
    update_geom_defaults("point", list(colour = "black"))
    p
}

rg.cdf <- function(results,
                   yscale="Y",
                   xscale="X",
                   keytitle="",
                   keys=NULL,
                   xmax=NULL,
                   xlimit=NULL,
                   rectangle=NULL) {
    theme_set(theme_grey(base_family="Times",base_size=18))
    #theme_set( theme_bw(base_size=24,base_family= "serif"))
    ## to plot results you can use something like
    reslong <- melt(results,variable.name="Type",
                    value.name="Entries",
                    measure.vars=colnames(results)[2:length(colnames(results))])
    ressum <- summarySE(reslong,"Entries",groupvars=c("Type",colnames(results)[1]))
###reslong$Type <- factor(reslong$Type,
###                       levels=c("BFwithTables","BFoneTable","L2switch"))
    update_geom_defaults("point", list(colour = NULL))
    if(is.null(keys))
        keys <- levels(reslong$Type)
    names(reslong)[1] <- "x"

    p <- ggplot(reslong,aes(x=Entries,colour=Type)) +
        stat_ecdf(size=1.2) +
            ## see http://colorbrewer2.org
            scale_colour_manual(labels=keys,
                                values=c("#e66101", "#fdb863", "#b2abd2","#5e3c99")) +
    ##scale_colour_brewer(palette="PuOr") +
    ##scale_colour_grey(start=0,end=0.5,
    ##                  labels=c("A","B")) +
    ##scale_colour_discrete(name="Implementation") +
    theme(legend.position=c(.85, .10),legend.title=element_blank()) +
        scale_x_continuous(name=xscale) +
            ##scale_y_continuous(limits=c(0,100),name=
            scale_y_continuous(name=yscale)
    if(!is.null(xlimit))
        p <- p+ coord_cartesian(xlim = xlimit)
    if(!is.null(rectangle))
        p <- p+annotate("rect",xmin=rectangle[1],xmax=rectangle[2],
                        ymin=rectangle[3],ymax=rectangle[4], alpha=0.05,
                        linetype=2,color="black",size=0.2)

    update_geom_defaults("point", list(colour = "black"))
    p
}
