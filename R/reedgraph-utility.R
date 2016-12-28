

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
