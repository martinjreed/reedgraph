###
###    Copyright (C) 2012 Martin J Reed              martin@reednet.org.uk
###    University of Essex, Colchester, UK
###
###    This program is free software; you can redistribute it and/or modify
###    it under the terms of the GNU General Public License as published by
###    the Free Software Foundation; either version 2 of the License, or
###    (at your option) any later version.
###
###    This program is distributed in the hope that it will be useful,
###    but WITHOUT ANY WARRANTY; without even the implied warranty of
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###    GNU General Public License for more details.
###
###    You should have received a copy of the GNU General Public License along
###    with this program; if not, write to the Free Software Foundation, Inc.,
###    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

require("bitops")

BLOOMLENGTH <- 256
BLOOMCHUNKS <- ceiling(BLOOMLENGTH/32)
ALLZEROS <- rep(0,BLOOMCHUNKS)

## CAREFUL this library uses a fragile encoding of length.
## you must not change the Bloom length and then use variables
## that were generated using the old Bloom length
rg.set.bloomlength <- function(m) {
  BLOOMLENGTH <<- m
  BLOOMCHUNKS <<- ceiling(BLOOMLENGTH/32)
  ALLZEROS <<- rep(0,BLOOMCHUNKS)
  
}

rg.prob.false.free <- function(m,alpha,beta) {
  return((1-0.618503^(m/alpha))^beta)
}

rg.exp.bloom.length <- function(alpha,beta) {
  notconverged <- TRUE
  Em <- (1-0.618503^(1/alpha))^beta
  m <- 2
  prod=1
  while(prod > 0) {
    prod=1
    for(i in 1:m-1) {
      prod <- prod * (1- (1-0.618503^(i/alpha))^beta)
    }
    Em <- Em + m * (1-0.618503^(m/alpha))^beta * prod
    m <- m + 1
  }
  cat("Expected is ",Em," converged at m=",m,"\n")
  return(Em)
}

rg.test.exp.bloom.length <- function(m,alpha,beta) {
  rg.set.bloomlength(m)
  k <- rg.bloom.optimum.k(m,alpha)
  rg.test.fp.averages()
  
}

## generate a random FID with n elements
rg.random.fid <-function(n,k=1) {
  fid <- ALLZEROS
  for(i in 1:n) {
    fid <- bitOr(fid,rg.createLid(k))
    #print(rg.bloom.to.binary(fid))
  }
  return(fid)
}

## count the number of false positives for
## an FID against t off path (random) LIDs
rg.test.fp <- function(n,t,k=1) {
  count <- 0
  fid <- rg.random.fid(n,k)
  #print("fid=")
  #print(rg.bloom.to.binary(fid))
  #print("lids")
  for(i in 1:t) {
    lid <- rg.createLid(k)
    #print(rg.bloom.to.binary(lid))
    if(rg.test.member(fid,lid))
      count <- count +1
  }
  return(count)
}

## find the false potisive rate for
## n (random) elements coded in a Bloom filter
## and compared against t off path (random) LIDs
## work this out over a averages and give the
## mean false postive rate
rg.test.fp.averages <-function(n,t,k=1,a) {
  sum <- 0
  for(i in 1:a) {
    sum <- sum + rg.test.fp(n,t,k)
  }
  return(sum/(a*t))
}

### example
### rg.test.fp.averages(10,20,10,100)
rg.test.fp.range.k <- function(n,t,k=c(1),a) {
  FP <- c()
  for(kay in k) {
    FP <- c(FP,rg.test.fp.averages(n,t,kay,a))
  }
  results <- data.frame(k=k,FP=FP)
  return(results)
    
}

### Function to print binary pattern.
### ONLY WRITTEN AS A DEMO, ONLY WORKS FOR BLOOMLENGTH <=32
rg.bloom.to.binary <- function(val) {
  result <- c()
  for(i in 1:BLOOMLENGTH) {
    mybit <- bitAnd(val,1)
    result <- c(mybit,result)
    val <- bitShiftR(val,1)
  }
  return(result)
}
rg.binary.to.bloom <- function(string) {
  result <- 0
  for(i in 1:BLOOMLENGTH) {
      result <- result *2

    result <- bitOr(result,string[i])
  }
  return(result)
}

### combinations possible with
### k hash functions, m length
rg.bloom.comb <- function(k,m) {

  prod <- 1.0
  for(i in m:(m-k+1)) {
    print(prod)
    prod <- i * prod
  }
  return(prod)
}

### Bloom filter false postive rate as estimated according to Bose 2007
### input: k hash functions, m length, n items in filter
### this is a strict lower bound for k>=2
### returns vector c(lower,upper)
### lower is a strict lower bound originally published by Bloom
### as exact figure, but Bose shows it is a strict lower bound.
### upper is approximate upper bound by Bose
rg.bloom.false <- function(k,m,n) {
  ## original Bloom figure is lower bound
  lbound <- (1-(1-1/m)^(k*n))^k

  ## from Bose
  p <- 1-(1-1/m)^(k*n)
  ubound <- lbound * (1 + k/p * sqrt((log(m) - k * log(p))/m ))


  return(c(lbound,ubound))
}



### optimium k for Bloom filter
### Input:
### n is number of elements bloom filter,
### m is length of bloom filter.
### returns k
rg.bloom.optimum.k <- function(m,n) {
  return(log(2) * m /n)
}

rg.createLid <- function(k=7) {
  lid <- rep(0,BLOOMCHUNKS)
  ## produce number 0-255 (zero offset)
  pos <- floor(runif(k) * BLOOMLENGTH)
  ## index in chunks (unit offset)
  ind <- pos %/% 32 +1
  ## set one bit in 0 ... 31
  val <- 2^((pos) %% 32)
  ## add to existing 
  for(i in 1:k) {
    lid[ind[i]] <- bitOr(lid[ind[i]],val[i])
  }
  return(lid)
}

### Test if element is member of bloom filter
### input:
### fid - Bloom filter
### lid - element to test if member of fid
### this could be simply written in-line but this
### call does not take much longer than the messy line
### it includes.
### returns - TRUE or FALSE as expected
rg.test.member <- function(fid,lid) {
  identical(bitXor(bitAnd(fid,lid),lid) ,ALLZEROS)
}

## THis look wrong! Should it be zero or unit offset?
rg.assign.lids.to.graph <- function(g,k=7) {
  for(e in 1:(length(E(g))) ) {
    ## have to assign a list to an attribute
    ## cannot assign a vector directly
    E(g)[e]$lid <- list(rg.createLid(k=k))
  }
  return(g)
}


### create fid along path, path is a sequence of node numbers
### as returned from   path <- get.shortest.paths(g,from,to)[[1]]
### graph - igraph
### value - vector representing fid
rg.create.fid <- function(graph,path) {
  fid <- ALLZEROS
  for(lid in E(graph,path=path)$lid) {
    fid <- bitOr(fid,lid)
  }
  return(fid)
}

### Count false positives on a path for given fid
### input
### graph - igraph with lid edge attributes set
### path - vector of vertex numbers for path
### fid - Bloom filter made up from lids on path
### returns list
### fpcount - count of false postives
### tpcount - count of true positives
### tncount - count of true negatives
###   ther are no false negatives for LIPSIN
#### WARNING may be broken, igraph stores the lids as seperate objects not concatonated list
rg.bloom.false.positive.on.path <- function(graph,path,fid) {
  if(!is.simple(graph)) {
    cat("Error in rg.bloom.false.positive, not a simple graph\n")
    return(NULL)
  }
  previous <- -1
  current <- path[1]
  fpcount <- 0
  tpcount <- 0
  tncount <- 0
  for(nxt in path[2:length(path)]) {
      n <- neighbors(graph,current)
      if(previous==-1) {
          tests <- n[n!=nxt]
      } else {
          tests <- n[n!=nxt & n!=previous]
      }
      for(t in tests) {
      if(rg.test.member(fid,E(graph)[ current %->% t ]$lid[[1]])) {
          fpcount <- fpcount + 1
        ##cat("false positive\n")
      } else {
          tncount <- tncount + 1
      }
    }
    previous <- current
    current <- nxt
    tpcount <- tpcount + 1
  }
  count <- list()
  count$fpcount <- fpcount
  count$tpcount <- tpcount
  count$tncount <- tncount
  return(count)
}

### Returns the false positive count on each path of length 1 ... diameter
### also returns true positive count and true negative count
rg.test.false.positives.all.paths <- function(graph) {
    
    fpcount <- c(0)
    tpcount <- c(0)
    tncount <- c(0)
    for(i in V(graph)) {
        paths <- get.shortest.paths(graph,i)
        for(path in paths$vpath) {
            nelements <- length(path)-1
            if(length(path)>1){
                fid <- rg.create.fid(graph,path)
                
                pcount <- rg.bloom.false.positive.on.path(graph,path,fid)
                if(is.na(fpcount[nelements])) {fpcount[nelements] <- 0}
                if(is.na(tpcount[nelements])) {tpcount[nelements] <- 0}
        if(is.na(tncount[nelements])) {tncount[nelements] <- 0}
        fpcount[nelements] <- fpcount[nelements] + pcount$fpcount
        tpcount[nelements] <- tpcount[nelements] + pcount$tpcount
        tncount[nelements] <- tncount[nelements] + pcount$tncount
      }
    }
}
        result <- list()
    
  result$fpcount <- fpcount
  result$tpcount <- tpcount
  result$tncount <- tncount
  return(result)
}

rg.test.false.posititives.loop <- function(graph,k=7) {
  rate <- c(0)
  count <- c(0)

  for(i in 1:50) {
    graph <- rg.assign.lids.to.graph(graph,k=k)
    result <- rg.test.false.positives.all.paths(graph)
    ## rate or count not long enough make them long
    ## enough and fill missing values with zero
    len <- length(result$rate)
    if(is.na(rate[len])) {
      rate[len] <- 0
      rate[is.na(rate)] <- 0
      count[len] <- 0
      count[is.na(count)] <- 0
    }
    rate <- rate + result$rate
    count <- count + result$count
  }
  count[is.na(count)] <- 0
  rate[is.na(rate)] <- 0
  range <- 1:length(count)
  bounds <- mapply(rg.bloom.false,k,BLOOMLENGTH,range)
  ## calculate mean rate from all results.
  result <- data.frame(N=range,Count=count,Rate=rate/count,Lower=bounds[1,],
                       Upper=bounds[2,])
  return(result)
}
