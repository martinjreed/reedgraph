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

BLOOMLENGTH <- 64
BLOOMCHUNKS <- ceiling(BLOOMLENGTH/32)
ALLZEROS <- rep(0,BLOOMCHUNKS)

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


rg.assign.lids.to.graph <- function(g,k=7) {
  for(e in 0:(length(E(g))-1) ) {
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
  for(lid in E(g,path=path)$lid) {
    fid <- bitOr(fid,lid)
  }
  return(fid)
}
### Count false postives on a path for given fid
### input
### graph - igraph with lid edge attributes set
### path - vector of vertex numbers for path
### fid - Bloom filter made up from lids on path
### returns list
### fpcount - count of false postives
### tpcount - count of true positives
### tncount - count of true negatives
###   ther are no false negatives for LIPSIN
rg.bloom.false.postive.on.path <- function(graph,path,fid) {
  if(!is.simple(graph)) {
    cat("Error in rg.bloom.false.postive, not a simple graph\n")
    return(NULL)
  }
  previous <- -1
  current <- path[1]
  fpcount <- 0
  tpcount <- 0
  tncount <- 0
  for(nxt in path[2:length(path)]) {
    n <- neighbors(graph,current)
    tests <- n[n!=nxt & n!=previous]
    for(t in tests) {
      if(rg.test.member(fid,E(graph)[ current %--% t ]$lid)) {
        fpcount <- fpcount + 1
        ##cat("false positive\n")
      }
      tncount <- tncount + 1
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

rg.test.false.postitives.all.paths <- function(graph) {
  
  rate <- c(0)
  count <- c(0)
  for(i in V(graph)) {
    paths <- get.shortest.paths(graph,i)
    for(path in paths) {
      nelements <- length(path)-1
      if(length(path)>1){
        fid <- rg.create.fid(graph,path)
        
        pcount <- rg.bloom.false.postive.on.path(graph,path,fid)
        if(is.na(rate[nelements])) {rate[nelements] <- 0}
        ## out of all of the non true positives what is the portion
        ## of false postitives?
        ## Sum this to running total so that we can create a mean later.
        rate[nelements] <- rate[nelements] +
          pcount$fpcount / (pcount$tncount + pcount$fpcount)
        ## this should not be:
        pcount$fpcount / (pcount$tpcount + pcount$tncount + pcount$fpcount)
        if(is.na(count[nelements])) {count[nelements] <- 0}
        ## the number of counts, saved for later to calculate mean
        count[nelements] <- count[nelements] + 1
      }
    }
    
  }
  result <- list()
  result$rate <- rate
  result$count <- count
  return(result)
}

rg.test.false.postitives.loop <- function(graph,k=7) {
  rate <- c(0)
  count <- c(0)

  for(i in 1:50) {
    graph <- rg.assign.lids.to.graph(graph,k=k)
    result <- rg.test.false.postitives.all.paths(graph)
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
