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

require("igraph",quietly=TRUE)
require("graph",quietly=TRUE,warn.conflicts=FALSE)
require("RBGL",quietly=TRUE,warn.conflicts=FALSE)
require("genalg",quietly=TRUE)
require("inlinedocs")


### This is rather dangerous and discouraged as it erases any state in
### the global environment that would mask objects from this
### package. DO NOT include in version submitted to CRAN
.onAttach <- function(libname, pkgname) {
  nummasked <- length(intersect(ls(".GlobalEnv"),ls("package:reedgraph")))
  if ( nummasked >0 ) {
    packageStartupMessage(
      paste("   WARNING ", nummasked, " objects from reedgraph\n",
            "            have overwritten the same objects in .GlobalEnv\n"))
  }
  rm(list=intersect(ls(".GlobalEnv"),ls("package:reedgraph")),
     pos=".GlobalEnv")
}

rg.update2package <- function # Update objects in .GlobalEnv
### Update objects in .GlobalEnv to the latest in the package
### Usefull if you have stale values in .Rdata files saved
### in a previous workspace.
() {
  print("Removing old definitions")
  rm(list=intersect(ls(".GlobalEnv"),ls("package:reedgraph")),
     pos=".GlobalEnv")
}

### return all the edges as a list of lists(head,tail)
### in graphNEL
### result is list of edges
### edgeMatrix might be better
rg.edgeL <- function(g) {
  ed <- list()
  for(i in nodes(g)) {
    for(el in graph::edges(g)[[i]]) {
      ed <- c(ed,list(c(i,el)))
    }
  }
  ed
}


### returns the attribute of the graphNEL as a double
### vector (in natural edge order)
### input
### g - graphNEL
### attr="weight" for the default attribute
### returns c(....)
rg.edgeVector <- function(g,attr="weight") {
  return(as.double(edgeData(g,attr=attr)))
}

rg.edgeVectorGamma <- function(g) {
  c <- as.double(edgeData(g,attr="capacity"))
  f <- as.double(edgeData(g,attr="weight"))
  gamma <- (c-f)/c
  return(gamma)
}

rg.demandsVector <- function(demands,attr="demand") {
  return(as.double(lapply(res$demands,"[[",attr)))
}

rg.count.paths <- function(demands,numlist=NULL) {
  num <- 0
  if(is.null(numlist))
    numlist <- seq(1,length(demands))
  for(i in numlist){
    for(p in names(demands[[i]]$paths)) {
      num <- num + 1
    }
  }
  return(num)
}


### Creates an adjacency matrix containing edge ids
### useful to get an edge number given two nodes
### g - graphNEL (with L edges)
### returns matrix (node/node) where element is edge
### number 1 ..... L between node/node

adjMatrix <- function(g) {
  em <- edgeMatrix(g)
  L <- length(em)/2
  N <- length(nodes(g))
  m <- matrix(nrow=N,ncol=N)
  for(i in seq(1:L)) {
    m[as.integer(em[1,i]),as.integer(em[2,i])] <- i
  }
  return(m)
}


### Calculates the number of edge disjoint paths
### between every pair of nodes
### Returns a matrix where element [u,v] is the
### number of paths between u and v where
### 0 means disconnected (or to itself)
### 1 for one path etc
### THIS IS INCORRECT
rg.num.edge.disjoint.paths <- function(g) {
  # set ip a matrix to count (will be symmetrical for undirected graph)
  count <- matrix(0,ncol=length(nodes(g)),
                      nrow=length(nodes(g)))

  
  for(i in nodes(g)) {
      splist <- dijkstra.sp(g,start=i)
      for(j in nodes(g)) {
      if ( i != j ) {
        if (is.infinite(splist$distances[[j]])) {
          count[as.integer(i),as.integer(j)] <- 0
        } else {
          gtmp <- g
          splisttmp <- splist
          while( is.finite(splisttmp$distances[[j]])) {
            count[as.integer(i),as.integer(j)] <-
              count[as.integer(i),as.integer(j)] + 1
            path <- extractPath(i,j,splisttmp$penult)

            # removing an edge takes longer, so using set.weight
            # instead
            #gtmp <- rg.remove.edges(gtmp.path)
            gtmp <- rg.set.weight.on.path(gtmp,path,Inf)
            splisttmp <- dijkstra.sp(gtmp,start=i)
          }
        }
      }
    }
  }
  return(count)
}

### Counts number of paths in demands
rg.count.paths <- function(demands) {
  count <- 0
  for(demand in demands) {
    count <- count + length(demand$paths)
  }
  return(count)
}

### Removes edges from a graph
### g is GraphNEL
### path is vector of the path vertex labels
### from source to tail
### returns GraphNEL with path removed
rg.remove.edges <- function(g,path) {
  fromlist <- path[1:{length(path)-1}]
  tolist <- path[2:{length(path)}]

  for(i in 1:length(fromlist)) {
    
    g <- removeEdge(as.character(fromlist[[i]]),
               as.character(tolist[[i]])
               ,g)
  }
  return(g)
}

### Simple way to set demands
### val the value to assign to ALL demands
### array of node pairs, as characters
rg.set.demand <- function(val=1.0,nodes) {
  if(length(nodes) %% 2 != 0) {
    cat("Error in rg.set.demands nodes must be even\n")
    return(NULL)
  }
  j <- 1
  count <- 1
  demand <- list()
  for(i in seq(1,length(nodes),2)) {
    j <- i+1
    demand[[as.character(count)]]$source <- as.character(nodes[i])
    demand[[as.character(count)]]$sink <- as.character(nodes[j])
    demand[[as.character(count)]]$demand <- val
    count <- count + 1
  }

  return(demand)

}

rg.gen.demands <- function # Generate some demands
### Generates a list of list("number",source=,sink=,demand=) where
### each demand is identified by a number in the order of the demand
### therefore probably not best as a mutuable structure unless care is taken
### that number and indexing is used correctly.
(g,         ##< graphNEL
 num=NULL,  ##< num number of demands to generate , it will be this many across a
            ##  uniformly random selection of node pairs (i,j) where i != j
            ##  If NULL then it will be between every pair
 val=1.0,   ##< value of demand. Accepts a vector and iterates through in
            ##  order and will continue back to beggining if necessary
            ##  use runif(num,min,max) if uniform random required
 nodes=NULL ##< list of nodes to generate demands between. If NULL it will
            ##  be every node in the graph
 ) {
  if( is.null(nodes) ) nodes <- nodes(g)
  comm <- list()
  n=1;
  k <- 1
  vallen <- length(val)
  if( is.null(num) ) { # all pairs uniform demands
    for(i in nodes) {
      ##cat("doing",i,"\n")
      for(j in nodes) {
        if ( i != j ) {
          comm[[as.character(n)]] <- list(source=i,sink=j,demand=val[k])
          k <- k + 1
          if(k > vallen) k <- 1
          n <- n+1
        }
      }
    }
  } else { # random selection of nodes
    numnodes <- length(nodes)
    for(i in seq(1:num)) {
      fromto <- floor(runif(2,min=1,max=numnodes+1))
      while(fromto[1] == fromto[2])
        fromto <- floor(runif(2,min=1,max=numnodes+1))
      comm[[as.character(n)]] <- list(
                                      source=
                                      nodes[fromto[1]],
                                      sink=
                                      nodes[fromto[2]]
                                      ,demand=val[k])
      k <- k + 1
      if(k > vallen) k <- 1
      n <- n+1
    }
  }
  comm
  ### a list of elements list("number",source=,sink=,demand=)
}

### Sets all demands to a constant value
### val - value to set
### returns new demands
rg.set.demands.demand <- function(demands,val) {
  for(n in names(demands)) {
    demands[[n]]$demand <- val
  }
  demands
}

rg.set.capacity <- function # Sets the capacity of the edges on the graph
(g,  ##< a graphNEL object
 val ##< scalar or vector length of edges
 ) {
  em <- edgeMatrix(g)
  if (match("capacity",names(edgeDataDefaults(g)),nomatch=1)) {
    edgeDataDefaults(g,"capacity") <- 1.0
  }
  nodes <- nodes(g)
  edgeData(g,
           from=nodes[em[1,]],
           to=nodes[em[2,]],
           attr="capacity") <- val
  g
  ### graphNEL with capacity attribute set
}

### Sets the weight of the edges on the graph
### val - scalar or vector length of edges
rg.set.weight <- function(g,val) {
  em <- edgeMatrix(g)
  if (match("weight",names(edgeDataDefaults(g)),nomatch=1)) {
    edgeDataDefaults(g,"weight") <- 1.0
  }
  edgeData(g,
           from=nodes(g)[em[1,]],
           to=nodes(g)[em[2,]],
           attr="weight") <- val
  g
}

### make sure the labels are the same as vertex numbers
### this seems to be needed else there are problems with
### dijkstra.sp which takes in "labels" and returns vertices
### labels as numbers
### use g <- rg.relabel(g)
rg.relabel <- function(myg) {
  nodes(myg) <- as.character(seq(1:length(nodes(myg))))
  myg
}

### adds a value val to the weight attribute on each edge
### in path
### input graphNEL, path vector of integer, val is  the amount to add
### this is very inefficient if called repeatedly, better to copy
### and paste "Inline" if done many times
rg.addto.weight.on.path <- function(g,path,val,attr="weight") {
  nodes <- nodes(g)
    if(length(path) == 1) {
      edgeData(g,nodes[path[1]],nodes[path[2]],attr) <-
      edgeData(g,nodes[path[1]],nodes[path[2]],attr) +val
  } else {
    
    fromlist <- path[1:{length(path)-1}]
    tolist <- path[2:{length(path)}]
    newvals <- as.double(edgeData(g,nodes[fromlist],nodes[tolist],attr)) + val
    edgeData(g,nodes[fromlist],nodes[tolist],attr) <- newvals
  }
  return(g)
}


### sets the weight attribute on each edge
### in path to val WARNING THIS WILL OVERIDE OTHER PATHS THAT HAVE SET
### THE EDGE WEIGHT : you probably want rg.addto.weight.on.path()
### input graphNEL, path vector of integer, val is  the amount to set
###
rg.set.weight.on.path <- function(g,path,val) {
  if(length(path) == 1) {
    edgeData(g,as.character(path[1]),as.character(path[2]),"weight") <- val
  } else {
    fromlist <- path[1:{length(path)-1}]
    tolist <- path[2:{length(path)}]
    edgeData(g,as.character(fromlist),as.character(tolist),"weight") <- val
  }
  return(g)
}

### demands: 
##' Compute the multicommodity flow with shortest path
##'
##' Compute the multicommodity flow in a graph using naive shortest
##' path g graph of type graphNEL with the edge weight variable used
##' for the shortest path and capacities ignored (ie it may overbook
##' an edge)
##' @title Shortest Path Max Concurrent Flow
##' @param g graphNEL object with capacity attribute on each edge
##' @param demands a list each element is a list [source, sink, demand]
##' @return a list(demands,gflow,lambda,gamma) where demands is a
##' list(demand,source,sink,paths,flow) where paths is a list("X|Y|Z")
##' of the path (only one path in this function, others have
##' more). lambda is the max concurrent flow primal value (very
##' unoptimal in this case), gamma is the prmal value for the min
##' congestion flow (again highly suboptimal)
##' @author Martin Reed
rg.sp.max.concurrent.flow <- function # Compute the multicommodity flow with shortest path
### Compute the multicommodity flow in a graph using naive shortest
### path g graph of type graphNEL with the edge weight variable used
### for the shortest path and capacities ignored (ie it may overbook
### an edge)
(g,      ##< graphNEL object with capacity attribute on each edge
 demands ##< a list each element is a list [source, sink, demand]
 ) {

  # convert demands node labels to the index
  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  for(c in names(demands)) {
    demands[[c]]$flow <- 0
  }
  gsol <- g
  g.sp <- list()

  g <- rg.set.weight(g,1.0)
  gsol <- rg.set.all.graph.edge.weights(gsol)
  
  ccount <- 1
  for(i in demands) {
    ## calculate dijkstra.sp if we do not have one for this vertex
    if(is.null(g.sp[[i$source]]))
      g.sp[[i$source]] <- dijkstra.sp(g,start=i$source)$penult
    path <- extractPath(i$source,i$sink,g.sp[[i$source]])
    gsol <- rg.addto.weight.on.path(gsol,path,i$demand)
    demands[[ccount]] <- updateExplicitFlow(g,g.sp[[i$source]],i$demand,i)
    demands[[ccount]]$flow <- demands[[ccount]]$flow +i$demand
    ccount <- ccount + 1
  }
  
  f <- as.double(edgeData(gsol,attr="weight"))
  c <- as.double(edgeData(gsol,attr="capacity"))
  lambda <- min(c/f)

  gamma <- min((c-f)/c)
    

    
  #gsol <- rg.max.concurrent.flow.graph(gsol,demands)

  ##put back the demands labels instead of indices
  nodes(gsol) <- nodelabels
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  retval <-  list(demands=demands,gflow=gsol,lambda=lambda,gamma=gamma)

    
  return(retval)
### A list(demands,gflow,lambda,gamma) where demands is a
### list(demand,source,sink,paths,flow) where paths is a list("X|Y|Z")
### of the path (only one path in this function, others have
### more). lambda is the max concurrent flow primal value (very
### unoptimal in this case), gamma is the prmal value for the min
### congestion flow (again highly suboptimal)
}

### Compute the multicommodity flow in a graph using naive shortest
### path with first fit for g, graph of type graphNEL, with the edge
### weight variable used for the shortest path. If the capacity is
### reached it will set the edge weight to infinity so that a fully
### booked edge is omitted next time demands: a list each element is a
### list [source, sink, demand] WARNING NOT WRITTEN YET
rg.sp.ff.max.concurrent.flow <- function(g,demands) {
  # convert demands node labels to the index
  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  for(c in names(demands)) {
    demands[[c]]$flow <- 0
    demands[[c]]$paths <- list()
  }
  gsol <- g
  g.sp <- list()

  g <- rg.set.weight(g,1.0)
  gsol <- rg.set.all.graph.edge.weights(gsol)
  
  edgesfrom <- as.character(
    sapply(names(edgeData(g,attr="capacity")),function(x) strsplit(x,"\\|")[[1]])[1,])
  edgesto <- as.character(
    sapply(names(edgeData(g,attr="capacity")),function(x) strsplit(x,"\\|")[[1]])[2,])

  for(d in names(demands)) {
    i <- demands[[d]]
    gtmp <- g
    ## find edges in gsol that can not accept demand and set them to infinity
    weights <- as.double(edgeData(gsol,att="weight"))
    capacities <- as.double(edgeData(gsol,att="capacity"))
    
    free <- (capacities-weights)

    gtmpw <- as.double(edgeData(gtmp,att="weight"))
    gtmpw[i$demand > free] <- Inf
    edgeData(gtmp,from = edgesfrom, to = edgesto, attr="weight") <- gtmpw
    dij <- dijkstra.sp(gtmp,start=i$source)
    if(dij$distance[i$sink] < Inf) {
      path <- extractPath(i$source,i$sink,dij$penult)
      fromlist <- path[1:length(path)-1]
      tolist <- path[2:length(path)]

      weights <- as.double(edgeData(gsol,from=as.character(fromlist),
                                    to=as.character(tolist),att="weight"))
      capacities <- as.double(edgeData(gsol,from=as.character(fromlist),
                                       to=as.character(tolist),att="capacity"))
    
      minfree <- min(capacities-weights)

      free <- (capacities-weights)
      gsol <- rg.addto.weight.on.path(gsol,path,i$demand)
      demands[[d]] <- updateExplicitFlow(g,dij$penult,i$demand,i)
      demands[[d]]$flow <- i$demand
    }
   } 
  f <- as.double(edgeData(gsol,attr="weight"))
  c <- as.double(edgeData(gsol,attr="capacity"))
  lambda <- min(c/f)
  
  gamma <- min((c-f)/c)
  
  
  ##put back the demands labels instead of indices
  nodes(gsol) <- nodelabels
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  retval <-  list(demands=demands,gflow=gsol,lambda=lambda,gamma=gamma)
  return(retval)
}


rg.demands.relable.to.indices <- function # Relable-demands-to-indices
### Relable a demands list using indices from labels
(demands,    ##<< demands list in format shown elsewhere
 nodelabels  ##<< vector of nodelables c("1","2","fred"...)
 ) {
  demands <- lapply(demands,function(x) {x$source <- as.character(which(nodelabels==x$source));
                                         return(x)})
  demands <- lapply(demands,function(x) {x$sink <- as.character(which(nodelabels==x$sink));
                                         return(x)})

  ## Need to fix the path labels here as well!
  ## below will create the correct key, just need to update a key
  for(i in names(demands)) {
    demand <- demands[[i]]
    if(!is.null(demand$paths)) {
      demands[[i]]$paths <- list()
      for(p in names(demand$paths)) {
        tmp <- as.numeric(demand$paths[[p]])
        demands[[i]]$paths[[paste0(as.character(sapply(
          strsplit(p, "\\|")[[1]],function(x) which(nodelabels==x))),collapse="|")]] <- tmp
      }
    }
  }
 
  return(demands)
  ### The list of demands with node labels converted to indices
}

rg.demands.relable.from.indices <- function(demands,nodelabels) {
  demands <- lapply(demands,function(x) {x$source <- nodelabels[as.integer(x$source)]; return(x)})
  demands <- lapply(demands,function(x) {x$sink <- nodelabels[as.integer(x$sink)]; return(x)})
  # Need to fix the path labels here as well!
  # below will create the correct key, just need to update a key
  for(i in names(demands)) {
    demand <- demands[[i]]
    demands[[i]]$paths <- list()
    for(p in names(demand$paths)) {
      tmp <- as.numeric(demand$paths[[p]])
      demands[[i]]$paths[[paste0(nodelabels[as.integer(strsplit(p,"\\|")[[1]])],collapse="|")]] <- tmp
    }
  }
  return(demands)
}


### Compute the multicommodity flow in a graph using naive shortest
### path g graph of type graphNEL with the edge weight variable used
### for the shortest path and capacities ignored (ie it may overbook
### an edge)
### demands: a list each element is a list [source, sink, demand]
rg.sp.max.concurrent.flow.residual <- function(g,demands) {


  for(c in names(demands)) {
    demands[[c]]$flow <- 0
  }

  
  gsol <- g
  g.sp <- list()
  
  gsol <- rg.set.all.graph.edge.weights(gsol)
  
  ccount <- 1
  for(i in demands) {
    ## calculate dijkstra.sp if we do not have one for this vertex
    if(is.null(g.sp[[i$source]]))
      g.sp[[i$source]] <- dijkstra.sp(g,start=i$source)$penult
      path <- extractPath(i$source,i$sink,g.sp[[i$source]])
      gsol <- rg.addto.weight.on.path(gsol,path,i$demand)
      demands[[ccount]]$flow <- demands[[ccount]]$flow +i$demand
      ccount <- ccount + 1
  }
  
  f <- as.double(edgeData(gsol,attr="weight"))
  c <- as.double(edgeData(g,attr="capacity"))
  lambda <- min(c/f)

  gamma <- min((c-f)/c)
    
  #demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,lambda)
    
  #gsol <- rg.max.concurrent.flow.graph(gsol,demands)
  
  retval <-  list(demands=demands,gflow=gsol,lambda=lambda,gamma=gamma)
    
  return(retval)
}



### Utility function for MCF algorithm
### at the end this is the dual solution value
calcBetaRestricted <- function(demands,gdual) {
  Alpha <- 0
  
  for(demand in demands) {
    lcost <- Inf
    for(p in names(demand$paths)) {
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- pv[1:length(pv)-1]
      tolist <- pv[2:length(pv)]
      weights <- as.double(edgeData(gdual,from=fromlist,to=tolist,att="weight"))
      distance <- sum(weights)
      if(lcost > distance)
        lcost <- distance
    }
    Alpha <- Alpha + demand$demand * lcost
    
  }
  D <- sum(as.double(edgeData(gdual,attr="capacity"))*
           as.double(edgeData(gdual,attr="weight")))

  Beta <- D / Alpha
  Beta
}


### Utility function for MCF algorithm
### at the end this is the dual solution value
calcBeta <- function(demands,gdual) {
  Alpha <- 0
  
  for(demand in demands) {
    sp <- dijkstra.sp(gdual,demand$source)
    Alpha <- Alpha + demand$demand * sp$distances[[demand$sink]]
  }
  D <- sum(as.double(edgeData(gdual,attr="capacity"))*
           as.double(edgeData(gdual,attr="weight")))

  Beta <- D / Alpha
  Beta
}


### Utility function for MCF algorithm
### at the end this is the primal solution value
calcLambda <- function(demands) {
  d <- as.double(lapply(demands,"[[","demand"))
  f <- as.double(lapply(demands,"[[","flow"))
  lambda <- min(f/d)
  lambda
}

### Calculates the maximum concurrent flow using
### rg.fleischer.max.concurrent.flow It uses the demand scaling
### algorithm suggested by "Faster approximation schemes for
### fractional multicommodity flow problems" G. Karakostas,
### Proceedings ACM/SIAM SODA 2002
### g: graphNEL weights ignored
### e: approximation value for a w-opt solution
###    where (1+w) = (1-e)^2 (default e=0.1)
### progress: display graphical progress bar (default false)
rg.max.concurrent.flow.prescaled <- function(g,demands,e=0.1,progress=FALSE,ccode=TRUE,updateflow=TRUE,permutation="random",deltaf=1.0) {

  savedemands <- demands
  ## first obtain a reasonable feasible flow using
  ## a poor shortest path flow
  estimate <- rg.sp.max.concurrent.flow(g,demands)
  ## we know that estimate$lambda < Beta
  ## so scale demands by this much to make sure now
  ## that beta > 1
  estlambdascale <- estimate$lambda
  demands <- rg.rescale.demands(demands,estimate$lambda)

  ## Beta might be very high so obtain 2-approximate solution
  ## (1+w) = 2 = (1-e)^-3
  e2= 1- 2^(-1/3)

  if(progress != FALSE)
    progress <- "Calculating 2 opt"
  
  if (ccode) {
    res.2opt <- rg.fleischer.max.concurrent.flow.c(g,demands,e=e2,
                                                   updateflow=FALSE,
                                                   progress=progress,permutation,
                                                   deltaf)
  }  else {
    res.2opt <- rg.fleischer.max.concurrent.flow(g,demands,e=e2,
                                                 updateflow=FALSE,progress=progress)
  }
  
  ## now have 2-opt solution as beta_est so  beta < beta_est < 2*beta
  ## scaling by beta_est/2 give best demand for reduced running time
  scalef <- res.2opt$beta /2
  opt2scale <- scalef
  
  demands <- rg.rescale.demands(demands,scalef)


  if(progress != FALSE)
    progress <- "Main calculation"

  if (ccode) {
    res <- rg.fleischer.max.concurrent.flow.c(g,demands,e=e,
                                              progress=progress,updateflow=updateflow,
                                              permutation,
                                              deltaf)
  } else {
    res <- rg.fleischer.max.concurrent.flow(g,demands,e=e,
                                            progress=progress,updateflow=updateflow)
  }
  


  for(n in names(demands)) {
    res$demands[[n]]$demand <- savedemands[[n]]$demand
  }
  
  res$lambda <- calcLambda(res$demands)
    
  res$beta <- calcBeta(res$demands,res$gdual)
    
  res$estlambdascale <- estlambdascale
  res$opt2scale <- opt2scale
  res

}

### Compute the maximum concurrent flow See Fleischer "Approximating
### fractional multicommodity flow independent of the number of
### commodities", SIAM J. Discrete Maths, Vol 13/4, pp 505-520
###
### g: graphNEL (vertex labels need to be same as vertex index)
### demands: vector
###
### g: graphNEL with capacity data on each edge, edge
### weight is ignored (traffic is routed over any way to meet
### concurrent flow
###
### e: accuracy ( 0 < e < 1 ) 0 is more accurate, 1 less accurate
###
### permutation: how the demands are chosen can be either
###              c(....) integers specifying demand order
###              "fixed" done in fixed order
###              "random" done in random order
###              "lowest" done in lowest cost (lowest dual path) order
###
### output: list:
###              graph: with weight set to edge flow
###              demands: each with met demand and paths
rg.fleischer.max.concurrent.flow.c <- function(g,
                                                   demands,
                                                   e=0.1,
                                                   updateflow=TRUE,
                                                   progress=FALSE,
                                                   permutation="random",
                                                   deltaf=1.0) {

  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")
  retlist <- .Call("rg_fleischer_max_concurrent_flow_c_new",
                   as.double(cap),
                   as.integer(em-1),
                   as.integer(nv),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e)
                 )

  demands <- retlist$demands

  delta <- (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  #demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  beta <- calcBeta(demands,gdual)
  
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3

  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound
                   )
}


rg.max.concurrent.flow.int <- function # Integer minimum congestion concurrent Flow
### Computes the integer minimum congestion concurrent flow
### Based on testing for best integer flow so far as part of fractional
### multicommodity flow solution.
(g, ##< graphNEL
 demands, ##< demands - list (index character count) of list("source","sink",demand)
 e=0.1, ##< approximation limit 0 < e < 1.0, smaller is more accurate but takes longer
 updateflow=TRUE, ##< update the flow
 progress=FALSE, ##< display progress bar (not working in latest version)
 permutation="random", ##< only random in this version
 deltaf=1.0, ##< not used in this version
 prescaled=TRUE ## scaling linear demands for faster run-time
 ) {
    
    nodelabels <- nodes(g)
    demands <- rg.demands.relable.to.indices(demands,nodelabels)
    g <- rg.relabel(g)
    
    em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")
  retlist <- .Call("max_concurrent_flow_int",
                   as.double(cap),
                   as.integer(em-1),
                   as.integer(nv),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e),
                   prescaled
                 )

  demands <- retlist$demands
  delta <- (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  #demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  beta <- retlist$beta
  
  lambda <- retlist$lambda
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3

  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound,
                   gamma=retlist$gamma
                   )
### a list of many things:
### demands - list of lists("source","sink",val,flow,path) <-format needs checking
### gflow - graphNEL with flow on each edge as a weight
### gdual -graphNEL with dual graph lengths - only useful for analysis
### beta - the linear dual problem limit
### lambda - the linear primal problem limit
### phases - the number of phases it took
### ratio - the found ratio between primal/dual, would like this to be near 1.0
### bound - the expected ratio between primal/dual
### gamma - the fractional free capacity on most constrained link (bigger better).
}

### As above but using prescaling to speed it up
rg.max.concurrent.flow.prescaled.c <- function(g,
                                                   demands,
                                                   e=0.1,
                                                   updateflow=TRUE,
                                                   progress=FALSE,
                                                   permutation="random",
                                                   deltaf=1.0) {

  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")
  retlist <- .Call("rg_fleischer_max_concurrent_flow_c_prescaled",
                   as.double(cap),
                   as.integer(em-1),
                   as.integer(nv),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e)
                 )

  demands <- retlist$demands

  delta <- (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))
  
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  beta <- calcBeta(demands,gdual)
  
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3

  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound
                   )
}

### As above but using prescaling to speed it up
rg.minimum.congestion.flow.c <- function(g,
                                     demands,
                                     e=0.1,
                                     updateflow=TRUE,
                                     progress=FALSE,
                                     permutation="random",
                                     deltaf=1.0) {

  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")
  retlist <- .Call("rg_min_congestion_flow",
                   as.double(cap),
                   as.integer(em-1),
                   as.integer(nv),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e)
                 )

  demands <- retlist$demands

  delta <- (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))
  
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  
  res$gamma <- rg.mcf.find.gamma(gflow)


  ##beta <- calcBeta(demands,gdual)
  beta <- retlist$beta
  ##lambda <- calcLambda(demands)
  lambda <- retlist$lambda
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3

  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound
                   )
}

rg.minimum.congestion.flow.int.c <- function(g,
                                     demands,
                                     e=0.1,
                                     updateflow=TRUE,
                                     progress=FALSE,
                                     permutation="random",
                                     deltaf=1.0) {

  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")
  retlist <- .Call("rg_min_congestion_flow_int",
                   as.double(cap),
                   as.integer(em-1),
                   as.integer(nv),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e)
                 )

  demands <- retlist$demands

  delta <- (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))
  
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  
  res$gamma <- rg.mcf.find.gamma(gflow)


  beta <- calcBeta(demands,gdual)
  
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3

  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound
                   )
}



rg.fleischer.max.concurrent.flow.c.old <- function(g,
                                               demands,
                                               e=0.1,
                                               updateflow=TRUE,
                                               progress=FALSE,
                                               permutation="random",
                                               deltaf=1.0) {


  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)
  
  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  cap <- as.double(edgeData(g,attr="capacity"))

  demands.sources <- as.integer(lapply(demands,"[[","source"))
  demands.sinks <- as.integer(lapply(demands,"[[","sink"))
  demands.demand <- lapply(demands,"[[","demand")

  doubreq <- 2/e * log(ne/(1-e),base=(1+e))

  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                        max = doubreq, style=3)
  } else {
    pb <- NULL
  }
  if(is.numeric(permutation))
    permutation <- permutation - 1
  else if(permutation == "fixed")
    permutation <- 0: (length(demands) -1)
  else if(permutation == "random")
    permutation <- -1
  else if(permutation == "lowest")
    permutation <- -2
  
  #permutation <- seq(0,length(demands)-1)

  retlist <- .Call("rg_fleischer_max_concurrent_flow_c_boost",
                   as.integer(nv),
                   as.integer(ne),
                   as.integer(em-1),
                   as.double(eW),
                   as.double(cap),
                   as.integer(length(demands)),
                   as.integer(demands.sources -1 ),
                   as.integer(demands.sinks -1),
                   as.double(demands.demand),
                   as.double(e),
                   as.logical(updateflow),
                   pb,
                   environment(),
                   progress,
                   permutation,
                   deltaf
                 )

  demands <- retlist$demands

  delta <- deltaf * (ne / (1-e)) ^ (-1/e) 

  scalef <- 1 / log(1/delta,base=(1+e))

  
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist$lengths

  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  beta <- calcBeta(demands,gdual)
  
  lambda=NULL
  if(updateflow) {
    lambda <- calcLambda(demands)
    foundratio <- beta / lambda
    ratiobound <- (1-e)^-3
  } else {
    foundratio <- NULL
    ratiobound <- NULL
  }
  

  if( progress != FALSE) {
    close(pb)
  }

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist$totalphases,e=e,
                   ratio=foundratio,
                   bound=ratiobound
                   )
}

### Utility function used by rg.fleischer.max.concurrent.flow()
### penult - vetor from Dijkstra.Sp
### mincap - mincapacity to set flow to this value
### demand - the specific demand to update
updateExplicitFlow <- function(g,penult,mincap,demand) {
  
  s <- demand$source
  t <- demand$sink
  p <- t
  t <- penult[[t]]
  tag <- paste(nodes(g)[t],"|",p,sep="")
  
  path <- tag
  
  while(s !=t ) {
    p <- t
    t <- penult[[t]]
    path <- paste(nodes(g)[t],"|",path,sep="")
  }
  
  if(is.null(demand$paths[[path]])) {
    demand$paths[[path]] <- mincap
  } else {
    demand$paths[[path]] <- demand$paths[[path]] + mincap
    
  }
  
  demand
}


### Native R version of rg.fleischer.max.concurrent.flow.c()
### Compute the maximum concurrent flow See Fleischer "Approximating
### fractional multicommodity flow independent of the number of
### commodities", SIAM J. Discrete Maths, Vol 13/4, pp 505-520
###
### g: graphNEL (vertex labels need to be same as vertex index)
### demands: vector
###
### g: graphNEL with capacity data on each edge, edge
### weight is ignored (traffic is routed over any way to meet
### concurrent flow
###
### e: accuracy ( 0 < e < 1 ) 0 is more accurate, 1 less accurate
###
### output: list:
###              gflow: with weight set to edge flow
###              gdual: with weight set to dual edge lengths
###              demands: each with met demand and paths and edges

rg.fleischer.max.concurrent.flow <- function(g,demands,e=0.1,updateflow=TRUE, progress=FALSE) {


  savedemands <- demands
  ## note this is not the dual value
  calcD <- function() {
    sum(as.double(edgeData(gdual,attr="capacity"))*
        as.double(edgeData(gdual,attr="weight")))
  }


  
  doubleCount <- 0
  doubleDemands <- function(demands) {
    for(n in names(demands)) {
      demands[[n]]$demand <- demands[[n]]$demand * 2
    }
    doubleCount <- doubleCount + 1
    demands
  }
  gdual <- g
  # number of arcs
  m <- length(rg.edgeL(g))
  # number of nodes
  N <- length(nodes(g))
  delta <- (m / (1-e)) ^ (-1/e)
  capacities <- as.double(edgeData(g,attr="capacity"))
  fromedges <- edgeMatrix(g)[1,]
  toedges <- edgeMatrix(g)[2,]

  edgeData(gdual,
           nodes(gdual)[fromedges],
           nodes(gdual)[toedges],
           "weight") <- delta / capacities

  for(c in names(demands)) {
    demands[[c]]$paths <- list()
    demands[[c]]$flow <- 0
  }

  
  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))
  ## updatepb used to measure progress as percent
  updatepb <- as.integer(ceiling(doubreq / 100.0))

  if(progress != FALSE)
    pb <- txtProgressBar(min = 0,
                         max = doubreq, style=3)
  phases <- 0
  totalphases <- 0
  D <- calcD()
  lambdav <- c()
  betav <- c()
  while(D < 1 ) {
    if(progress != FALSE && totalphases %% updatepb == 0) {
      setTxtProgressBar(pb,totalphases)
    }
    if(phases > doubreq) {
      cat("doubling\n");
      demands <- doubleDemands(demands)
      phases <- 0
    }

    ccount <- 1
    for(c in demands) {
      demand <- c$demand
      D <- calcD()
      while( D < 1 && demand > 0 ) {
        sp <- dijkstra.sp(gdual,c$source)
        p <- extractPath(which(nodes(gdual)==c$source),
                               which(nodes(gdual)==c$sink),sp$penult)
        caponpath <- as.double(edgeData(gdual,
                                        from=nodes(gdual)[p[1:length(p)-1]],
                                        to=nodes(gdual)[p[2:length(p)]],
                                        attr="capacity"))

        mincap <- min(demand,caponpath)

        demand <- demand - mincap
        
        lengths <- as.double(edgeData(gdual,
                                      from=nodes(gdual)[p[1:length(p)-1]],
                                      to=nodes(gdual)[p[2:length(p)]],
                                      attr="weight"))
        lengths <- lengths * (1 + (e*mincap) / caponpath)
        edgeData(gdual,
                 from=nodes(gdual)[p[1:length(p)-1]],
                 to=nodes(gdual)[p[2:length(p)]],
                 attr="weight") <- lengths
        if(updateflow)
          c <- updateExplicitFlow(gdual,sp$penult,mincap,c)
        c$flow <- c$flow + mincap
        D <-  calcD()
      }
      demands[[ccount]] <- c
      ccount <- ccount + 1
    }
    phases <- phases +1
    totalphases <- totalphases + 1

  }
  ## If there had been demands doubling we need to
  ## fix things for the returned demands
  for(n in names(demands)) {
    demands[[n]]$demand <- savedemands[[n]]$demand
  }

  scalef <- 1 / log(1/delta,base=(1+e))
  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)

  beta <- calcBeta(demands,gdual)

  lambda=NULL
  if(updateflow) {
    lambda <- calcLambda(demands)
      foundratio <- beta / lambda
    ratiobound <- (1-e)^-3
    
  }
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,phases=totalphases,e=e)

  if(progress != FALSE)
    close(pb)

  retval
}

### Calculates the maximum concurrent flow as per Fleischer, but with the
### modification that demand paths have been specified
### g - graphNEL with edge attributes capacity and weight (weight not used)
### demands$flow
### demands$paths - list of paths each element list[["1|2|3"]] = flow on path
###                 only the path is used not the amount of flow that
###                 may have been calculated on a previous run
### e - the approximation limit
### updateflow - default TRUE, may run a bit faster if FALSE but only
###              lambda, bega and demand flow are given (path values not
###              updated.
### progress - prints out the value of D() (finishes when D >=1.0)
### returns list:
### demands paths - list of paths each element list[["1|2|3"]] = flow on path
### gflow - graph as g but with weight set to traffic flow
### gdual - the dual graphNEL with weight set to edge lengths
### beta - the dual value
### lambda - the primal value
### phases - the total number of phases
### e - as per the input
rg.fleischer.max.concurrent.flow.restricted <- function(g,demands,e=0.1,updateflow=TRUE, progress=FALSE) {


  ## this is actually slower than below!
  findshortestpathtest <- function(paths,vlength) {
    sums <- sapply(paths,function(x) sum(vlength[x]))
    return(which.min(sums))
  }

  findshortestpath <- function(paths,vlength) {
    min <- Inf
    minp <- NULL
    minpcount <- NULL
    for(i in 1:length(paths)) {
      if ( sum(vlength[paths[[i]]]) < min) {
        minp <- paths[[i]]
        minpcount <- i
        min <- sum(vlength[paths[[i]]])
      }
    }
    
    return(minpcount)
  }
  
  savedemands <- demands
  
  vdemands <- as.double(lapply(demands,"[[","demand"))
  vcapacity <- as.double(edgeData(g,attr="capacity"))
  vweight <- rep(0.0,length(vcapacity))
  vlength <- rep(0.0,length(vcapacity))
  vdemandflow <- rep(0.0,length(vcapacity))
  L <- length(vlength)
  edgeMap <- adjMatrix(g)
  demandpaths <- list()
  demandpathflows <- list()
  j <- 1
  for(d in demands) {
    demandpaths[[j]] <- list()
    demandpathflows[[j]] <- list()
    k <- 1
    for(p in names(d$paths)) {
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- as.integer(pv[1:length(pv)-1])
      tolist <- as.integer(pv[2:length(pv)])
      me <- rbind(fromlist,tolist)
      pv <- c()
      for(i in 1:length(fromlist)) {
        v <- as.vector(me[,i])
        em <- edgeMap[v[1],v[2]]
        pv <- append(pv,em)
      }
      demandpaths[[j]][[k]] <- pv
      demandpathflows[[j]][[k]] <- 0.0
      k <- k + 1
    }
    j <- j + 1
  }
  ## note this is not the dual value
  calcD <- function() {
    sum(vcapacity*vlength)
  }
  
  doubleCount <- 0
  doubleDemands <- function(demands) {
    vdemands <- vdemands * 2.0
    doubleCount <- doubleCount + 1
    vdemands
  }
  gdual <- g

  ## number of arcs
  m <- length(rg.edgeL(g))
  ## number of nodes
  N <- length(nodes(g))
  delta <- (m / (1-e)) ^ (-1/e)
  vlength <- delta / vcapacity

  for(c in names(demands)) {
    demands[[c]]$paths <- list()
    demands[[c]]$flow <- 0
  }

  
  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))
  ## updatepb used to measure progress as percent
  updatepb <- as.integer(ceiling(doubreq / 100.0))

  if(progress != FALSE)
    pb <- txtProgressBar(min = 0,
                        max = doubreq, style=3)
  phases <- 0
  totalphases <- 0
  D <- calcD()

  ## main algorithm
  while(D < 1 ) {
    if(progress != FALSE && totalphases %% updatepb == 0) {
      setTxtProgressBar(pb,totalphases)
    }
    if(phases > doubreq) {
      cat("doubling",doubreq,totalphases,"\n");
      vdemands <- doubleDemands(vdemands)
      phases <- 0
    }

    for(i in 1:length(vdemands)) {
      demand <- vdemands[i]
      D <- calcD()
      while( D < 1 && demand > 0 ) {
        p <- findshortestpath(demandpaths[[i]],vlength)
        caponpath <- vcapacity[ demandpaths[[i]][[p]] ]
        mincap <- min(demand,caponpath)
        demand <- demand - mincap
        lengths <- vlength[ demandpaths[[i]][[p]] ]
        lengths <- lengths * (1 + (e*mincap) / caponpath)
        vlength[ demandpaths[[i]][[p]] ] <- lengths
        demandpathflows[[i]][[p]] <- demandpathflows[[i]][[p]] + mincap
        vdemandflow[i] <- vdemandflow[i] + mincap
        D <-  calcD()
      }
    }
    phases <- phases +1
    totalphases <- totalphases + 1
  }
  ## end of main algorithm - now need to pack results

  ## If there had been demands doubling we need to
  ## fix things for the returned demands
  for(n in 1:length(demands)) {
    demands[[n]]$demand <- savedemands[[n]]$demand
    demands[[n]]$paths <- savedemands[[n]]$paths
    demands[[n]]$flow <- vdemandflow[n]
    for(i in 1:length(demands[[n]]$paths)) {
      demands[[n]]$paths[[i]] <- demandpathflows[[n]][[i]]
    }
  }
  gdual <- rg.set.weight(gdual,vlength)
  scalef <- 1 / log(1/delta,base=(1+e))
  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)

  beta <- calcBeta(demands,gdual)
  betar <- calcBetaRestricted(demands,gdual)
  lambda=NULL
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3
    
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,betar=betar,
                 lambda=lambda,phases=totalphases,e=e,vlength=vlength)

  if(progress != FALSE)
    close(pb)

  retval
}

### As the R version see above Tested 21/5/2010 performs the same,
### however some rounding differences mean there will be differences
### in actual result from R version. A small rounding error can
### influence the route (shortest path) so differences could be
### substantial, however, lambda and beta should be within the bounds
### set by e
rg.fleischer.max.concurrent.flow.restricted.c <- function(g,demands,e=0.1,updateflow=TRUE,progress=FALSE,bestgamma=-1) {

  ## note this is not the dual value
  calcD <- function() {
    sum(vcapacity*vlength)
  }

  savedemands <- demands
  
  vdemands <- as.double(lapply(demands,"[[","demand"))
  vcapacity <- as.double(edgeData(g,attr="capacity"))
  vlength <- rep(0.0,length(vcapacity))
  vdemandflow <- rep(0.0,length(vcapacity))
  L <- length(vlength)
  delta <- (L / (1-e)) ^ (-1/e)
  vlength <- delta / vcapacity
  edgeMap <- adjMatrix(g)

  demands.paths <- as.integer(c())

  ## code the paths as integers
  ## [demandno(e.g=1) numdempaths lengthpath1 path1[1] path1[2] ...
  ##  lengthpath2 path2[1] path2[2]....demandno(e.g=2)....]
  rdemandpaths <- list();
  for(d in 1:length(demands)) {
    rdemandpaths[[d]] <- list()
    demands.paths <- append(demands.paths,d-1)
    demands.paths <- append(demands.paths,length(demands[[d]]$paths))
    for(p in 1:length(demands[[d]]$paths)) {
      pathi <-
        as.integer(strsplit
                   (names(demands[[d]]$paths[p]),"|",fixed=TRUE)[[1]])
      pathi <- pathi -1
      rdemandpaths[[d]][[p]] <- pathi
      demands.paths <- append(demands.paths,length(pathi))
      demands.paths <- append(demands.paths,pathi)
    }

  }

  demandpaths <- list()
  j <- 1
  for(d in demands) {
    demandpaths[[j]] <- list()
    k <- 1
    for(p in names(d$paths)) {
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- as.integer(pv[1:length(pv)-1])
      tolist <- as.integer(pv[2:length(pv)])
      me <- rbind(fromlist,tolist)
      pv <- c()
      for(i in 1:length(fromlist)) {
        v <- as.vector(me[,i])
        em <- edgeMap[v[1],v[2]]
        pv <- append(pv,em)
      }
      demandpaths[[j]][[k]] <- pv -1
      k <- k + 1
    }
    j <- j + 1
  }

  
  
  m <- length(rg.edgeL(g))

  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))

  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = doubreq, style=3)
  } else {
    pb <- NULL
  }

  retlist <- .Call("rg_fleischer_max_concurrent_flow_restricted_c",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   e,
                   progress,
                   pb,
                   bestgamma,
                   environment()
                   );

  if( progress != FALSE) {
    close(pb)
  }

  ## now need to unpack results
  demands <- savedemands
  
  for(n in 1:length(demands)) {
    flow <- 0
    for(i in 1:length(demands[[n]]$paths)) {
      flow <- flow + retlist$pathflows[[n]][[i]]
      demands[[n]]$paths[[i]] <- retlist$pathflows[[n]][[i]]
    }
    demands[[n]]$flow <- flow
  }

  gdual <- g
  gdual <- rg.set.weight(gdual,retlist$vlengths)
  scalef <- 1 / log(1/delta,base=(1+e))
  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  beta <- calcBeta(demands,gdual)
  betar <- calcBetaRestricted(demands,gdual)
  lambda=NULL
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3
    
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,betar=betar,
                 lambda=lambda,phases=retlist$totalphases,e=e,vlength=retlist$vlengths,
                 countgamma=retlist$countgamma,
                 bestgamma=retlist$bestgamma,
                 bestpaths=retlist$bestpaths+1)

  
  return(retval)
}



rg.max.concurrent.flow.analyse <- function(res) {

  demands <- res$demands
  for(d in demands) {
    numpaths = length(d$paths)
    cat("demand source =",d$source,"sink =",d$sink,"paths = ",numpaths,"\n")
  }

}

rg.fleischer.max.concurrent.flow.stats <- function(reslist) {

  gdual <- reslist$gdual
  D <- sum(as.double(edgeData(gdual,attr="capacity"))*
           as.double(edgeData(gdual,attr="weight")))

  cat("D=",D,"\n",sep="")
  cat("in",reslist$phases,"phases\n")
  e <- reslist$e
  

  beta <- calcBeta(reslist$demands,gdual)
  cat("beta=",beta,"\n",sep="")

  lambda <- calcLambda(reslist$demands)
  cat("lambda=",lambda,"\n",sep="")
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3
  cat("Ratio dual/primal found =",foundratio,
      " Ratio should be less than =",ratiobound,"\n",sep="")
  m <- length(rg.edgeL(gdual))

  doubreq <- 2/e * log(m/(1-e),base=(1+e))
  if(doubreq < reslist$phases)
    cat("Doubling was required!\n")
}

### Rescale the demands by scalef
rg.max.concurrent.flow.rescale.demands.flow <- function(demands,scalef) {
  for(dn in names(demands)) {
    demands[[dn]]$flow <- demands[[dn]]$flow * scalef
    for(en in names(demands[[dn]]$paths)) {
      demands[[dn]]$paths[[en]] <- demands[[dn]]$paths[[en]] * scalef
    }
  }
  demands
}

rg.rescale.demands <- function(demands,scalef,integer=FALSE) {

  for(dn in names(demands)) {
    if(integer) {
      demands[[dn]]$demand <- round(demands[[dn]]$demand * scalef)
    } else {
      demands[[dn]]$demand <- demands[[dn]]$demand * scalef
    }
  }
  demands
}

rg.set.all.graph.edge.weights <- function(g,val=0.0)  {
  fromedges <- edgeMatrix(g)[1,]
  toedges <- edgeMatrix(g)[2,]
  nodes=nodes(g)
  edgeData(g,
           nodes[fromedges],
           nodes[toedges],
           "weight") <- val
  g
}


### Utililty function used by rg.fleischer.max.concurrent.flow
### given a graph with no edge weights set, set the edge
### weights given edge flows in demand
###
### graph - GraphNEL with empty edge weights (or ignored if they are set)
### demands - each demand path value used to update graph edge weight
### return - new graph with edge weights set

rg.max.concurrent.flow.graph <- function(g,demands) {
  g <- rg.set.all.graph.edge.weights(g)
  for(d in demands) {
    for(p in names(d$paths)) {
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- pv[1:length(pv)-1]
      tolist <- pv[2:length(pv)]
      
      weights <- as.double(edgeData(g,from=fromlist,
                                    to=tolist,att="weight"))
      weights <- weights + d$paths[[p]]
      edgeData(g,from=fromlist,
               to=tolist,
               att="weight") <- weights

    }
  }
  g
}

### Calculate the gamma of a graph (minimum proportion of free edge capacity)
### gflow - GraphNEL with edge attr weight and capacity set
### lambda - defaults to 1 else divide weight by lambda
###          this is useful for output of max concurrent flow
###          use lambda=1 for gflow from minimim.congestion.flow
rg.mcf.find.gamma <- function(gflow,lambda=1.0) {
  flows <- as.double(edgeData(gflow,attr="weight"))/lambda
  capacities <- as.double(edgeData(gflow,attr="capacity"))

  gamma <- min ( (capacities - flows)/ capacities)
  gamma
}

### Compute the minimum congestion flow using the max concurrent flow
### and rescale demands by 1/lambda
### g - GraphNEL with attr capacity
### demands - list with source, sink, demand
### e - approximation 0 < e < 1.0
### progress - bool show progress bar default FALSE
### return list
###            gflow - GraphNEL with attr weight set
###            demands - list with paths and edges of solution
###            gamma - minimum proportion of free edge capacity
rg.minimum.congestion.flow <- function(g,demands,e=0.1,progress=FALSE,permutation="random",deltaf=1.0,ccode=TRUE) {
    if(!isDirected(g)) {
        print("g must be directed")
        return(NULL)
    }
  res <- rg.max.concurrent.flow.prescaled(g,demands,e,progress=progress,ccode,permutation=permutation,deltaf=1.0)
  res$demands <- rg.max.concurrent.flow.rescale.demands.flow(res$demands,1/res$lambda)
  
  res$gflow <- rg.max.concurrent.flow.graph(res$gflow,res$demands)
  
  res$gamma <- rg.mcf.find.gamma(res$gflow)
  
  res
}

### Truncated quantile function for an arbitrary distribution
### spec = distribution name
### a = negative truncation (must be greater than this)
### b = positive truncation (must be less/= to this)
### .... arbitrary arguments to send to the distribution
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
    tt <- p
    G <- get(paste("p", spec, sep = ""), mode = "function")
    Gin <- get(paste("q", spec, sep = ""), mode = "function")
    tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
    return(tt)
}

### Random number generator using a random distribution which
### is optionally truncated
### spec = name of distribution
### a = lower truncation
### b = upper truncation
### ... arguments to be sent to probability distibution
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
    x <- u <- runif(n, min = 0, max = 1)
    x <- qtrunc(u, spec, a = a, b = b,...)
    return(x)
}

### Generate a random graph with a range of node degree
### uses igraph routines
### n - number of nodes
### mindeg - minimum node degree default = 2
### maxdeg - maximum node degree default = 20
### dist - the name of the distrubution for the outdegree
### dist1 - the first parameter used by dist
### dist2 - the second parameter used by dist
### return GraphNEL
rg.generate.random.graph <- function(n=10,mindeg=2,maxdeg=20,
                                     dist="weibull",
                                     dist1=0.42,dist2=1,retry=10) {
  if(n< maxdeg -1) {
    maxdeg=n-2
  }
  ## note this needs a bit more work - not quite right
  ## as the truncation gives an error - but good enough
  ## for generating some "random" graphs
  if(maxdeg==Inf)
      degrees <- ceiling(dist(n,dist1,dist2))
  else
      degrees <- floor(rtrunc(n,dist,a=mindeg,b=maxdeg,dist1,dist2))

  ## note sum(degrees) has to be even, add one to smallest
  ## node degree of first node found.
  if(sum(degrees) %% 2 != 0) {
    minval <- min(degrees)
    pos <- match(minval,degrees)
    degrees[pos] <- degrees[pos] +1
  }

  gi <- NULL
  tryCatch(gi <- degree.sequence.game(degrees,method="vl"),
           error=function(e) { print("trying graph generation again") })
  if(is.null(gi)) {
    for(i in 1:10) {
      degrees <- floor(rtrunc(n,dist,a=mindeg,b=maxdeg,dist1,dist2))
      
      ## note sum(degrees) has to be even, add one to smallest
      ## node degree of first node found.
      if(sum(degrees) %% 2 != 0) {
        minval <- min(degrees)
        pos <- match(minval,degrees)
        degrees[pos] <- degrees[pos] +1
      }
      tryCatch(gi <- degree.sequence.game(degrees,method="vl"),
               error=function(e) { print("trying graph generation again") })
      if(!is.null(gi)) break
    }
  }
  gi <- as.directed(gi)
  # Finally convert to RBGL - my preferred library for graphs
  G <- rg.relabel(igraph.to.graphNEL(gi))
  return(G)
}


rg.augment.graph.for.wavelengths <- function # Augment graph for wavelengths
### Augments a directed graph with a copy for each wavelength. Source and sink
### nodes are created to attach the demands to.
(g,    ##<< graphNEL object to augment
 count ##<< number of wavelengths (same on each link)
 ) {
  ## create as many graph as there are wavelengths as simple copies,
  ## they are not connected at the moment
  edgeL <- list()
  nodeL <- c()
  linkgroupmap <- list()
  for( i in 1:count ) {
    inedgeL <- mapply(paste0,graph::edges(g),MoreArgs=list("L",i),SIMPLIFY=FALSE)
    innodeL <- as.vector(paste0(nodes(g),"L",i))
    names(inedgeL) <- innodeL
    edgeL <- append(edgeL,inedgeL)
    nodeL <- append(nodeL,innodeL)
  }
  edgenames <- names(edgeData(g))
  for(e in edgenames) {
    for( i in 1:count ) {
      ed <- sapply(strsplit(e,"\\|"),function(x) paste0(x,"L",i))[,1]
      circuit.name <- paste0(ed[[1]],"|",ed[2])
      linkgroupmap[[circuit.name]] <- e
    }
  }

  #nodeL <- append(nodeL,nodes(g))
  ga <- new("graphNEL",nodeL,edgeL,"directed")
  edgeDataDefaults(ga,attr="capacity") <- 1.0
  ga <- rg.set.weight(ga,1.0)
  ga <- rg.set.capacity(ga,1.0)
  ## now connect the graphs to nodes of the same name as the original nodes
  ## but with nodes split into source and sink by appending s and t to names
  ga <- addNode(as.vector(paste0(nodes(g),"s")),ga)
  ga <- addNode(as.vector(paste0(nodes(g),"t")),ga)
  for(v in nodes(g)) {
    for( i in 1:count ) {
      node <- paste0(v,"L",i)
      ga <- addEdge(paste0(v,"s"),node,ga)
      edgeData(ga,paste0(v,"s"),node,attr="capacity") <- as.integer(graph::degree(g,v)$inDegree)
      ga <- addEdge(node,paste0(v,"t"),ga)
      edgeData(ga,node,paste0(v,"t"),attr="capacity") <- as.integer(graph::degree(g,v)$inDegree)
    }
  }
  return(list(ga=ga,linkgroupmap=linkgroupmap))
  ### graphNEL object augmented with count copies of the original
  ### graph, an original node n for the wavelength l is labelled as
  ### nLl. Additional nodes nt and ns are created for the source/sink
  ### at the orignal node n. Now we have edges ns -> nLl and nLl -> nt
  ### as well as the original edges between the nodes (albeit in
  ### multiple copies of the graph) This is better understood using a
  ### picture...
}




### Calculate the path cost in the a graph
### gdual - GraphNEL edge attr weight used for path cost
### path - path as "|" separated nodes e.g. "14|7|5"
### return - path cost as double
rg.path.cost <- function(gdual,path) {

  pv <- as.vector(strsplit(path,"|",fixed=TRUE)[[1]])
  from <- pv[1:length(pv)-1]
  to <- pv[2:length(pv)]
  cost <- sum(as.double(edgeData(gdual,from=from,to=to,att="weight")))
  return(cost)
  
}

### Obtain the edge attribute along the path cin the a graph
### g - GraphNEL 
### path - path as "|" separated nodes e.g. "14|7|5"
### return - attribute vector same length as path
rg.path.attr <- function(g,path,attr="weight") {

  pv <- as.vector(strsplit(path,"|",fixed=TRUE)[[1]])
  from <- pv[1:length(pv)-1]
  to <- pv[2:length(pv)]
  return(as.double(edgeData(g,from=from,to=to,att=attr)))
}

### Attempt to calculate integer min congestion flow
### does not work!
### res - output from rg.integer.min.congestion.flow.ga()
rg.integer.iterative <- function(res) {
  dem <- rg.max.concurrent.flow.rescale.demands.flow(res$splitflow$demands,res$splitflow$lambda)

  g <- res$splitflow$gflow

  
  e <- res$splitflow$e
  
  fromedges <- edgeMatrix(g)[1,]
  toedges <- edgeMatrix(g)[2,]

  edgeData(g,
           as.character(fromedges),
           as.character(toedges),
           "weight") <- 0.0
  
  outdem <- list()
  demsav <- dem
  gdual <- res$splitflow$gdual

  gflow <- rg.max.concurrent.flow.graph(res$splitflow$gflow,dem)
  
  beta <- calcBeta(dem,gdual)
  betaa <- calcBetaRestricted(dem,gdual)
  lambda <- calcLambda(dem)

  gap <- lambda / beta

  mgap <- (1-res$splitflow$e)^3
  
  cat("Initial: beta=",beta,betaa,"lambda=",lambda,"gap=",gap,"mgap=",mgap,"\n")

  multFlows <- FALSE

  for(d in names(dem)) {
    if(length(dem[[d]]$paths) > 1 ) {
      multFlows <- TRUE
      break
    }
  }

  for(d in names(dem)) {
    dem[[d]]$original <- d
    dem[[d]]$paths <- NULL
  }
  
  while(length(dem) > 0) {
## calculate mcf
##    browser()
    mcf <- rg.fleischer.max.concurrent.flow.c(g,dem,e)
    demands <- mcf$demands
    gdual <-  mcf$gdual
    lowp <- ""
    val <- Inf
    lowd <- ""
    
    for(d in names(demands)) {
      if(length(demands[[d]]$paths) > 1 ) {
        for(p in names(demands[[d]]$paths)) {
          cost <- rg.path.cost(gdual,p)
          if (val > cost) {
            val <- cost
            lowp <- p
            lowd <- d
          }
          flow <- demands[[d]]$paths[[p]]
        }
      }
    }

    cat("lowest demand=",lowd,"lowest path=",lowp,"\n")

    ## nail chosen demand on chosen path of output dem

    orignum <- dem[[lowd]]$original
    
    outdem[[orignum]]$source <- dem[[lowd]]$source
    outdem[[orignum]]$sink <- dem[[lowd]]$sink
    outdem[[orignum]]$demand <- dem[[lowd]]$demand
    outdem[[orignum]]$paths <- list()
    outdem[[orignum]]$paths[[lowp]] <- dem[[lowd]]$demand

    ## DO WE NEED TO FILL IN EDGES AS WELL?
    ## update graph

    pv <- as.vector(strsplit(lowp,"|",fixed=TRUE)[[1]])

    fromlist <- pv[1:{length(pv)-1}]
    tolist <- pv[2:{length(pv)}]
#    print(fromlist)
#    print(tolist)
    oldvals <- as.double(edgeData(g,
                                  as.character(fromlist),
                                  as.character(tolist),
                                  "capacity"))
    newvals <- oldvals - outdem[[orignum]]$demand

    if( min(newvals) < 0 ) {
      print("WARNING SET EDGE LESS THAN ZERO")
    }
    
    ## remove from working demands
    dem[[lowd]] <- NULL

    if (length(dem) > 0) {
      names(dem) <- as.character(1:length(dem))
    }
  }

  gflow <- rg.set.all.graph.edge.weights(res$splitflow$gflow)
  
  
  for(d in names(outdem)) {
    p <- names(outdem[[d]]$paths[1])
    pathv <- strsplit(p,"|",fixed=TRUE)
    from <- pathv[[1]][1:length(pathv[[1]])-1]
    to <- pathv[[1]][2:length(pathv[[1]])]
    
    edgeData(gflow,from=from,to=to,attr="weight") <-
      as.double(edgeData(gflow,from=from,to=to,attr="weight")) + outdem[[d]]$demand
  }


  ga.gamma <- rg.mcf.find.gamma(res$gflow)
  gamma <- rg.mcf.find.gamma(gflow)

  cat("ga.gamma=",ga.gamma,"gamma=",gamma,"\n")
  res$demands <- outdem
  res$gflow <- gflow
  
  return(res)
}

### Attempt to calculate integer min congestion flow
### does not work!
### res - output from rg.integer.min.congestion.flow.ga()
rg.integer.rounding <- function(res) {


  dem <- rg.max.concurrent.flow.rescale.demands.flow(res$splitflow$demands,res$splitflow$lambda)


  gdual <- res$splitflow$gdual

  gflow <- rg.max.concurrent.flow.graph(res$splitflow$gflow,dem)
  
  beta <- calcBeta(dem,gdual)
  betaa <- calcBetaRestricted(dem,gdual)
  lambda <- calcLambda(dem)

  gap <- lambda / beta

  mgap <- (1-res$splitflow$e)^3
  
  cat("Initial: beta=",beta,betaa,"lambda=",lambda,"gap=",gap,"mgap=",mgap,"\n")

  multFlows <- FALSE

  for(d in names(dem)) {
    if(length(dem[[d]]$paths) > 1 ) {
      multFlows <- TRUE
      break
    }
  }
  
  while(multFlows) {
    
    lowp <- ""
    val <- Inf
    lowd <- ""
    
    for(d in names(dem)) {
      if(length(dem[[d]]$paths) > 1 ) {
        data <- data.frame(Paths=character(0),cost=numeric(0),flow=numeric(0),highest=numeric(0))
        count <-  1
        lowc <- 1
        for(p in names(dem[[d]]$paths)) {
          cost <- rg.path.cost(gdual,p)
          if (val > cost) {
            val <- cost
            lowp <- p
            lowc <- count
            lowd <- d
          }
          flow <- dem[[d]]$paths[[p]]
          ##        data <- rbind(data,data.frame(Paths=p,cost=cost,flow=flow,highest=0))
          count <- count + 1
        }
        
        ##      data[[lowc,4]] <- 1
        ##      data <- data[do.call(order,data[3]),]
        ##    print(data)
      }
    }

    cat("lowest demand=",lowd,"lowest path=",lowp,"\n")
    
    ## work out where to put flow
    ## the one that make gamma highest
    gt <- Inf
    gp <- ""
    lowestload <- 0
    free <- 0
    for(p in names(dem[[lowd]]$paths)) {
      if (p != lowp) {
        pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
        from <- pv[1:length(pv)-1]
        to <- pv[2:length(pv)]
                                        #        load <- as.double(edgeData(gflow,from=from,to=to,att="weight"))
                                        #        capacities <- as.double(edgeData(gflow,from=from,to=to,att="capacity"))
                                        #        print(p)
        maxratio <- 0
        free <- 0
        for(i in seq(1:length(from))) {
                                        #          cat(from[i],to[i],"\n")
          load <- as.double(edgeData(gflow,from=from[i],to=to[i],att="weight"))
          capacity <- as.double(edgeData(gflow,from=from[i],to=to[i],att="capacity"))
          if (load/capacity > maxratio) {
            maxratio <- load/capacity
            free <- capacity - load
 #           if (free < 0) {
 #             print("Error free less than zero")
 #             return(0)
 #           }
          }
        }
        
                                        #lowestload <- max(load/capacities)
        lowestload <- maxratio
        if(lowestload < gt) {
          gt <- lowestload
          gp <- p
        }
      }
    }
    
    ## update flow
    
    lowflow <- dem[[lowd]]$paths[[lowp]]
    
    addval <- lowflow
#    if (addval > free) {
#      addval <- free
#      dem[[lowd]]$flow <- dem[[lowd]]$flow - lowflow + addval
#      
#    }
    
    ## first add flow to lowest constrained path in demand
    pv <- as.vector(strsplit(gp,"|",fixed=TRUE)[[1]])
    from <- pv[1:length(pv)-1]
    to <- pv[2:length(pv)]
    
    weights <- as.double(edgeData(gflow,
                                  from=from,
                                  to=to,
                                  attr="weight"))

    dem[[lowd]]$paths[[gp]] <- dem[[lowd]]$paths[[gp]] + addval
    
    weights <- weights + addval
    
    edgeData(gflow,
             from=from,
             to=to,
             attr="weight") <- weights
    

    ## Now take it away from the one that was the lowest path cost
    pv <- as.vector(strsplit(lowp,"|",fixed=TRUE)[[1]])
    from <- pv[1:length(pv)-1]
    to <- pv[2:length(pv)]
    
    weights <- as.double(edgeData(gflow,
                                  from=from,
                                  to=to,
                                  attr="weight"))
    
    weights <- weights - lowflow
    edgeData(gflow,
             from=from,
             to=to,
             attr="weight") <- weights
    
    dem[[lowd]]$paths[[lowp]] <- NULL
    
    ## update lengths on the lowest path cost
    ##        lengths <- lengths * (1 + (e*mincap) / caponpath)
    if(FALSE) {
    lengths <- as.double(edgeData(gdual,
                                  from=from,
                                  to=to,
                                  attr="weight"))
    caponpath <- as.double(edgeData(gdual,
                                    from=from,
                                    to=to,
                                    attr="capacity"))
    
    lengths <- lengths * ( 1 + addval/caponpath)
    
    edgeData(gdual,
             from=from,
             to=to,
             attr="weight") <- lengths
  }
    
    ## calculate beta, gamma and lambda
    beta <- calcBetaRestricted(dem,gdual)
    
    lambda <- calcLambda(dem)
    
    gap <- lambda / beta
    
    mgap <- (1-res$splitflow$e)^3
    
    cat("beta=",beta,"lambda=",lambda,"gap=",gap,"mgap=",mgap,"\n")

    multFlows <- FALSE
    for(d in names(dem)) {
      if(length(dem[[d]]$paths) > 1 ) {
        multFlows <- TRUE
        break
      }
    }
  }

  for(d in names(dem)) {
    dem[[d]]$flow <- dem[[d]]$demand
    dem[[d]]$paths[1] <- dem[[d]]$demand
  }

  
  gflow <- rg.max.concurrent.flow.graph(gflow,dem)

  load <- as.double(edgeData(gflow,attr="weight"))
  capacity <- as.double(edgeData(gflow,attr="capacity"))
  print(load/capacity)
  gamma(1- max(load/capacity))
  
        
  dem$gflow <- gflow
  return(dem)
}

### Working function to calculate stats from a run of results
### from rg.integer.min.congestion.flow.ga()
### path - file path
rg.test.idea.stats <- function(path) {
  files <- Sys.glob(paste(path,"run[0-9]*",sep=""))
  files <- files[order(as.numeric(gsub("^.*run([0-9]*).*$","\\1",files,files)))]
  for(file in files) {
    load(file)
    e <- res$splitflow$e
    n <- length(nodes(res$splitflow$gflow))
    nd <- length(res$demands)
    lambda <- res$splitflow$lambda
    gamma <- res$splitflow$gamma

    
    capacities <- as.double(edgeData(res$splitflow$gflow,attr="capacity"))

    dual <- res$splitflow$gdual

    

    lengths <- as.double(edgeData(dual,attr="weight"))
    
    lengths <- lengths / (1 - gamma)

    dual <- rg.set.all.graph.edge.weights(dual,lengths)

    numpaths <- 0
    for( i in lapply(res$splitflow$demands,"[[","paths")) {
      numpaths <- numpaths + length(names(i))
    }
    
    
    scap <- capacities / lambda

    min <- 1- max( scap /capacities)

    igamma <- rg.mcf.find.gamma(res$gflow)

    minEval <-  min(res$ga$evaluations)

    beta <- res$splitflow$beta
    
#    cat("n=",n,"e=",e,"nd=",nd,"lambda=",lambda,"\n")
    print(file)
    cat("n=",n,"nd=",nd,"beta=",beta,"lambda=",lambda,"gamma=",igamma,"numpaths=",numpaths,"\n")
#    cat("maxl=",maxl,"lambda=",lambda,"lambda2",res2$lambda,"gamma=",gamma,"igamma=",igamma,"min=",min,"\n")
  }
}


### Solve the integer min congestion flow using a GA
### uses the rg.minimum.congestion.flow() to give starting
### value. GA permutes various paths from this non-integer solution
### rather than using "all possible paths".
### uses rbga.bin() from the genalg package
### g - GraphNEL attr capacity set
### demands - list with source, sink, demand set for each demand
### e - approximation value used for the internal rg.minimum.congestion.flow()
### return - list
###        - gflow integer flow GraphNEL with attr weight set
###        - splitflow output from rg.minimum.congestion.flow
###        - val 1- minEval (minEval is minimum ga evaluation)
###        - ga output from rbga.bin
###        - demands integer solution with paths set (but not edges)
rg.integer.min.congestion.flow.ga <- function(g,demands,e=0.1) {


  monitor <- function(obj) {
##    minEval = min(1/obj$evaluations);
##    plot(obj, type="hist");
    print(min(obj$evaluations))
  }
  
  evaluate <- function(string=c()) {
    chromasome$string <- string
    paths <- rg.demands.paths.from.chromasome(chromasome)
    
    g <- rg.set.all.graph.edge.weights(chromasome$g)
    
    j <- 1
    valid <- TRUE
    for(p in paths) {
      demands[[j]]$paths <- list()
      demands[[j]]$paths[[p]] <- demands[[j]]$demand
      if ( p == "NULL" ) {
        valid <- FALSE
        break
      } else {
        d <- chromasome$demands[[j]]$demand
        j <- j + 1
        pathv <- strsplit(p,"|",fixed=TRUE)
        from <- pathv[[1]][1:length(pathv[[1]])-1]
        to <- pathv[[1]][2:length(pathv[[1]])]
        
        edgeData(g,from=from,to=to,attr="weight") <-
          as.double(edgeData(g,from=from,to=to,attr="weight")) + d
      }
    }
    if (valid) {
      load <- as.double(edgeData(g,attr="weight"))
      capacity <- as.double(edgeData(g,attr="capacity"))
      val <- max(load / capacity)
    } else {
      val <- Inf
    }
    val
  }
  cat("Calculating split flow\n")
  splitflow <- rg.minimum.congestion.flow(g,demands,e,progress=FALSE)

  chromasome <- rg.demands.paths.as.chromasome(splitflow$demands)
  
  chromasome$g <- g
  res <- list()

  cat("Starting GA\n")

  res$ga <- rbga.bin(size=chromasome$totallength,
                  popSize=200,
                  iters=100,
                  mutationChance=1/(chromasome$totallength + 1),
                  elitism=as.integer(chromasome$totallength * 0.2),
                  evalFunc=evaluate,
                  ##monitorFunc=monitor,
                  verbose=TRUE)
  cat("Done GA\n")

  save(list=ls(all=TRUE),file="rg.integer.min.congestion.flow.ga.robj")
  minEval <- min(res$ga$evaluations)
  filter <- res$ga$evaluations == minEval
  bestSolutions <- unique(res$ga$population[filter,])
  res$bestSolutions <- bestSolutions
  res$splitflow <- splitflow

  chromasome$string <- bestSolutions[1,]
  paths <- rg.demands.paths.from.chromasome(chromasome)
  g <- rg.set.all.graph.edge.weights(g)
    
  j <- 1
  valid <- TRUE
  for(p in paths) {
    demands[[j]]$paths <- list()
    demands[[j]]$paths[[p]] <- demands[[j]]$demand
    j <- j + 1
  }
  res$demands <- demands
  gflow <- rg.integer.min.congestion.flow.ga.graph(bestSolutions[1,],splitflow)
    
  res$gflow <- gflow
  res$val <- 1 - minEval
  res
}

### Utility function used by rg.integer.min.congestion.flow.ga()
### generates a blank chromasome with gene lengths correct
### for the demands paths
### demands - input demands list from max concurrent flow
###           paths is converted to chromasome (each a gene)
### return - chromasome each gene represents one of a possible
###          path for each demand. Gene length has to be power
###          of two - so might be more bits than paths
rg.demands.paths.as.chromasome <- function(demands) {
  chromasome <- list()
  chromasome$numgenes <- length(demands)
  chromasome$genelengths <- c()
  chromasome$pathlengths <- c()
  chromasome$string <- c()
  j <- 1
  totallength <- 0
  for(i in demands) {
    chromasome$pathlengths[j] <- length(i$paths)
    chromasome$genelengths[j] <- ceiling(log(length(i$paths),2))
    totallength <- totallength + ceiling(log(length(i$paths),2))
    j <- j + 1
  }
  chromasome$totallength <- totallength
  chromasome$demands <- demands
  chromasome
}

### Utility function used by rg.integer.min.congestion.flow.ga()
### given a populated chromasome return the paths represented
### by this chromasome
### chromosome - chromasome generated using rg.demands.paths.as.chromasome()
###              chromasome$string populated randomly by rbga.bin
### return - list each element a path (for each demand)
rg.demands.paths.from.chromasome <- function(chromasome) {
  start <- 1
  j <- 1
  paths <- list()
  for(i in names(chromasome$demands)) {
    finish <- chromasome$genelengths[j] + start - 1
    gene <- chromasome$string[start:finish]
    pathselector <- rg.binary.to.integer(gene) + 1
    ## possible problem here, once got
    ## "Error in if (pathselector <= chromasome$pathlengths[j]) { : 
    ##  missing value where TRUE/FALSE needed"

    ##cat("pathselector",pathselector,"\n")
    if ( pathselector <= chromasome$pathlengths[j] ) {
      paths[[i]] <- names(chromasome$demands[[i]]$paths)[[pathselector]]
    } else {
      paths[[i]] <- "NULL"
    }
    start <- finish + 1
    j <- j + 1
  }
  paths
}

### Utility function used by rg.demands.paths.from.chromasome()
### string - a vector of binary digits
### return - value of vector of binary digits as an integer
rg.binary.to.integer <- function(string=c(0)) {
  sum(string*2^(rev(seq(along.with=string))-1))
}

### Utility function used by rg.integer.min.congestion.flow.ga()
### given a chromasome string and the output from rg.minimum.congestion.flow()
### calculate the graph represented by the string.
### string - a chromasome$string - string is generated by rbga.bin()
### split - the output from rg.minimum.congestion.flow() which
###         supplies the demand paths and graph
### return - GraphNEL with attr weight set with integer flows for the demand
rg.integer.min.congestion.flow.ga.graph <- function(string,split) {
  chromasome <- rg.demands.paths.as.chromasome(split$demands)
  chromasome$string <- string

  paths <- rg.demands.paths.from.chromasome(chromasome)
  
  g <- rg.set.all.graph.edge.weights(split$gflow)
  
  j <- 1
  for(p in paths) {
    d <- chromasome$demands[[j]]$demand
    j <- j + 1
    pathv <- strsplit(p,"|",fixed=TRUE)
    from <- pathv[[1]][1:length(pathv[[1]])-1]
    to <- pathv[[1]][2:length(pathv[[1]])]
    
    edgeData(g,from=from,to=to,attr="weight") <-
      as.double(edgeData(g,from=from,to=to,attr="weight")) + d
  }
  g
}

rg.count.accepted.demands <- function # Count the number of accepted demands
### Note that due to numerical accuracy comparing flow == demand
### is potentially risky. This actually checks if flow >= demand - 1e-03 to
### allow for rounding errors. This limit can be set
(demands,     ##<< demands to check
 limit=1e-03  ##<< tolerance on equality of flow == demand to allow for rounding errors
 ) {
  accepted <- sapply(demands,function(i) {i$flow >= i$demand - limit})
  return(length(accepted[accepted==TRUE]))
  ### number of accepted demands
}

### Returns a graph with number of flows in each edge
rg.count.edge.flows <- function(g,demands) {
  for(i in 1:length(demands)) {
    for(j in names(demands[[i]]$paths)) {
      ##results$demands[[i]]$paths <- 1.0
      demands[[i]]$paths[[j]] <- 1.0
    }
  }
  gcount <- rg.max.concurrent.flow.graph(g,demands)
  return(gcount)
}

### Calculate the integer minimum congestion flow. It uses the max
### concurrent flow (fractional) to try various combinations and records the best
### NOTE! Probably better to use rg.max.concurrent.flow.int() instead
### if you want to find the best integer flows that meet the maximum
### capacity. This function will overbook traffic (for negative gamma and lambda)
rg.max.concurrent.flow.int.c <- function # Integer minimum congestion flow
(g, ##<< graphNEL to optimise
 demands, ##<< demands in a list ("1"=list("source","sink","demand","flow","paths")
          ##   paths is a list of paths from a fractional flow result
 e=0.1, ##<< approximation limit 0 means exact, infinitely slow, 1.0 inexact but fast
        ##   best to choose a value between 0.05-0.1
 eInternal=NULL, ##<< change e value used internally, best left NULL
 updateflow=TRUE, ##<< update the flow (slightly fast if FALSE, but not useful!)
 progress=FALSE, ##<< create a progress bar, not working in this version
 scenario=NULL, ##<< use a previous scenario to seed the paths used to attempt
 permutation="random", ##<< how the demands are chosen can be
                       ## either c(....) integers specifying demand order
                       ## or "fixed" done in fixed order
                       ## "random" done in random order
                       ## "lowest" done in lowest cost (lowest dual path) order
 deltaf=1.0, ##<< fix the update to start with a higher value to try speeding up
             ##   best left alone, but if used make it higher than 1
 linkgroupmap=NULL ##<< used to group links which are multiple wavelengths on a fibre
                   ## not working yet, leave NULL
 ) {

  ## This is a bit of a mess. RBGL only understands indexes for nodes
  ## so all the names need to be mapped to indices. This is for nodes, edges
  ## and the linkgroup map. They will be remapped at the end.
    link2linkgroup <- c(-1)
    linkgroupcap <- c(-1)
  if(!is.null(linkgroupmap)) {
    link2name <- names(edgeData(g))
    
    linkgroup2name <- unique(as.character(linkgroupmap))
    
    link2linkgroup <- rep(NA,length(link2name))
    linkgroupcap <- rep(0.0,length(linkgroup2name))
    for(n in names(linkgroupmap)) {
      i <- linkgroupmap[[n]]
      linkgroup <- which(linkgroup2name == i)
      edgepair <- strsplit(n,"\\|")[[1]]
      
      linkgroupcap[linkgroup] <- linkgroupcap[linkgroup] +
        as.double(edgeData(g,from=edgepair[[1]],to=edgepair[[2]],attr="capacity"))
    link2linkgroup[which(link2name == n)] <-
      linkgroup
    }
    link2linkgroup[is.na(link2linkgroup)] <- 0
    ## Not sure about below!
    link2linkgroup <- link2linkgroup - 1
  } 
  nodelabels <- nodes(g)
  demands <- rg.demands.relable.to.indices(demands,nodelabels)
  g <- rg.relabel(g)

  ## OK all the relabeling is done, phew
  
  ## note this is not the dual value
  calcD <- function() {
    sum(vcapacity*vlength)
  }

  if(is.null(eInternal))
    eInternal <- e

  if(is.numeric(permutation))
    permutation <- permutation - 1
  else if(permutation == "fixed")
    permutation <- 0:(length(demands)-1)
  else if(permutation == "random")
    permutation <- -1
  else if(permutation == "lowest")
    permutation <- -2

  savedemands <- demands
  vdemands <- as.double(lapply(demands,"[[","demand"))
  vcapacity <- as.double(edgeData(g,attr="capacity"))
  vlength <- rep(0.0,length(vcapacity))
  vdemandflow <- rep(0.0,length(vcapacity))
  L <- length(vlength)
  delta <- deltaf * (L / (1-e)) ^ (-1/e)
  vlength <- delta / vcapacity
  edgeMap <- adjMatrix(g)

  demands.paths <- as.integer(c())
  bestpaths <- as.integer(c())
  rdemandpaths <- list();

  if(!is.null(scenario)) {
    bestgamma=scenario$intgammalist[[1]]
    for(d in 1:length(demands)) {
      rdemandpaths[[d]] <- list()
      demands.paths <- append(demands.paths,d-1)
      demands.paths <- append(demands.paths,length(demands[[d]]$paths))
      missing <- TRUE
      
      for(p in 1:length(demands[[d]]$paths)) {
        pathi <-
          as.integer(strsplit
                     (names(demands[[d]]$paths[p]),"|",fixed=TRUE)[[1]])
        pathi <- pathi -1
        if(identical(names(demands[[d]]$paths[p]),
                     names(scenario$intdemands[[d]]$paths))) {
          cat("best path in",d,names(scenario$intdemands[[d]]$paths),", ",p,"\n")
          missing <- FALSE
          bestpaths[d] <- p
        }
        rdemandpaths[[d]][[p]] <- pathi
        demands.paths <- append(demands.paths,length(pathi))
        demands.paths <- append(demands.paths,pathi)
      }
      if(missing)
        cat("best path in",d,"is MISSING!\n")
      
    }
    print(bestpaths)
  } else {
    bestgamma <- Inf
  }
  ## code the paths as integers
  ## [demandno(e.g=1) numdempaths lengthpath1 path1[1] path1[2] ...
  ##  lengthpath2 path2[1] path2[2]....demandno(e.g=2)....]

  demandpaths <- list()
  j <- 1
  for(d in demands) {
    demandpaths[[j]] <- list()
    k <- 1
    for(p in names(d$paths)) {
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- as.integer(pv[1:length(pv)-1])
      tolist <- as.integer(pv[2:length(pv)])
      me <- rbind(fromlist,tolist)
      pv <- c()
      for(i in 1:length(fromlist)) {
        v <- as.vector(me[,i])
        em <- edgeMap[v[1],v[2]]
        pv <- append(pv,em)
      }
      demandpaths[[j]][[k]] <- pv -1
      k <- k + 1
    }
    j <- j + 1
  }


  
  m <- length(rg.edgeL(g))

  # old for ori
  # doubreq <- 2 * ceiling(1/e * log(m/(1-e),base=(1+e)))
  doubreq <- 2* ceiling(log(1.0/delta,1+eInternal))
  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = doubreq, style=3)
  } else {
    pb <- NULL
  }
  retlist <- .Call("rg_max_concurrent_flow_int_c",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   e,
                   eInternal,
                   progress,
                   pb,
                   bestgamma,
                   environment(),
                   permutation,
                   deltaf,
                   link2linkgroup,
                   linkgroupcap
                   );
  
  if( progress != FALSE) {
    close(pb)
  }

  foundbestpaths <- retlist$bestpaths + 1
  ## now need to unpack results
  demands <- savedemands
  
  for(n in 1:length(demands)) {
    flow <- 0
    for(i in 1:length(demands[[n]]$paths)) {
      flow <- flow + retlist$pathflows[[n]][[i]]
      demands[[n]]$paths[[i]] <- retlist$pathflows[[n]][[i]]
    }
    demands[[n]]$flow <- flow
  }

  gdual <- g
  gdual <- rg.set.weight(gdual,retlist$vlengths)
  scalef <- 1 / log(1/delta,base=(1+e))
  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  beta <- calcBeta(demands,gdual)
  betar <- calcBetaRestricted(demands,gdual)
  lambda=NULL
  lambda <- calcLambda(demands)
  foundratio <- beta / lambda
  ratiobound <- (1-e)^-3
    
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  if(FALSE) {
    for(d in 1:length(demandpaths)) {
      cat("demand ",d,": ")
      for(p in 1:length(demandpaths[[d]])) {
        mincap = Inf
        for(ed in demandpaths[[d]][[p]]) {
          if(vcapacity[ed+1] < mincap) {
            mincap = vcapacity[ed+1]
          }
        }
        cat(" ",retlist$pathcount[[d]][[p]])
        if(mincap < savedemands[[d]]$demand) {
          cat("u")
        } else {
          cat("o")
        }
        
        if(foundbestpaths[[d]] == p)
          cat("F")
        
        if(identical(names(demands[[d]]$paths[p]),
                     names(scenario$intdemands[[d]]$paths))) {
                                        #if(bestpaths[[d]] == p)
          cat("X")
        }
      }
      cat("\n")
    }
  }

  intdemands <- list()
  bestpaths=retlist$bestpaths+1

  for(i in names(demands)) {
    intdemands[[i]]$source <- demands[[i]]$source
    intdemands[[i]]$sink <- demands[[i]]$sink
    intdemands[[i]]$demand <- demands[[i]]$demand
    intdemands[[i]]$flow <- demands[[i]]$demand
    intdemands[[i]]$paths <- list()
    intdemands[[i]]$paths[[names(demands[[i]]$paths[bestpaths[as.integer(i)]])]] <-
      demands[[i]]$demand      
  }




  ##put back the demands labels instead of indices
  demands <- rg.demands.relable.from.indices(demands,nodelabels)
  intdemands <- rg.demands.relable.from.indices(intdemands,nodelabels)
  nodes(gflow) <- nodelabels
  nodes(gdual) <- nodelabels

  gflowint <- rg.max.concurrent.flow.graph(gdual,intdemands)
    
    ##value<< Returns a list with many elements, the most useful are
  retval <- list(demands=demands,
                 gflow=gflow,
                 gdual=gdual,
                 beta=beta,
                 betar=betar,
                 lambda=lambda,
                 phases=retlist$totalphases,
                 e=e,
                 vlength=retlist$vlengths,
                 countgamma=retlist$countgamma,
                 bestgamma=retlist$bestgamma,
                 bestpaths=bestpaths,
                 pathdiffcount=retlist$pathdiffcount,
                 phasepathdiffcount=retlist$phasepathdiffcount,
                 gammavals=retlist$gammavals,
                 betavals=retlist$betavals,
                 lambdavals=retlist$lambdavals,
                 intdemands=intdemands,     ##<< intdemands - the output demands list with the "optimal" path for each demand
                 gflowint=gflowint     ##<< gflowint - graphNEL with edge weight equal to assigned flow
)
  return(retval)
}
