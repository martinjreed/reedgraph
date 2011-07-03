require("igraph",quietly=TRUE)
require("graph",quietly=TRUE,warn.conflicts=FALSE)
require("RBGL",quietly=TRUE,warn.conflicts=FALSE)
require("genalg",quietly=TRUE)
### return all the edges as a list of lists(head,tail)
### in graphNEL
### result is list of edges
### edgeMatrix might be better
rg.edgeL <- function(g) {
  ed <- list()
  for(i in nodes(g)) {
    for(el in edges(g)[[i]]) {
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

### Calculates the number of edge disjoint paths
### between every pair of nodes
### Returns a matrix where element [u,v] is the
### number of paths between u and v where
### 0 means disconnected (or to itself)
### 1 for one path etc
### NOTE this is an R stub calling C++ code
### NOTE THIS IS INCORRECT!
rg.num.edge.disjoint.paths.c <- function(g) {
  em <- edgeMatrix(g)
  nN <- nodes(g)
  nv <- length(nN)

  ne <- ncol(em)
  eW <- unlist(edgeWeights(g))

  count <- .Call("rg_num_edge_disjoint_paths_c",
     as.integer(nv), as.integer(ne),
     as.integer(em-1), as.double(eW))
  count <- matrix(count,nrow=nv,ncol=nv)
  count
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

### generate demands

### num: number of demands to generate, it will be this many across a
###      uniformly random selection of node pairs (i,j) where i != j
###      If NULL then it will be between every pair

### val: value of demand. Accepts a vector and iterates through in
###      order and will continue back to beggining if necessary
###      use runif(num,min,max) if uniform random required
### returns a list of list("number",source=,sink=,demand=) where
### each demand is identified by a number in the order of the demand
### therefore probably not best as mutuable unless care is taken
### that number and indexing is used correctly
rg.gen.demands <- function(g,num=NULL,val=1.0) {
  comm <- list()
  n=1;
  j <- 1
  vallen <- length(val)
  if( is.null(num) ) { # all pairs uniform demands
    for(i in nodes(g)) {
      for(j in nodes(g)) {
        if ( i != j ) {
          comm[[as.character(n)]] <- list(source=i,sink=j,demand=val[j])
          j <- j + 1
          if(j > vallen) j <- 1
          n <- n+1
        }
      }
    }
  } else { # random selection of nodes
    numnodes <- length(nodes(g))
    for(i in seq(1:num)) {
      fromto <- floor(runif(2,min=1,max=numnodes+1))
      while(fromto[1] == fromto[2])
        fromto <- floor(runif(2,min=1,max=numnodes+1))
      comm[[as.character(n)]] <- list(
                                       source=
                                       as.character(fromto[1]),
                                       sink=
                                       as.character(fromto[2])
                                       ,demand=val[j])
      j <- j + 1
      if(j > vallen) j <- 1
      n <- n+1
  }
}
comm
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

### Sets the capacity of the edges on the graph
### val - scalar or vector length of edges
rg.set.capacity <- function(g,val) {
  em <- edgeMatrix(g)
  if (match("capacity",names(edgeDataDefaults(g)),nomatch=1)) {
    edgeDataDefaults(g,"capacity") <- 1.0
  }
  edgeData(g,
           from=as.character(em[1,]),
           to=as.character(em[2,]),
           attr="capacity") <- val
  g
}

### Sets the weight of the edges on the graph
### val - scalar or vector length of edges
rg.set.weight <- function(g,val) {
  em <- edgeMatrix(g)
  if (match("weight",names(edgeDataDefaults(g)),nomatch=1)) {
    edgeDataDefaults(g,"weight") <- 1.0
  }
  edgeData(g,
           from=as.character(em[1,]),
           to=as.character(em[2,]),
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
    if(length(path) == 1) {
      edgeData(g,as.character(path[1]),as.character(path[2]),attr) <-
      edgeData(g,as.character(path[1]),as.character(path[2]),attr) +val
  } else {
    
    fromlist <- path[1:{length(path)-1}]
    tolist <- path[2:{length(path)}]
    newvals <- as.double(edgeData(g,as.character(fromlist),as.character(tolist),attr)) + val
    edgeData(g,as.character(fromlist),as.character(tolist),attr) <- newvals
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




### Compute the multicommodity flow in a graph using naive shortest
### path g graph of type graphNEL with the edge weight variable used
### for the shortest path and capacities ignored (ie it may overbook
### an edge)
### demands: a list each element is a list [source, sink, demand]
rg.sp.max.concurrent.flow <- function(g,demands) {


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
  c <- as.double(edgeData(gsol,attr="capacity"))
  lambda <- min(c/f)

  gamma <- min((c-f)/c)
    
  #demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,lambda)
    
  #gsol <- rg.max.concurrent.flow.graph(gsol,demands)
  
  retval <-  list(demands=demands,gflow=gsol,lambda=lambda,gamma=gamma)
    
  return(retval)
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
rg.max.concurrent.flow.prescaled <- function(g,demands,e=0.1,progress=FALSE,ccode=TRUE,updateflow=TRUE,permutation="fixed") {

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
                                                   updateflow=FALSE,progress=progress,permutation)
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
                                              progress=progress,updateflow=updateflow,permutation)
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


rg.fleischer.max.concurrent.flow.c <- function(g,demands,e=0.1,updateflow=TRUE,progress=FALSE,permutation="lowest") {


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
  retlist <- .Call("rg_fleischer_max_concurrent_flow_c",
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
                   permutation
                 )
  
  
  demflow <- retlist[[2]]
  for(c in names(demands)) {
    demands[[c]]$flow <- demflow[[as.integer(c)]]
    demands[[c]]$paths <- list()
  }


  if(retlist[[3]][[1]] > 0) {
    retdemkey <- retlist[[4]]
    retdemval <- retlist[[5]]
    i <- 1
    for(n in seq(1:retlist[[3]][[1]])) {
      demand <- retdemkey[[i]]
      i <- i+1
      path <- retdemkey[[i]]
      i <- i+1
      demands[[demand]]$paths[[path]] <- retdemval[[n]]
    }
  }
  delta <- (ne / (1-e)) ^ (-1/e)

  scalef <- 1 / log(1/delta,base=(1+e))

  
  gdual <- g

  fromlist <- edgeMatrix(g)[1,]
  tolist <- edgeMatrix(g)[2,]

  edgeData(gdual,from=as.character(fromlist),to=as.character(tolist),attr="weight") <- retlist[[1]]

  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,scalef)
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)

  beta <- calcBeta(demands,gdual)
  
  lambda=NULL
  if(updateflow) {
    lambda <- calcLambda(demands)
    foundratio <- beta / lambda
    ratiobound <- (1-e)^-3
  }

  if( progress != FALSE) {
    close(pb)
  }

  retlist2 <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda,
                   phases=retlist[[3]][[2]],e=e)
}

### Utility function used by rg.fleischer.max.concurrent.flow()
### penult - vetor from Dijkstra.Sp
### mincap - mincapacity to set flow to this value
### demand - the specific demand to update
updateExplicitFlow <- function(penult,mincap,demand) {
  
  s <- demand$source
  t <- demand$sink
  p <- t
  t <- penult[[t]]
  tag <- paste(t,"|",p,sep="")
  
  path <- tag
  
  while(s !=t ) {
    p <- t
    t <- penult[[t]]
    path <- paste(t,"|",path,sep="")
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
           as.character(fromedges),
           as.character(toedges),
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
        p <- extractPath(c$source,c$sink,sp$penult)
        caponpath <- as.double(edgeData(gdual,
                                        from=as.character(p[1:length(p)-1]),
                                        to=as.character(p[2:length(p)]),
                                        attr="capacity"))

        mincap <- min(demand,caponpath)

        demand <- demand - mincap
        
        lengths <- as.double(edgeData(gdual,
                                      from=as.character(p[1:length(p)-1]),
                                      to=as.character(p[2:length(p)]),
                                      attr="weight"))
        lengths <- lengths * (1 + (e*mincap) / caponpath)
        edgeData(gdual,from=as.character(p[1:length(p)-1]),
                 to=as.character(p[2:length(p)]),
                 attr="weight") <- lengths
        if(updateflow)
          c <- updateExplicitFlow(sp$penult,mincap,c)
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

rg.max.concurrent.flow.rescale.demands.flow <- function(demands,scalef) {
  for(dn in names(demands)) {
    demands[[dn]]$flow <- demands[[dn]]$flow * scalef
    for(en in names(demands[[dn]]$paths)) {
      demands[[dn]]$paths[[en]] <- demands[[dn]]$paths[[en]] * scalef
    }
  }
  demands
}

rg.rescale.demands <- function(demands,scalef) {

  for(dn in names(demands)) {
    demands[[dn]]$demand <- demands[[dn]]$demand * scalef
  }
  demands
}

rg.set.all.graph.edge.weights <- function(g,val=0.0)  {
  fromedges <- edgeMatrix(g)[1,]
  toedges <- edgeMatrix(g)[2,]

  edgeData(g,
           as.character(fromedges),
           as.character(toedges),
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
      weights <- as.double(edgeData(g,from=fromlist,to=tolist,att="weight"))
      weights <- weights + d$paths[[p]]
      edgeData(g,from=fromlist,to=tolist,att="weight") <- weights

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

rg.minimum.congestion.flow <- function(g,demands,e=0.1,progress=FALSE,permutation="random") {

  res <- rg.max.concurrent.flow.prescaled(g,demands,e,progress=progress,ccode=TRUE,permutation=permutation)
  res$demands <- rg.max.concurrent.flow.rescale.demands.flow(res$demands,1/res$lambda)
  
  res$gflow <- rg.max.concurrent.flow.graph(res$gflow,res$demands)
  
  res$gamma <- rg.mcf.find.gamma(res$gflow)
  
  res
}

### Generate a random graph with a range of node degree
### uses igraph routines
### n - number of nodes
### mindeg - minimum node degree default = 3
### maxdeg - maximum node degree default = 7
### return GraphNEL
rg.generate.random.graph <- function(n=10,mindeg=3,maxdeg=7) {
  degrees <- floor(runif(n,mindeg,maxdeg+1))

  ## note sum(degrees) has to be even, add one to smallest
  ## node degree of first node found.
  if(sum(degrees) %% 2 != 0) {
    minval <- min(degrees)
    pos <- match(minval,degrees)
    degrees[pos] <- degrees[pos] +1
  }

  gi <- degree.sequence.game(degrees,method="vl")
  gi <- as.directed(gi)
  G <- rg.relabel(igraph.to.graphNEL(gi))
  G
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

