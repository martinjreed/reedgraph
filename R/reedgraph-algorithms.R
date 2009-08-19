require("RBGL")

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

rg.set.demands.demand <- function(demands,val) {
  for(n in names(demands)) {
    demands[[n]]$demand <- val
  }
  demands
}

### Sets the capacity of the edges on the graph
rg.set.capacity <- function(g,val) {
  em <- edgeMatrix(g)
  if (match("capacity",names(edgeDataDefaults(G)),nomatch=0)) {
    edgeDataDefaults(g,"capacity") <- 1.0
  }
  edgeData(g,
           from=as.character(em[1,]),
           to=as.character(em[2,]),
           attr="capacity") <- val
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
###
rg.addto.weight.on.path <- function(g,path,val) {
  if(length(path) == 1) {
    edgeData(g,as.character(path[1]),as.character(path[2]),"weight") <-
      edgeData(g,as.character(path[1]),as.character(path[2]),"weight") +val
  } else {
    fromlist <- path[1:{length(path)-1}]
    tolist <- path[2:{length(path)}]
    newvals <- as.vector(pathWeights(g,as.character(path))) + val
    edgeData(g,as.character(fromlist),as.character(tolist),"weight") <- newvals
  }
  return(g)
}

### adds a value val to the weight attribute on each edge
### in path
### input graphNEL, path vector of integer, val is  the amount to add
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

  updateFlow <- function(penult,mincap,demand) {
    s <- demand$source
    t <- demand$sink
    p <- t
    t <- penult[[t]]
    tag <- paste(t,"|",p,sep="")
    if(is.null(demand$edges[[tag]])) {
      demand$edges[[tag]] <- mincap
    } else {
      demand$edges[[tag]] <-
        demand$edges[[tag]] + mincap
    }
    
    while(s !=t ) {
      p <- t
      t <- penult[[t]]
      tag <- paste(t,"|",p,sep="")
      if(is.null(demand$edges[[tag]])) {
        demand$edges[[tag]] <- mincap
      } else {
        demand$edges[[tag]] <-
          demand$edges[[tag]] + mincap
      }
      
    }
    
    demand$flow <- demand$flow + mincap
    demand
  }

  for(c in names(demands)) {
    demands[[c]]$edges <- list()
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
    demands[[ccount]] <- updateFlow(g.sp[[i$source]],i$demand,i)
    ccount <- ccount + 1
  }

  f <- as.double(edgeData(gsol,attr="weight"))
  c <- as.double(edgeData(gsol,attr="capacity"))
  lambda <- min(c/f)
  
  demands <- rg.max.concurrent.flow.rescale.demands.flow(demands,lambda)
  gsol <- rg.max.concurrent.flow.graph(gsol,demands)

  retval <-  list(demands=demands,gflow=gsol,lambda=lambda)
  
  return(retval)
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
rg.max.concurrent.flow.prescaled <- function(g,demands,e=0.1,progress=FALSE,ccode=TRUE,
                                             updateflow=TRUE) {

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
                                                   updateflow=FALSE,progress=progress)
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
                                              progress=progress,updateflow=updateflow)
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
### output: list:
###              graph: with weight set to edge flow
###              demands: each with met demand and paths


rg.fleischer.max.concurrent.flow.c <- function(g,demands,e=0.1,updateflow=TRUE,
                                               progress=FALSE) {

  calcD <- function() {
    sum(as.double(edgeData(gdual,attr="capacity"))*
        as.double(edgeData(gdual,attr="weight")))
  }


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
    pb <- tkProgressBar(title = "progress bar", min = 0,
                        max = doubreq, width = 300)
  } else {
    pb <- NULL
  }
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
                   progress
                 )
  
  
  demflow <- retlist[[2]]
  for(c in names(demands)) {
    demands[[c]]$edges <- list()
    demands[[c]]$flow <- demflow[[as.integer(c)]]
    demands[[c]]$paths <- list()
  }
  if(retlist[[3]][[1]] > 0) {
    retdemkey <- retlist[[4]]
    retdemval <- retlist[[5]]
    i <- 1
    j <- 2
    k <- 3
    for(n in seq(1:retlist[[3]][[1]])) {
      tag <- paste(retdemkey[[j]],"|",retdemkey[[k]],sep="")
      demand <- retdemkey[[i]]
      demands[[demand]]$edges[[tag]] <- retdemval[[n]]
      i <- i+3
      j <- j+3
      k <- k+3
    }
  }




  if(retlist[[3]][[2]] > 0) {
    retdemkey <- retlist[[6]]
    retdemval <- retlist[[7]]
    i <- 1
    for(n in seq(1:retlist[[3]][[2]])) {
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
                   phases=retlist[[3]][[3]],e=e,retlist=retlist)
}

## at the end this is the dual solution value
##
updateExplicitFlow <- function(penult,mincap,demand) {
  
  s <- demand$source
  t <- demand$sink
  p <- t
  t <- penult[[t]]
  tag <- paste(t,"|",p,sep="")
  
  path <- tag
  
  if(is.null(demand$edges[[tag]])) {
    demand$edges[[tag]] <- mincap
  } else {
    demand$edges[[tag]] <-
      demand$edges[[tag]] + mincap
  }
  while(s !=t ) {
    p <- t
    t <- penult[[t]]
    tag <- paste(t,"|",p,sep="")
    if(is.null(demand$edges[[tag]])) {
      demand$edges[[tag]] <- mincap
    } else {
      demand$edges[[tag]] <-
        demand$edges[[tag]] + mincap
    }
    path <- paste(t,"|",path,sep="")
  }
  
  if(is.null(demand$paths[[path]])) {
    demand$paths[[path]] <- mincap
  } else {
    demand$paths[[path]] <- demand$paths[[path]] + mincap
    
  }
  
  demand
}



rg.fleischer.max.concurrent.flow <- function(g,demands,e=0.1,updateflow=TRUE,
                                             progress=FALSE) {


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
    demands[[c]]$edges <- list()
    demands[[c]]$paths <- list()
    demands[[c]]$flow <- 0
  }

  
  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))

  if(progress != FALSE)
    pb <- tkProgressBar(title = "progress bar", min = 0,
                        max = doubreq, width = 300)
  phases <- 0
  totalphases <- 0
  D <- calcD()
  while(D < 1 ) {
    if(progress != FALSE) {

      lab = ""
      if ( is.character(progress) )
        lab = progress
      setTkProgressBar(pb,totalphases,
                       label=paste( round(totalphases/doubreq *100, 0),
                                        lab))
    }
    if(phases > doubreq) {
      cat("doubling",doubreq,totalphases,"\n");
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

rg.max.concurrent.flow.demand.graph <- function(graph,demand) {

  graph <- rg.set.all.graph.edge.weights(graph)

    for(en in names(demand$edges)) {
      from <- strsplit(en,"|",fixed=TRUE)[[1]][1]
      to <- strsplit(en,"|",fixed=TRUE)[[1]][2]
      edgeData(graph,from=from,to=to,attr="weight") <- demand$edges[[en]]
      
    }
  graph
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
    for(en in names(demands[[dn]]$edges)) {
      demands[[dn]]$edges[[en]] <- demands[[dn]]$edges[[en]] * scalef
    }
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
rg.max.concurrent.flow.graph <- function(g,demands) {
  g <- rg.set.all.graph.edge.weights(g)
  for(d in demands) {
    for(en in names(d$edges)) {
      from <- strsplit(en,"|",fixed=TRUE)[[1]][1]
      to <- strsplit(en,"|",fixed=TRUE)[[1]][2]
      edgeData(g,from=from,to=to,attr="weight") <-
        as.double(edgeData(g,from=from,to=to,attr="weight")) + d$edges[[en]]
    }
  }
  g
}

rg.mcf.find.gamma <- function(gflow,lambda=1.0) {
  flows <- as.double(edgeData(gflow,attr="weight"))/lambda
  capacities <- as.double(edgeData(gflow,attr="capacity"))

  gamma <- min ( (capacities - flows)/ capacities)
  gamma
}

rg.minumum.congestion.flow <- function(g,demands,e=0.1,prgress=FALSE) {

  res <- rg.max.concurrent.flow.prescaled(g,demands,e,progress=progress)
  res$demands <- rg.max.concurrent.flow.rescale.demands.flow(res$demands,1/res$lambda)
  res$gflow <- rg.max.concurrent.flow.graph(res$gflow,res$demands)
    
  res
}

rg.try.max <- function(g,demands,e=0.1,progress=FALSE) {

  # lets try to route gamma times capacity on each edge
  # maximise gamma
  
  res <- rg.max.concurrent.flow.prescaled(g,demands,e,progress=progress)
  resmcf <- res
  mcflambda <- res$lambda
  maxlambda <- res$beta

  if ( res$lambda < 1.0 ) {
    res$gamma <- 0
    return(res)
  }

  mingap <- (1-e)^(-3)


  # need better estimates of these, my guess is that a reasonable lower
  # lower bound on gamma is (from basic MCF)
  # min ( ( c_e - f_e / lambda )/ c_e )
  # as we can do this well, but might do better
  edgeflows <- (as.double(edgeData(res$gflow,attr="weight"))) 
  capacities <- as.double(edgeData(res$gflow,attr="capacity"))
  lbound <- min ( (capacities - edgeflows/res$lambda) / capacities)

  lbound <- lbound / mingap
  ilbound <- lbound

#  lbound <- 0.0
  
  # a maximum bound is when the smallest flow goes over the maximum capacity path
  # below is a bad estimate of this

  demval <- as.double(lapply(demands,"[[","demand"))

  ubound <- ( max(capacities) - min(demval) )/ max(capacities)

  # now try to find it by iteratively searching for best gamma

  boundgap <- ubound / lbound


  lambda <- 0.0
  resaug <- list()
  while ( boundgap > mingap ) {
    boundgap <- ubound/lbound
    gamma <- (ubound - lbound) /2 + lbound
      ## augment demands
    demaug <- demands
    ## demands are numbered in order so find where to start adding
    i <- length(demaug) + 1
    for(ed in rg.edgeL(g)) {
      ndem <- list()
      ndem$source <- ed[[1]]
      ndem$sink <- ed[[2]]
      capacity <- as.double(edgeData(g,from=ed[[1]],to=ed[[2]],attr="capacity"))
      ndem$demand <- capacity * gamma
      demaug[[as.character(i)]] <- ndem
      i <- i + 1
    }
    resaug <- rg.max.concurrent.flow.prescaled(g,demaug,e,progress=progress,
                                               updateflow=TRUE)
    
    testdem <- as.vector(resaug$demands)[seq(1:10)]
    lambda <- calcLambda(testdem)
    
    if(lambda < 1) {
      ## gamma too large
      ubound <- gamma
    } else {
      ## gamma too small
      lbound <- gamma
    }
  }

  if( lambda < 1.0 ) {
    gamma <- lbound
    ## augment demands
    demaug <- demands
    ## demands are numbered in order so find where to start adding
    i <- length(demaug) + 1
    for(ed in rg.edgeL(g)) {
      ndem <- list()
      ndem$source <- ed[[1]]
      ndem$sink <- ed[[2]]
      capacity <- as.double(edgeData(g,from=ed[[1]],to=ed[[2]],attr="capacity"))
      ndem$demand <- capacity * gamma
      demaug[[as.character(i)]] <- ndem
      i <- i + 1
    }
    resaug <- rg.max.concurrent.flow.prescaled(g,demaug,e,progress=progress,
                                               updateflow=TRUE)
    
    testdem <- as.vector(resaug$demands)[seq(1:10)]
    lambda <- calcLambda(testdem)

  }
  
  resaug$gamma <- gamma
  resaug$ilbound <- ilbound
  # get rid of the artificially generated demands
  resaug$demands <- resaug$demands[1:length(demands)]
  resaug$mcflambda <- mcflambda
  resaug$resmcf <- resmcf
  resaug
}


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


# test the rg.try.max function across many graphs/demands seems to
# show that res$ilbound <= res$gamma <= res$ilbound*(1+w) in other
# words we do not need to run rg.try.max but simply run MCF and then
# set flows to the demand! Needs prooving.
rg.test.idea <- function() {


  run <- 1

  while(run < 50) {
    n <- floor(runif(1,10,30))
    
    nd <- floor(runif(1,10,100))
    
    e <- runif(1,0.005,0.05)
    
    g <- rg.generate.random.graph(n)
    edgeDataDefaults(g,"capacity") <- 1.0
    fromlist <- as.character(edgeMatrix(g)[1,])
    m <- length(fromlist)
    tolist <- as.character(edgeMatrix(g)[2,])
    
    
    edgeData(g,from=fromlist,to=tolist,attr="capacity") <- runif(m,3,20)
    dem <- rg.gen.demands(g,nd,runif(nd,0.5,5))
    
                                        # find the demension of the demand flow (use lambda to work this out
    res <- rg.max.concurrent.flow.prescaled(g,dem,e=e,updateflow=FALSE)
    
    for(fi in seq(1:5) ) {
      scalef <- fi/5 
      
      dems <- rg.rescale.demands(dem,scalef*res$lambda * 0.9)
      res <- rg.try.max(g,dems,e=e)
      cat("lambda =",res$lambda,"\n")
      cat("gamma =",res$gamma,"\n")
      cat("ilbound =",res$ilbound,"\n")
      
      filename <- paste("/home/mjreed/tmp/run",run,"n",n,"nd",nd,"fi",fi,"e",e,".robj",sep="")
      run <- run + 1
      save(res,file=filename)
    }
    
  }  
    
    res
}

rg.test.idea.stats <- function() {
  files <- Sys.glob("~/tmp/run[0-9]*")
  files <- files[order(as.numeric(gsub("^.*run([0-9]*).*$","\\1",files,files)))]
  for(file in files) {
    load(file)
    print(file)
    e <- res$e
    mingap <- (1-e)^(-3)
#    tmpres <- rg.max.concurrent.flow.prescaled(res$gflow,res$demands,e)
#    cat(res$gamma,res$ilbound,res$mingap,"\n")
#    cat(res$gamma >= res$ilbound && res$gamma <= res$ilbound*mingap)
    if( res$gamma !=0 && is.numeric(res$gamma) && is.numeric(res$ilbound)) {
      if( res$gamma >= res$ilbound && res$gamma <= res$ilbound*mingap) {
        ##        cat(res$gamma/res$ilbound, mingap, res$gamma,res$ilbound,res$ilbound*mingap,"\n")
 #       cat(res$gamma/res$ilbound, mingap, res$gamma,res$ilbound,res$ilbound*mingap,"\n")
        gamma <- rg.mcf.find.gamma(res$resmcf$gflow,res$resmcf$lambda)
        cat(res$gamma,gamma,1- 1/res$mcflambda,"\n")
        if(res$gamma > gamma )
          print("Warning congecture false! res$gamma > gamma  ")
        if(gamma < 1- 1/res$mcflambda ) {
          print("Warning lambda is not a lower limit")
        }
        caplimit <- FALSE
        for(d in names(res$demands)) {
          for(p in names(res$demands[[d]]$paths)) {
            pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
            fromlist <- pv[1:length(pv)-1]
            tolist <- pv[2:length(pv)]
            mincap <- min(as.double(edgeData(res$gflow,from=fromlist,to=tolist,att="capacity")))
            if(res$demands[[d]]$demand > mincap) {
  ##            cat(p,":",res$demands[[d]]$demand,mincap,"\n")
              caplimit <- TRUE
            }

          }
        }
        if(caplimit)
          print("Capacity limited!")
#        cat("  ",1/(1-res$gamma),res$mcflambda,"\n")
      } else {
        cat("FALSE",res$gamma,res$ilbound,res$ilbound*mingap,"\n")
      }
    }
  }
}

# idea obtained in flash may not work! Try first then try to proove
# why it may work

rg.integer.min.congestion.flow <- function(g,demands,e=0.1) {

  sres <- rg.max.concurrent.flow.prescaled(g,demands,e)
  sres$mcfdemands <- sres$demands
  sres$mcfgflow <- sres$gflow
  m <- length(rg.edgeL(g))

  nd <- length(demands)
  delta <- (m / (1-e)) ^ (-1/e)

  capacities <- as.double(edgeData(res$gflow,attr="capacity"))

  # or should this be the gdual lengths?
  weights <- edgeData(sres$gflow,attr="weight")
  for(w in seq(1:length(weights))) 
    if(weights[w] < delta/capacities[w])
      weights[w] <- nd * delta/capacities[w] 

  for(d in names(sres$demands)) {
    demands[[d]]$edges <- list()
    demands[[d]]$paths <- list()
    demand <- sres$demands[[d]]
    
    dgraph <- rg.max.concurrent.flow.demand.graph(sres$gflow,
                                                  demand)
    dflow <- edgeData(dgraph,attr="weight")

    for(w in seq(1:length(dflow))) 
      if(dflow[w] < delta/capacities[w])
        dflow[w] <- delta/capacities[w] 


    dweights <- as.double(weights)/as.double(dflow)
    
    sp <- dijkstra.sp(dgraph,demand$source,dweights)

    sres$demands[[d]] <- updateExplicitFlow(sp$penult,demand$demand,demands[[d]])
  }

  
  sres$gflow <- rg.max.concurrent.flow.graph(sres$gflow,sres$demands)

  sres
}


rg.find.demand.paths <- function(demands,source=NULL,sink=NULL,edge=FALSE) {

  res <- list()
  if(edge) {
    mtext <- paste("(\\||^)",source,"\\|",sink,"(\\||$)",sep="")
  } else {
    mtext <- paste("(\\||^)",source,"(\\|.*\\|)|(\\|)",sink,"(\\||$)",sep="")
  }
  
  
  for(n in names(demands)) {
    for(p in names(demands[[n]]$paths)) {
#      cat("looking in demand ",n,p,"\n")
      if(regexpr(mtext,p) != -1) {
#        print("found")
        if( ! (n %in% names(res)) )
          {
            res[[n]] <- list()
            res[[n]]$source <- demands[[n]]$source
            res[[n]]$sink <- demands[[n]]$sink
            res[[n]]$demand <- demands[[n]]$demand
            res[[n]]$paths <- list()
          }
        res[[n]]$paths[[p]] <- demands[[n]]$paths[[p]]
      } 
    }
  }
  res
}
    
rg.integer.min.congestion.flow.analyse <- function(res) {

  mcfdemands <- rg.max.concurrent.flow.rescale.demands.flow(res$mcfdemands,1/res$lambda)

  for(d in mcfdemands) {
    for(p in mcfdemands$paths) {
      
    }
  }

}
