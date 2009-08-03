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
  
### generates a quick and dirty demand list
### of unit demand between each pair of nodes
### each element is [source, sink, demand]
rg.gen.demands <- function(g,num=NULL,val=1.0) {
  comm <- list()
  n=0;
  if( is.null(num) ) { # all pairs uniform demands
    for(i in nodes(g)) {
      for(j in nodes(g)) {
        if ( i != j ) {
          comm[[as.character(n)]] <- list(source=i,sink=j,demand=val)
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
                                       ,demand=val)
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
rg.sp.multicomm.flow <- function(g,demands) {

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
  edgeDataDefaults(gsol,"weight") <- 0

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
  
  demands <- rg.max.concurrent.flow.rescale.demands(demands,lambda)
  gsol <- rg.max.concurrent.flow.graph(gsol,demands)
  cat("lambda=",lambda,"\n",sep="")

  retval <-  list(demands=demands,gflow=gsol,lambda=lambda)
  
  return(retval)
}

### incomplete!
rg.max.concurrent.flow.prescaled <- function(g,demands,e=0.1) {
  estimate <- rg.sp.multicomm.flow(g,demands)
  

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

rg.max.concurrent.flow <- function(g,demands,e=0.1,updateflow=TRUE) {

  savedemands <- demands
  ## note this is not the dual value
  calcD <- function() {
    sum(as.double(edgeData(gdual,attr="capacity"))*
        as.double(edgeData(gdual,attr="weight")))
  }

  ## at the end this is the dual solution value
  ##
  calcBeta <- function(demands,gdual) {
    Alpha <- 0

    for(demand in demands) {
      sp <- dijkstra.sp(gdual,demand$source)
      Alpha <- Alpha + demand$demand * sp$distances[[demand$sink]]
    }
    Beta <- calcD() / Alpha
    Beta
  }

  ## at the end this is the primal solution value
  calcLambda <- function(demands) {
    d <- as.double(lapply(demands,"[[","demand"))
    f <- as.double(lapply(demands,"[[","flow"))
    lambda <- min(f/d)
    lambda
  }
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
    demands[[c]]$flow <- 0
  }
  print(calcD())

  doubreq <- 2/e * log(m/(1-e),base=(1+e))
  phases <- 0
  totalphases <- 0
  while(calcD() < 1 ) {
    print(calcD())

    if(phases > doubreq) {
      print("Doubling required")
      demands <- doubleDemands(demands)
      phases <- 0
    }

    ccount <- 1
    for(c in demands) {
      demand <- c$demand

      while( calcD() < 1 && demand > 0 ) {
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
          c <- updateFlow(sp$penult,mincap,c)

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

  print("finished")
  cat("D=",calcD(),"\n",sep="")
  cat("in",totalphases,"phases\n")

  scalef <- 1 / log(1/delta,base=(1+e))
  demands <- rg.max.concurrent.flow.rescale.demands(demands,scalef)

  beta <- calcBeta(demands,gdual)
  cat("beta=",beta,"\n",sep="")

  lambda=NULL
  if(updateflow) {
    lambda <- calcLambda(demands)
    cat("lambda=",lambda,"\n",sep="")
      foundratio <- beta / lambda
    ratiobound <- (1-e)^-3
    cat("Ratio dual/primal=",foundratio,
        " Ratio bound=",ratiobound,"\n",sep="")
    
  }
  gflow <- rg.max.concurrent.flow.graph(gdual,demands)
  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,lambda=lambda)
  retval
}


rg.max.concurrent.flow.rescale.demands <- function(demands,scalef) {

  for(dn in names(demands)) {
    demands[[dn]]$flow <- demands[[dn]]$flow * scalef
    for(en in names(demands[[dn]]$edges)) {
      demands[[dn]]$edges[[en]] <- demands[[dn]]$edges[[en]] * scalef
    }
  }
  demands
}

rg.max.concurrent.flow.graph <- function(g,demands) {
  edgeDataDefaults(g,"weight") <- 0
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
