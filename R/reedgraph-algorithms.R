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
  
### generates a quick and dirty commodity list
### of unit demand between each pair of nodes
rggencomm <- function(g) {
  comm <- list()
  for(i in nodes(g)) {
      for(j in nodes(g)) {
        if ( i != j )
          comm <- c(comm,list(list(head=i,tail=j,demand=1.0)))
      }
  }
  comm
}

### Sets the capacity of the edges on the graph
rg.set.capacity <- function(g,val) {
  em <- edgeMatrix(g)
  if (match("capacity",names(edgeDataDefaults(G)),nomatch=0) == 0) {
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
  nodes(myg) <- as.character(seq(1:length(nodes(G))))
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



### Compute the multicommodity flow in a graph using naive shortest path
### g graph of type graphNEL with the edge weight variable being capacity
### demands a matrix each row is a commodity vector [source, sink, demand]
### !!! At the moment this just computes a niave flow !!! ie the
### flows all follow the shortest path - very much indevelopment!
rg.sp.multicomm.flow <- function(g,demands) {
  g <- rg.relabel(g)

  gdual <- g
  edgeDataDefaults(gdual,"weight") <- 1

  gprim <- g
  edgeDataDefaults(gprim,"weight") <- 0

  g.sp <- list()
  
  for(i in demands) {
    ## calculate dijkstra.sp if we do not have one for this vertex
    if(is.null(g.sp[[i$head]])) g.sp[[i$head]] <- dijkstra.sp(gdual,start=i$head)$penult
    path <- extractPath(i$head,i$tail,g.sp[[i$head]])
    gprim <- rg.addto.weight.on.path(gprim,path,i$demand)
  }
  return(gprim)
}

rg.sp.multicomm.flow <- function(g,demands) {
  

}

