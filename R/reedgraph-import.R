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

### read gml files such as internet-toplogy-zoo
### input
### path directory to search for files, (default is harcoded to a place on my system)
### Make sure path is not relative to "~" in Mac as it fails!
### returns list of names as keys each element is a graph in igraph format
rg.import.multi.gml <- function(path="~/main/R/topologies/internet-toplogy-zoo/") {
  files <- list.files(path,pattern="\\.gml",recursive=TRUE)

  graphs <- list()
  for(file in files) {
    ##cat("importing ",file,"\n")
    graph.name <- strsplit(file,"\\.")[[1]][1]
    ##cat("created ",graph.name,"\n")
    g <- read.graph(paste(path,"/",file,sep=""),format="gml")
    graphs[[graph.name]] <- g
  }
  return(graphs)
}

### read rocketfuel files
### result igraph
rg.import.rocketfuel <- function(file) {
  myfile <- file(file, "r")
  lines <- readLines(myfile)
  closeAllConnections()
  edgem <- matrix(nrow=0,ncol=2)
  nodecount <- 0
  nodes <- c()
  nodesbb <- c()

  ## need two passes to find the nodes
  ## find only nodes in lines with bb

  ## find only nodes with more than degree one
  for(line in lines) {
    ##cat(line,"\n")
    if(grepl("^[0-9]+\\s+@.*\\s\\+\\sbb.*<.*?>",line)) {
      # first get the node from the first parameter
      node <- paste("n",sub("^([0-9]+)\\s.*","\\1",line),sep="")
      # convergt output of regexp to matrix (pos, length)
      match <- lapply(gregexpr("<.*?>",line), function(x) cbind(x, attr(x, "match.length")))[[1]]
      # get the text in another matrix
      connected <- apply(match,1,function(x) return(c(node,paste("n",substr(line,x[1]+1,x[1]+x[2]-2),sep=""))))
      connected <- t(connected)
      # only add nodes if there is more than one link
      if(length(connected) > 2) {
        nodesbb <- c(nodesbb,node)
       }
    }
  }
  ##print(nodesbb)

  for(line in lines) {
    ##cat(line,"\n")
    ## create matrix of matchpos,length
    if(grepl("^[0-9]+\\s+@.*\\s\\+\\sbb.*<.*?>",line)) {
      node <- paste("n",sub("^([0-9]+)\\s.*","\\1",line),sep="")
      match <- lapply(gregexpr("<.*?>",line), function(x) cbind(x, attr(x, "match.length")))[[1]]
      connected <- apply(match,1,function(x) return(c(node,paste("n",substr(line,x[1]+1,x[1]+x[2]-2),sep=""))))
      connected <- t(connected)
      ##print(connected)
      bbcount <- 0
      for(i in 1:(length(connected)/2)) {
        if(match(connected[i,2],nodesbb,nomatch=0) != 0) {
          bbcount <- bbcount + 1
        }
        
      }
      ##print(bbcount)
      if(bbcount > 1) {
        nodes <- c(nodes,node)
        nodecount <- nodecount +1
        ##print(nodes)
      }
    }
  }

  ##print("nodes")
  ##print(nodes)
  for(line in lines) {
    ##cat(line,"\n")
    ## create matrix of matchpos,length
    if(grepl("^[0-9]+\\s+@.*\\s\\+\\sbb.*<.*?>",line)) {
      node <- sub("^([0-9]+)\\s.*","\\1",line)
      match <- lapply(gregexpr("<.*?>",line), function(x) cbind(x, attr(x, "match.length")))[[1]]
      connected <- apply(match,1,function(x) return(c(paste("n",node,sep=""),paste("n",substr(line,x[1]+1,x[1]+x[2]-2),sep=""))))
      connected <- t(connected)
      ##print(connected)
      bbcount <- 0
      for(i in 1:(length(connected)/2)) {
        if(match(connected[i,2],nodes,nomatch=0) != 0) {
          bbcount <- bbcount + 1
        }
        
      }
      ##print(bbcount)
      for(i in 1:(length(connected)/2)) {
        if(match(connected[i,2],nodes,nomatch=0) != 0 &
           match(connected[i,1],nodes,nomatch=0) != 0 &
           bbcount > 1) {
          ##print(i)
          ##print(connected[i,2])
          edgem <- rbind(edgem,connected[i,])
        } 
        
      }

      
    }
  }

  g <- graph.edgelist(edgem)

  ## there still might be some nodes with one out or in edge so remove them
  
  vertices <- c()
  tmp <- get.adjedgelist(g,mode="in")
  for(i in 1:length(V(g)) ) {
    if(length(tmp[[i]]) <= 1) {
      vertices <- c(vertices,i-1)
    }
  }
  g <- delete.vertices(g,vertices)

  vertices <- c()
  tmp <- get.adjedgelist(g,mode="out")
  for(i in 1:length(V(g)) ) {
    if(length(tmp[[i]]) <= 1) {
      vertices <- c(vertices,i-1)
    }
  }
  g <- delete.vertices(g,vertices)

  g <- igraph.to.graphNEL(g)
  g <- rg.relabel(g)
}


### read brite format file
### result igraph
### Read in the standard brite format see:
### www.cs.bu.edu/brite/user_manual/node29.html
### create a graph "g" which is an igraph class object

rgimportbrite <- function(file) {
  myfile <- file(file, "r")
  myskip <- 0

  ## read the file until it says "Nodes:" then stop
  line <- readLines(myfile,1)
  while ( length(grep("Nodes:", line)) == 0 ) {
    line <- readLines(myfile,1)
    myskip <- myskip + 1
  }

  ## get from the line read how many nodes we have
  numnodes <- as.integer(sub("Nodes: \\((.*)\\).*","\\1",line))
  ## read in the whole Nodes spec as a table (list)
  nodes <- read.table(myfile, sep = " ", nrows=numnodes)

  ## now do the same with edges
  while ( length(grep("Edges:", line)) == 0 ) {
    line <- readLines(myfile,1)
    myskip <- myskip + 1
  }

  numedges <- as.integer(sub("Edges: \\((.*)\\).*","\\1",line))
  edges <- read.table(myfile, sep = " ", nrows=numedges)
  ##now we have edges and nodes

  ## done with reading
  close(myfile)

  ## create an empty vector
  edgev <- c()

  ## need this to go 1,2 then 3,4 then 5,6 in each loop of i
  j <- 1

  gcoords <- matrix(,numnodes,2)

  for(i in 1:numnodes) {
    gcoords[i,1] <- nodes[[2]][i]
    gcoords[i,2] <- nodes[[3]][i]

  }
  ##get the edge from the file and put in vector edgev as pairs
  for(i in 1:numedges) {
    ## from
    edgev[j] <- edges[[2]][i]
    j <- j+1
    ## to
    edgev[j] <- edges[[3]][i]
    j <- j+1
  }

  ## this creates a graph in igraph class taking the edge vector
  ## as input
  g <- graph(edgev)

  ## this assigns the layout to be actual graph parameters

  g <- set.graph.attribute(g,"layout",gcoords)

  ## now assign the "weight" attribute as 1.0

  E(g)$weight <- 1.0

  g
}

rg.write.graphNEL <- function(g,layout=NULL,fileprefix=NULL) {
  gi <- igraph.from.graphNEL(g)
  gfile <- paste(fileprefix,".graph.gml",sep="")
  lfile <- paste(fileprefix,".layout.rtb",sep="")
  write.graph(gi,gfile,format="graphml")
  write.table(layout,lfile)
}

