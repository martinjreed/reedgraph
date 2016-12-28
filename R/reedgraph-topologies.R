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


## For GraphNEL add extra nodes to give single link
## fanout of order n (n is in addition to any other
## links that connect to other nodes)
rg.addFanOut <- function # Augment graphNEL with fanout
(g, ##<< graphNEL to augment
    n=2 ##<< fanout at each node
) 
{
    startindex <- length(nodes(g)) + 1
    for(d in nodes(g)) {
        if(length(edges(g)[[d]])==1)
            ## then this is already a fanout node so
            next
        ##print(d)
        ## count the number of single ended links
        count <- n - sum(lengths(edges(g)[edges(g)[[d]]])==1)
        if(count > 0)
            for(i in 1:count) {
                ##cat("adding node ",as.character(startindex),"to",d,"\n")
                g <- addNode(as.character(startindex),
                             g,edges=list(d))
                ##cat("adding edge from",as.character(d),"to",
                ##    as.character(startindex),"\n")
                g <- addEdge(from=as.character(d),to=as.character(startindex),g)
                startindex <- startindex +1
            }
    }
    ##value<< Returns a new graphNEL with augmented fanout
    return(g)
}

## generate a fully connected core graphNEL as directed graph
## 
rg.generate.connected.core <- function ##generate a fully connected core graphNEL
(n=4 ##<< the number of core nodes to interconnect
) {
    V <- as.character(1:n)
    edL1 <- vector("list",length=n)
    names(edL1) <- V
    for(i in 1:n)
        edL1[[i]] <- list(edges=V[V!=i])
    ##value<< new graphNEL as fully connected core
    g <- graphNEL(nodes=V,edgeL=edL1,edgemode='directed')
    
}



### Creates a list with Latitude and Longitude infrormation from
### Toplogy Zoo data in igraph format
### Input: g an igraph object obtained from importing from Internet
### Topology zoo (needs attributes Latitude, Longitude and label)
### result: list of all of the nodes indexed as igraph vertex index
### (as character) each node has attribute $label, $Latitude,
### $Longitude  and $population (NULL placeholder)
### WARNING: some Internet Topology zoo Lat/Long are set to zero and
### not all graphs have unique "labels" for each vertex
rg.create.graph.pop.mapping <- function(g) {
  popMapping <- list()
  indices <- as.integer(V(g))
  ##if(any(duplicated(get.vertex.attribute(g,"label")))) {
    ##cat("WARNING in rg.create.graph.pop.mapping, some labels in graph are not unique\n")
  ##}
  for(i in indices) {
    node <- list(index=i,
                 label=get.vertex.attribute(g,"label",i),
                 Longitude=get.vertex.attribute(g,"Latitude",i),
                 Latitude=get.vertex.attribute(g,"Longitude",i),
                 population=NULL)
    popMapping[[as.character(i)]] <- node
  }
  return(popMapping)
}

### Creates a list of graph mappings (see rg.create.graph.pop.mapping)
### input: graph a list of graphs each an igraph object indexed by
### name (see rg.import.multi.gml for suitable input)
### result: list of lists, each list indexed by graph name (as input)
### each list contains a list object created by
### rg.create.graph.pop.mapping
### WARNING: some Internet Topology zoo Lat/Long are set to zero and
### not all graphs have unique "labels" for each vertex
rg.create.all.graphs.pop.mappings <- function(graphs) {
  popMappings <- list()
  for(gname in names(graphs)) {
    popMappings[[gname]] <- rg.create.graph.pop.mapping(graphs[[gname]])
    
  }
  return(popMappings)
}


