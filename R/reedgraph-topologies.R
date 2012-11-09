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


### Demo for students to simulate a "dice" loaded with exponential
### like distribution.
### input: n - number of throws to return
### rate - the parameter in the (negative expontential)
### returns vector of n throw results
loadeddice <- function(n,rate=0.7) {
  x <- floor(rexp(n,rate))+1
  x <- x[ 1 <= x & x <=6]
  return(x)
}

### Demo for students to simulate a "dice"
### input: n - number of throws to return
### returns vector of n throw results
dice <- function(n) {
  x <- floor(runif(n,min=1,max=7))
  x <- x[ 1 <= x & x <=6]
  return(x)
}
