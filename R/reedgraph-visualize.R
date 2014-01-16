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

require(tcltk,quietly=TRUE)
require(igraph,quietly=TRUE)
require("graph",quietly=TRUE,warn.conflicts=FALSE)


rg.show.paths <- function(demands,numlist=NULL) {
  num <- 0
  if(is.null(numlist))
    numlist <- seq(1,length(demands))
  for(i in numlist){
    for(p in names(demands[[i]]$paths)) {
      cat(i,": ",p,": ",demands[[i]]$paths[[p]],"\n")
      num <- num + 1
    }
  }
  cat("total number of paths is ",num,"\n")
}


createLayout <- function(graph) {
    gi <- igraph.from.graphNEL(graph)
    layout <- layout.fruchterman.reingold(gi)
    return(layout)
}


rg.plot.gammas <- function(g,layout) {
  flows <- as.double(edgeData(g,attr="weight"))
  caps <- as.double(edgeData(g,attr="capacity"))
  gammas <- (caps-flows)/caps
  plot.graphNEL(g,layout=layout,edgeweights=gammas)

}

### Plot directed igraph using tcltk
### args (with defaults if not given)
### graph (an igraph object)
### width=800
### height=800
### margin=0.1 (portion of window to leave as a margin)
### edgedata vector of same length as the number of edges
###          this will be used as the edge color data

plot.graphNEL <- function(x,y,width=800,height=800,
                          margin=0.1,edgeweights=NULL,layout=NULL,...) {
  graph <- x
  if ( match("graphNEL",class(graph),nomatch=0) == 0 ) {
    print("Error in rgGraph@plot: graph is not class graphNEL")
    return(NULL)
  }
  nodelabels=nodes(graph)
  graph <- rg.relabel(graph)
  if ( is.null(layout) ) {
    gi <- igraph.from.graphNEL(graph)
    layout <- layout.fruchterman.reingold(gi)
  }
  
  edgecolours <- c("black",
                   "grey",
                   "brown",
                   "orchid",
                   "violet",
                   "purple",
                   "blue",
                   "cyan",
                   "aquamarine",
                   "seagreen3",
                   "green",
                   "yellowgreen",
                   "greenyellow",
                   "yellow",
                   "gold",
                   "orange",
                   "red")

  valtocolour <- function(val,min,max) {
    if ( (max-min == 0) ) {
      index <- 1
    } else {
      index <- floor((val - min) * (length(edgecolours) -1 ) / (max - min) + 1)
    }
    return(edgecolours[index])
  }
  
  xrange <- ( max(layout[,2]) - min(layout[,2]) ) * 
    ( 1 + 2 * margin )
  yrange <- ( max(layout[,1]) - min(layout[,1]) ) * 
    ( 1 + 2 * margin )

  xscale <- width / xrange
  yscale <- height / yrange
  
  minx <- min(layout[,2]) - margin * xrange
  miny <- min(layout[,1]) - margin * yrange
  maxx <- max(layout[,2]) + margin * xrange
  maxy <- max(layout[,1]) + margin * yrange

  tktop <- tktoplevel()
  
  tkcan <- tkcanvas(tktop,width=width,height=height)

  tkpack(tkcan)

  em <- edgeMatrix(graph)
  numedges <- length(rg.edgeL(graph))
  
  at <- names(edgeData(graph)[[1]])

  ew <- NULL
  if( is.null(edgeweights)) {
    if ( match("weight",at,nomatch=0 )) {
      ew <- as.double(edgeData(graph,attr="weight"))
    } else {
      ew <- rep(1.0, times=numedges)
    }
  } else {
    ew <- edgeweights
  }
  minedgew = min(ew)
  cat("minedgew ",minedgew,"\n")
  maxedgew = max(ew)
  cat("maxedgew ",maxedgew,"\n")
  
  
  ## draw key
  keyframe <- tkframe(tkcan)
  numedgecolours <- length(edgecolours)
  i <- minedgew
  increment <- (maxedgew - minedgew) /(length(edgecolours)-1)
  for( color in edgecolours) {
    val <- format(i,digit=2)
    fr <- tkframe(keyframe,background=color,width=50,heigh=2)
    lb <- tklabel(keyframe,text=val)
    tkgrid(fr,lb)
    i <- i + increment
  }

  tkcreate(tkcan,"window",0,0,anchor="nw",window=keyframe)
  ## do edges
  icount <- 1
  for (i in rg.edgeL(graph)) {
    v1 <- as.integer(i[1])
    v2 <- as.integer(i[2])
    
    x1 <- (layout[v1,2] - minx) * xscale
    x2 <- (layout[v2,2] - minx) * xscale
    y1 <- (layout[v1,1] - miny) * yscale
    y2 <- (layout[v2,1] - miny) * yscale
    
    theta <- atan( (y2-y1) / (x2-x1) )
    
    dy <- abs ( 7 * sin(theta) )
    
    dx <- abs ( 7 * cos(theta) )
    
    x2 <- x1 + (x2-x1)/2 + sign(x1-x2) * dx / 3.0
    y2 <- y1 + (y2-y1)/2 + sign(y1-y2) * dy / 3.0
    
    x1 <- x1 - sign(x1-x2) * dx
    y1 <- y1 - sign(y1-y2) * dy
    
    text=paste("E=",paste(i,collapse="-"),sep="")
    colour="black"
    
    
    for(n in at){
      text <- paste(text," ",n,"=",
                    format(edgeData(graph,from=i[1],to=i[2])[[1]][[n]],digit=3),
                    ", ",sep="",collapse="")
    }
    ##weight <- edgeData(graph,from=i[1],to=i[2],attr="weight")
    if( !is.null(edgeweights)) {
      text <- paste(text," inw=",format(ew[icount],digit=3),sep="",collapse="")
    }
    weight <- ew[icount]
    colour=valtocolour(as.double(weight),minedgew, maxedgew)
    
    tag <- tkcreate(tkcan,"line",x1,y1,x2,y2,fill=colour,width=2)
    tkitembind(tkcan,tag,"<Enter>",eval(substitute(function(x="%x",y="%y") {
      tkdelete(tkcan,"tooltip")
      texttip <- tkcreate(tkcan,"text",x,{as.integer(y)-10},text=text,tag="tooltip")
      whback <- tkcreate(tkcan,"rectangle",tkbbox(tkcan,texttip),fill="white",tag="tooltip")
      tkitemlower(tkcan,whback,texttip)
    })))
    tkitembind(tkcan,tag,"<Leave>",eval(substitute(function(x="%x",y="%y") {
      tkdelete(tkcan,"tooltip")
    })))
    icount <- icount +1
  }


  ## do vertices
  for (i in as.integer(nodes(graph))) {
    xpos <- (layout[i,2] - minx) * xscale
    
    ypos <- (layout[i,1] - miny) * yscale
                                        #tag <- tkcreate(tkcan,"oval",xpos-5,ypos-5,xpos+5,ypos+5,fill="blue",width="2")
    ## substitute is necessary to replace the value of i each time a new inner
    ## function is created. Then this has to be actually evaluated!

    nodelab <- tkcreate(tkcan,"text",xpos,ypos,text=nodelabels[i])
    nodelabbk <- tkcreate(tkcan,"rectangle",tkbbox(tkcan,nodelab),fill="white")
    tkitemlower(tkcan,nodelabbk,nodelab)

    text=paste("V=",as.character(i),sep="")
    
    tkitembind(tkcan,nodelab,"<Enter>",eval(substitute(function(x="%x",y="%y") {
      tkdelete(tkcan,"tooltip")
      texttip <- tkcreate(tkcan,"text",x,{as.integer(y)-10},text=text,tag="tooltip")
      whback <- tkcreate(tkcan,"rectangle",tkbbox(tkcan,texttip),fill="white",tag="tooltip")
      tkitemlower(tkcan,whback,texttip)

    })))
    tkitembind(tkcan,nodelab,"<Leave>",eval(substitute(function(x="%x",y="%y") {
      tkdelete(tkcan,"tooltip")

    })))
    
  }
  tktop
}

setMethod("plot",
          signature(x="graphNEL"),plot.graphNEL)
