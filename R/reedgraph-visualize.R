
require(igraph)
require(tcltk)

### Plot directed igraph using tcltk
### args (with defaults if not given)
### graph (an igraph object)
### width=800
### height=800
### margin=0.1 (portion of window to leave as a margin)
### vertexattr=NULL (if vector or list each element will be displayed in tooltip)
### edgeattr=NULL (as for vertexattr but displays edges)

rgplot <- function(graph,width=800,height=800,margin=0.1,vertexattr=NULL,edgeattr=NULL,layout=NULL) {


  valtocolour <- function(val,min,max) {
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
    if ( (max-min == 0) ) {
      index <- 1
    } else {
      index <- val * (length(edgecolours) -1 ) / (max - min) + 1
    }
    return(edgecolours[index])
  }
  
  if ( is.igraph(graph) ) {
    if ( is.null(layout) )
      layout <- get.graph.attribute(graph,"layout")
  } else {
    graph <- igraph.from.graphNEL(graph)
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

  edgeweights = get.edge.attribute(graph,"weight")
  minedgew = min(edgeweights)
  maxedgew = max(edgeweights)
  
  ## do edges
  for (i in E(graph)) {
    v1 <- get.edge(graph,i)[1]
    v2 <- get.edge(graph,i)[2]

    x1 <- (layout[v1+1,2] - minx) * xscale
    x2 <- (layout[v2+1,2] - minx) * xscale
    y1 <- (layout[v1+1,1] - miny) * yscale
    y2 <- (layout[v2+1,1] - miny) * yscale

    theta <- atan( (y2-y1) / (x2-x1) )

    dy <- abs ( 7 * sin(theta) )

    dx <- abs ( 7 * cos(theta) )

    x2 <- x1 + (x2-x1)/2 + sign(x1-x2) * dx / 3.0
    y2 <- y1 + (y2-y1)/2 + sign(y1-y2) * dy / 3.0

    x1 <- x1 - sign(x1-x2) * dx
    y1 <- y1 - sign(y1-y2) * dy

    text=c("E=",i)
    colour="black"
    if ( ! is.null(E(graph)[i]$weight) ) {
      text=c(text,",",E(graph)[i]$weight)
      colour=valtocolour(E(graph)[i]$weight,minedgew, maxedgew)
    }
    if ( ! is.null(edgeattr[i] ) ){
      text=c(text,",",edgeattr[i]) } 
    
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
  }


  ## do vertices
  for (i in V(graph)) {
    xpos <- (layout[i+1,2] - minx) * xscale
    
    ypos <- (layout[i+1,1] - miny) * yscale
    #tag <- tkcreate(tkcan,"oval",xpos-5,ypos-5,xpos+5,ypos+5,fill="blue",width="2")
    ## substitute is necessary to replace the value of i each time a new inner
    ## function is created. Then this has to be actually evaluated!

    nodelab <- tkcreate(tkcan,"text",xpos,ypos,text=V(graph)[i]$name)
    nodelabbk <- tkcreate(tkcan,"rectangle",tkbbox(tkcan,nodelab),fill="white")
    tkitemlower(tkcan,nodelabbk,nodelab)

    if ( ! is.null(vertexattr[i] ) ){
      text=c("V=",V(graph)[i]$name,",",vertexattr[i]) } else {
        text=c("V=",V(graph)[i]$name)
      }
    
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


