require("igraph")

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
