library(doParallel)
registerDoParallel(cores=(detectCores()-1))

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

rg.generate.documentation <- function() {
    package.skeleton.dx(pkgdir="/Users/mjreed/main/R/reedgraph",
                        excludePattern=paste0("(reedgraph-scratch.R)|",
                            "(reedgraph-altbloom.R)|",
                            "(reedgraph-bloom.R)|",
                            "(ddlpn.R)"
                            ))
}

### NOTE THIS IS USING INLINEDOCS
### EXAMPLE
rg.inlinedocs.example <- function # Example function
### Here is some maths in the longer description \eqn{n} 
### which probably has more lines
##references<< Abelson, Hal; Jerry Sussman, and Julie
##Sussman. Structure and Interpretation of Computer
##Programs. Cambridge: MIT Press, 1984.
(n,    ##<< description of vairiables 
 times
### here is a different way: Number of Fermat tests to perform. More
### tests are more likely to give accurate results.
 ){
  if(times==0)TRUE
  ##seealso<< \code{\link{fermat.test}}
  else if(fermat.test(n)) is.pseudoprime(n,times-1)
  else FALSE
### the return value logical TRUE if n is probably prime.
}


### count number of packets to get one PPM packet from each router
rg.ppm.count.required <- function(N,p) {

  all <- FALSE
  count <- 0
  received_marks <- rep(FALSE,N)
  while(!all) {
    count <- count + 1
    mark <- NULL
    for(i in 1:N) {
      if(runif(1) < p) {
        mark <- i
      }
    }
    if(!is.null(mark)) {
      received_marks[mark] <- TRUE
    }
    if(sum(received_marks) == N) {
      break
    }
  }

  return(count)
}

rg.ppm.count.required.averaged <- function(N,p,averages) {

  sum <- 0
  for(i in 1:averages) {
    sum <- rg.ppm.count.required(N,p) + sum
  }
  return(sum/averages)
}

rg.ppm.savage <- function(N,p) {
  e <- (log(N)) / ( p*(1-p)^(N-1) )
  return(e)
}

rg.ppm.reed <- function(N,p) {
  sum <- 1
  for(i in (N-1):0) {
    print(i)
    sum <- sum  + 1/ ( p * (1-p)^i )
  }
  
  return(sum)
}

rg.ppm.test <- function(N,p) {

  prob <- p * (1-p) * (1-p)^2 * (1-p * (1-p))

  return(prob)
}
### Generate graph to simulate attack graph
### Warning this only produces a tree at the moment
rg.gen.attack.graph <- function(depth=4,mindeg=2,maxdeg=5,lower=0.3) {

  nodes <- list()
  alledges <- list()
  nodes[[1]] <- c("1")
  allnodes <- c(nodes[[1]])
  ## for each level
  for(i in 2:depth) {
    cat("creating nodes level ",i,"\n")
    count <- 1
    nodes[[i]] <- c()


    ## for each node at higher level
    cat("nodes[i-1]length=",length(nodes[[i-1]]),"\n")
    first <- TRUE
    for(j in 1:length(nodes[[i-1]])) {
      edges <- c()


      num <- floor(runif(1,mindeg,maxdeg+1))
      edges <- c(paste(i,".",as.character(count),sep=""))
      count <- count +1
      for(k in 2:num) {
        edges <- append(edges,paste(i,".",as.character(count),sep=""))
        count <- count +1
      }
      if(first) {
        nodes[[i]] <- edges
        first <- FALSE
      } else {
        nodes[[i]] <- append(nodes[[i]],edges)
      }

      ### this does not do anything yet. needs moving up to where each
      ### new node name is created and add an edge
      if(runif(1) < lower ) {
        prevd <- floor(runif(1,1,i+1))
        prev <- floor(runif(1,1,length(nodes[[prevd]])+1))
        print("i would connect to")
        print(nodes[[prevd]][prev])
      }

      alledges[[nodes[[i-1]][j]]]<- edges
      
    }
    print(nodes[[i]])
    allnodes <- append(allnodes,nodes[[i]])

  }
  allnodes <- append(allnodes,"t")
  for(j in 1:length(nodes[[depth]])) {
    alledges[[nodes[[depth]][j]]] <- "t"
  }

  
  print(alledges)
  g <- new("graphNEL",nodes=allnodes,edgeL=alledges,"directed")

  ## should not need to do this but some bug seems to be evident if we do not
  g <- addEdge("t","1",g,1)
  g <- removeEdge("t","1",g)

  return(g)
}

### THIS MAY HAVE BUGS!
### Not a full implementation! This uses the Fleischer maximum commodity
### flow algorithm to compute an approximate max-flow. THIS IS NOT THE
### BEST WAY TO DO THIS! but it is all I had "off the shelf"
### Demands must be less than
### a valid flow. Once the max-flow is computed it can determine the best
### s-t cut (the dual problem).
### g - graphNEL with capacities set
### s the starting node (character)
### t the finishing node (characeter)
### progress - printed progress bar if TRUE
### e - default=0.05, bound on approximate max-flow
### tolerance = 0.05 anything with less than this portion of the capacity
###                  in the flow will be assumed to be a congesting flow
###                  and will be used to form the cut
rg.st.cut <- function(g,s,t,progress=TRUE,e=0.05,tolerance=0.05,res=NULL) {

  ### need to do this better as this will be very slow
  validflow <- min(rg.edgeVector(g,"capacity"))
  demand <- rg.set.demand(validflow,c(s,t))
  if(is.null(res)) {
    res <- rg.fleischer.max.concurrent.flow(g,progress=progress,demands=demand,
                                            e=e)
  }
  flow <- as.double(edgeData(res$gflow,attr="weight"))
  cap <- as.double(edgeData(res$gflow,attr="capacity"))
  #cat("(cap-flow)/cap=",(cap-flow)/cap,"\n")
                   
  cuts <- NULL
  gcut <- g
  for(e in edgeNames(g)) {
    from <- strsplit(e,"~")[[1]][1]
    to <- strsplit(e,"~")[[1]][2]
    cap <- as.double(edgeData(g,from=from,to=to,attr="capacity"))
    flow <- as.double(edgeData(res$gflow,from=from,to=to,attr="weight"))
    if( (cap-flow)/cap <= tolerance ) {
      cuts <- c(cuts,from,to)
      gcut <- removeEdge(from,to,gcut)
    }
  }
  numcomp <- length(connComp(gcut))
  if(numcomp !=2 ) {
    cat("Error, Number of components=",numcomp," it should be 2\n")
    cat("Returning res, consider passing this back in with lower tolerance\n")
    return(res)
  } else {
    if("1" %in% connComp(gcut)[[1]] ) {
      if ("t" %in% connComp(gcut)[[2]] ) {

      } else {
        cat("Something is wrong expecting 1 and t to be in ",
            "in different components but they are not\n")
      }
    } else if ("t" %in% connComp(gcut)[[1]] ) {
      if ("1" %in% connComp(gcut)[[2]] ) {
        
      } else {
        cat("Something is wrong expecting 1 and t to be in ",
            "in different components but they are not\n")
      }

    } else {
      cat("Something is wrong expecting 1 and t to be in ",
          "in different components but they are not\n")
    }
  }
  
  return(cuts)
}


rg.analyse.tests <- function(res) {

  var <- data.frame()
  for(N in unique(res$N)) {
    var <- rbind(var,lapply(res,subset,res$N==N & res$L==res$N * res$N))
  }
  return(var)
}

rg.test.diameter <- function() {
    Size <- seq(100,800,100)
    Diameter <- c()
    for(N in Size) {
        g <- rg.generate.random.graph(N,2,25)
        Diameter <- c(Diameter,diameter(igraph.from.graphNEL(g)))
    }
    res <- data.frame(Size=Size,Diameter=Diameter)
    return(res)
}

rg.test.integer.plot <- function(res) {
  ## convert to long format
  resl <- melt(res,measure.vars=c("MF","FF"),variable.name="Algorithm")
  ## summarise
  ress <- summarySE(resl,"value",groupvars=c("Algorithm","Load"))
  ## plot!
  pd <- position_dodge(0.5) # offest so errorbars do not obscure
  ggplot(ress, aes(x=Load,y=value, colour=Algorithm)) + geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=1, position=pd) + geom_line(position=pd) + geom_point(position=pd, size=3, shape=21, fill="white") + ylab("Blocking Probability") + expand_limits(y=0)

  ggplot(res, aes(x=factor(Load))) + geom_boxplot(aes(y=FF),position=pd,width=0.25) + geom_boxplot(aes(y=MF),position=pd,width=0.25)
  
}

### NOT COMPLETE YET, the augmentation works but  it crashes in 
rg.test.integer <- function(emin,emax,estep,repeats,N=10,wavelengths=4,e=0.1) {
  g <- rg.generate.random.graph(N,2,25)
  M <- length(edgeMatrix(g))/2
  g <- rg.set.capacity(g,1.0)

  au <- rg.augment.graph.for.wavelengths(g,wavelengths)
  ga <- au$ga
  linkgroupmap <- au$linkgroupmap
  
  ga <- rg.set.weight(ga,0.0)
  sp.blocking <- c()
  flow.blocking <- c()
  flowint.blocking <- c()
  ff.blocking <- c()
  ff.gamma <- c()
  int.gamma <- c()
  int.mgamma <- c()
  ff.mgamma <- c()
  dcount <- c()
  int.results <- list()
  for(i in seq(emin,emax,estep) ) {
    for(j in 1:repeats) {
      demands <- rg.gen.demands(g,i)
      ## adapt the demands to the new source/sink nodes.
      for(j in names(demands) ){
        demands[[j]]$source <- paste0(demands[[j]]$source,"s")
        demands[[j]]$sink <- paste0(demands[[j]]$sink,"t")
      }
      numdems <- length(demands)
      dcount <- c(dcount, numdems)
      
      
      ff.results <- rg.sp.ff.max.concurrent.flow(ga,demands)
      accepted <- rg.count.accepted.demands(
        rg.check.blocked.demands(ff.results$demands,ga))
      ff.blocking <- c(ff.blocking,(numdems-accepted)/numdems)
      ff.gamma <- c(ff.gamma,ff.results$gamma)
      weights <- as.double(edgeData(ff.results$gflow,attr="weight"))
      capacities <- as.double(edgeData(ff.results$gflow,attr="capacity"))
      gamma <- (capacities - weights) / capacities
      mgamma <- mean(gamma)
      ff.mgamma <- c(ff.mgamma,mgamma)
      
      int.results <- rg.max.concurrent.flow.int(ga,demands,e=e)
      ch.demands <- rg.check.blocked.demands(int.results$demands,ga)
      accepted <- rg.count.accepted.demands(ch.demands)
      gtmp <- rg.max.concurrent.flow.graph(ga,ch.demands)
      weights <- as.double(edgeData(gtmp,attr="weight"))
      capacities <- as.double(edgeData(gtmp,attr="capacity"))
      gamma <- (capacities - weights) / capacities
      mgamma <- mean(gamma)
      int.mgamma <- c(int.mgamma, mgamma)
      flowint.blocking <- c(flowint.blocking,(numdems-accepted)/numdems)
      int.gamma <- c(int.gamma, min(gamma))
      
      results <- data.frame("Load"=dcount,
                            "FFB" = ff.blocking,
                            "IntB"=flowint.blocking)
      
      print(results[nrow(results),],digits=4)
    }
  }
  results <- data.frame("Load"=dcount,
                        "FF" = ff.blocking,
                        "MF"=flowint.blocking)
  return(results)
}




# Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean,
##   and confidence interval (default 95%).  data: a data frame.
##   measurevar: the name of a column that contains the variable to be
##   summariezed groupvars: a vector containing names of columns that
##   contain grouping variables na.rm: a boolean that indicates
##   whether to ignore NA's conf.interval: the percent range of the
##   confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm)
                              )
                          },
                    measurevar,
                    na.rm
             )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

rg.check.blocked.demands <- function # Checks demands in order to lower flow if blocked
### The results from a max-flow (min-congestion) algorithm may assign
### more traffic than is actually possible. This function will test
### each demand (and path in demand) in turn and only allocate the flow that is possible
### thus it is possible to see which demands are blocked by comparing flows against demands
### Note that trying the demands in order is not the optimum. It should be randomised
### rather trying it in order. Another possibility is ordering the demands/paths in
### order of length, shortest first. This way more demands are likely to be accepted.
(demands, ##<< demand list in standard format with paths
  g       ##<< graphNEL object with capacity attribute on edges
 ) {
  
  g <- rg.set.all.graph.edge.weights(g)
  for(nd in names(demands)) {
    d <- demands[[nd]]
    ##reassign the flow to the actual allocated
    flow <- 0
    for(p in names(d$paths)) {
      #cat(nd,":",p,"=",d$paths[[p]],"\n")
      pv <- as.vector(strsplit(p,"|",fixed=TRUE)[[1]])
      fromlist <- pv[1:length(pv)-1]
      tolist <- pv[2:length(pv)]
      
      weights <- as.double(edgeData(g,from=fromlist,
                                    to=tolist,att="weight"))
      capacities <- as.double(edgeData(g,from=fromlist,
                                       to=tolist,att="capacity"))

      minfree <- min(capacities-weights)

      free <- (capacities-weights)
      ## only for debugging
      ##minind <- which(free == minfree)
      ##for(m in minind) {
      ##  cat("Most constrained ",fromlist[m],"->",tolist[m],
      ##      weights[m],capacities[m]," ")
      ##}
      ## end only for debugging
      
      if (minfree < d$paths[[p]]) {
        ##cat("blocked")
        demands[[nd]]$paths[[p]] <- minfree
        flow <- flow + minfree
      } else {
        flow <- flow + d$paths[[p]]
      }
      ##cat("\n")
      weights <- weights + demands[[nd]]$paths[[p]]
      edgeData(g,from=fromlist,
               to=tolist,
               att="weight") <- weights
      
    }
    demands[[nd]]$flow <- flow
    
  }
  demands
  ### Demands with only the accepted traffic flow up to the limit of the graph capacities
  ### other demands will have their flow set to the maximum allowable.
}


rg.path.lengths <- function(demands) {
  lengths <- c()
  for(i in demands) {
    for(j in names(i$paths)) {
      lengths <- c(lengths,length(strsplit(j,"\\|")[[1]])-1)
    }
  }
  return(lengths)
}

## This was the routine used for the icc2012 paper
rg.run.tests <- function(min=8,max=20,averages=10,e=0.05,target=0.1,dirname=NULL) {
  vecN <- c()
  vecL <- c()
  vecSPgamma <- c()
  vecFLOWgamma <- c()
  vecINTgamma <- c()
  vecFLOWtime <- c()
  vecMAXmpls <- c()
  vecMAXbf <- c()
  run <- 1
  
  for(N in min:max ) {
    for(a in 1:averages) {
      cat("\nN=",N, "\n")

      scenario <- rg.test.int.versus.nonint.flow(N=N,L=NULL,target=target,e=e)
      if( ! is.null(dirname) ) {
        save(scenario,file=paste(dirname,"/",run,".robj",sep=""))
        cat("writing",paste(dirname,"/",run,".robj",sep=""),"\n")
      }
      vecN <- c(vecN,N)
      vecSPgamma <- c(vecSPgamma,scenario$sp.results$gamma)
      vecFLOWgamma <- c(vecFLOWgamma,scenario$flow.results$gamma)
      vecINTgamma <- c(vecINTgamma,scenario$int.results$bestgamma)
      vecFLOWtime <- c(vecFLOWtime,scenario$flow.results$runtime)
      vecMAXmpls <- c(vecMAXmpls,scenario$flow.results$max.mpls)
      vecMAXbf <- c(vecMAXbf,scenario$flow.results$max.bf)
      run <- run +1
    }
    
  }
  cat("\n")
  results <- data.frame(N=vecN,
                        SPgamma=vecSPgamma,
                        INTgamma=vecINTgamma,
                        FLOWgamma=vecFLOWgamma,
                        FLOWtime=vecFLOWtime,
                        MAXmpls=vecMAXmpls,
                        MAXbf=vecMAXbf
                        )
  
  return(results)
}

rg.run.tests.target <- function(min=0,max=1,step=0.1,averages=10,e=0.05,N=20,dirname=NULL) {
  vecTarget <- c()
  vecL <- c()
  vecSPgamma <- c()
  vecFLOWgamma <- c()
  vecINTgamma <- c()
  vecFLOWtime <- c()
  vecMAXmpls <- c()
  vecMAXbf <- c()
  run <- 1
  
  for(target in seq(min,max,step)) {
    for(a in 1:averages) {
      cat("\ntarget=",target, "\n")

      scenario <- rg.test.int.versus.nonint.flow(N=N,L=NULL,target=target,e=e)
      if( ! is.null(dirname) ) {
        save(scenario,file=paste(dirname,"/",run,".robj",sep=""))
        cat("writing",paste(dirname,"/",run,".robj",sep=""),"\n")
      }
      vecTarget <- c(vecTarget,target)
      vecSPgamma <- c(vecSPgamma,scenario$sp.results$gamma)
      vecFLOWgamma <- c(vecFLOWgamma,scenario$flow.results$gamma)
      vecINTgamma <- c(vecINTgamma,scenario$int.results$bestgamma)
      vecFLOWtime <- c(vecFLOWtime,scenario$flow.results$runtime)
      vecMAXmpls <- c(vecMAXmpls,scenario$flow.results$max.mpls)
      vecMAXbf <- c(vecMAXbf,scenario$flow.results$max.bf)
      run <- run +1
    }
    
  }
  cat("\n")
  results <- data.frame(target=vecTarget,
                        SPgamma=vecSPgamma,
                        INTgamma=vecINTgamma,
                        FLOWgamma=vecFLOWgamma,
                        FLOWtime=vecFLOWtime,
                        MAXmpls=vecMAXmpls,
                        MAXbf=vecMAXbf
                        )
  
  return(results)
}


rg.run.tests.process <- function(dirname) {
  vecN <- c()
  vecL <- c()
  vecSPgamma <- c()
  vecSPmaxflow <- c()
  vecFLOWgamma <- c()
  vecFLOWmaxflow <- c()
  vecINTgamma <- c()
  vecINTmaxflow <- c()
  vecFLOWtime <- c()
  vecMAXmpls <- c()
  vecMAXbf <- c()
  vecRatio <- c()
  run <- 1
  file <- paste(dirname,"/",run,".robj",sep="")
  while(file.exists(file)) {
    load(file)
    N <- length(nodes(scenario$g))
    vecN <- c(vecN,N)
    vecSPgamma <- c(vecSPgamma,scenario$sp.results$gamma)
    vecSPmaxflow <- c(vecSPmaxflow,
                      rg.calc.max.flow(scenario$sp.results$demands,scenario$sp.results$gamma))
    vecFLOWgamma <- c(vecFLOWgamma,scenario$flow.results$gamma)
    vecFLOWmaxflow <- c(vecFLOWmaxflow,
                        rg.calc.max.flow(scenario$flow.results$demands,scenario$flow.results$gamma))
#    vecFLOWmaxflow <- c(vecFLOWmaxflow,
#                        rg.calc.max.flow(scenario$flow.results$demands,0))
    vecINTgamma <- c(vecINTgamma,scenario$int.results$bestgamma)
    vecINTmaxflow <- c(vecINTmaxflow,
                       rg.calc.max.flow(scenario$int.results$intdemands,
                                        scenario$int.results$bestgamma))
#    vecINTmaxflow <- c(vecINTmaxflow,rg.calc.max.flow(scenario$int.results$intdemands,0))
    vecFLOWtime <- c(vecFLOWtime,scenario$flow.results$runtime)
    vecMAXmpls <- c(vecMAXmpls,scenario$flow.results$max.mpls)
    vecMAXbf <- c(vecMAXbf,scenario$flow.results$max.bf)
    vecRatio <- c(vecRatio,scenario$flow.results$ratio)
    run <- run+1
    file <- paste(dirname,"/",run,".robj",sep="")
  }
  print(length(vecSPmaxflow))
  print(length(vecFLOWmaxflow))
  print(length(vecINTmaxflow))
  results <- data.frame(N=vecN,
                        SPgamma=vecSPgamma,
                        SPmaxflow=vecSPmaxflow,
                        INTgamma=vecINTgamma,
                        INTmaxflow=vecINTmaxflow,
                        FLOWgamma=vecFLOWgamma,
                        FLOWmaxflow=vecFLOWmaxflow,
                        FLOWtime=vecFLOWtime,
                        MAXmpls=vecMAXmpls,
                        MAXbf=vecMAXbf,
                        Ratio=vecRatio
                        )
}

rg.calc.max.flow <- function(demands,gamma) {
  flow <- 0.0
  for(i in demands) {
    flow <- flow + i$flow / (1 - gamma)
  }
  return(flow)
}

rg.run.tests.process.maxflow <- function(dirname) {
  vecN <- c()
  vecL <- c()
  vecSPgamma <- c()
  vecFLOWgamma <- c()
  vecINTgamma <- c()
  vecFLOWtime <- c()
  vecMAXmpls <- c()
  vecMAXbf <- c()
  vecRatio <- c()
  run <- 1
  file <- paste(dirname,"/",run,".robj",sep="")
  while(file.exists(file)) {
    load(file)
    N <- length(nodes(scenario$g))
    vecN <- c(vecN,N)
    vecSPgamma <- c(vecSPgamma,scenario$sp.results$gamma)
    vecFLOWgamma <- c(vecFLOWgamma,scenario$flow.results$gamma)
    vecINTgamma <- c(vecINTgamma,scenario$int.results$bestgamma)
    vecFLOWtime <- c(vecFLOWtime,scenario$flow.results$runtime)
    vecMAXmpls <- c(vecMAXmpls,scenario$flow.results$max.mpls)
    vecMAXbf <- c(vecMAXbf,scenario$flow.results$max.bf)
    vecRatio <- c(vecRatio,scenario$flow.results$ratio)
    run <- run+1
    file <- paste(dirname,"/",run,".robj",sep="")
  }
  results <- data.frame(N=vecN,
                        SPgamma=vecSPgamma,
                        INTgamma=vecINTgamma,
                        FLOWgamma=vecFLOWgamma,
                        FLOWtime=vecFLOWtime,
                        MAXmpls=vecMAXmpls,
                        MAXbf=vecMAXbf,
                        Ratio=vecRatio
                        )
}

rg.demands.to.csv <- function(demands,file,paths=TRUE) {
  cat("",file=file)
  if(paths) {
    for(i in demands) {
      for(j in names(i$paths)) {
        cat(i$paths[[j]],",",gsub("\\|",",",j),"\n",sep="",file=file,append=TRUE)
      }
    }
  } else {
    for(i in demands) {
        cat(i$demand,",",i$source,",",i$sink,"\n",sep="",file=file,append=TRUE)
    }
  }
}

rg.smart.round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}


rg.karcius.run.all <- function(gnames=NULL,targets=NULL,zoo=zoo,dir="./",incomm=NULL,progress=FALSE,e=0.1) {
  g <- NULL
  for(gname in gnames) {
    if(gname == "lattice") {
      g <- NULL
      gi <- make_lattice(c(2,3))
      g <- igraph.to.graphNEL(as.directed(gi))
    } else {
      gi <- zoo[[gname]]
      g <- igraph.to.graphNEL(as.directed(simplify(gi)))
    }
    N <- length(nodes(g))
    if(is.null(incomm)) {
      comm <- rg.gen.demands(g,val=10*rlnorm(N*N,16,1)/exp(16.5))
    } else {
      comm <- incomm
    }
    foreach(target=targets) %dopar% {
      #for(target in targets){
      network.name <- paste0(dir,gname,"-",as.character(target))
      cat(network.name,"\n")
      results <- rg.karcius(g=g,comm=comm,network.name=network.name,target=target,e=e,progress=progress)
    }
  }
}

rg.karcius <- function(g=NULL,comm=NULL,dir="./",network.name="unknown-net",target=0.3,e=0.1,progress=FALSE) {
  N <- length(nodes(g))
  g <- rg.set.capacity(g,10000)
  #comm <- rg.gen.demands(g,val=10)
  g <- rg.set.weight(g,0.0)
  results <- rg.minimum.congestion.flow(g,comm,e=e,progress=progress)
  ## now work out how much to scale demand to meet a max-flow target
  multiplier <- (1-target)/(1-results$gamma)
  ## and now scale the results
  comm <- rg.rescale.demands(comm,multiplier,integer=TRUE)
  flow.results <- rg.minimum.congestion.flow(g,comm,e=e,progress=progress)
  int.results <- rg.max.concurrent.flow.int(g,flow.results$demands,e=e,progress=progress)
  sp.results <- rg.sp.max.concurrent.flow(g,comm)
  for(i in 1:length(flow.results$demands)) {
    rounded.flows <- rg.smart.round(as.numeric(flow.results$demands[[i]]$paths))
    if(sum(rounded.flows) > flow.results$demands[[i]]$demand ) {
      rounded.flows[which.max(rounded.flows)] <- rounded.flows[which.max(rounded.flows)] -
        sum(rounded.flows) + flow.results$demands[[i]]$demand
    } else if(sum(rounded.flows) < flow.results$demands[[i]]$demand ) {
      rounded.flows[which.max(rounded.flows)] <- rounded.flows[which.max(rounded.flows)] +
        sum(rounded.flows) - flow.results$demands[[i]]$demand
    }
    if(sum(rounded.flows) > flow.results$demands[[i]]$demand ) {
      cat(i,"too high","\n")
    } else if(sum(rounded.flows) > flow.results$demands[[i]]$demand ) {
      cat(i,"too low","\n")
    }
    count <- 1
    for(j in names(flow.results$demands[[i]]$paths)) {
      flow.results$demands[[i]]$paths[[j]] <- rounded.flows[count]
      count <- count+1
    }
  }
  flow.results$gflow <- rg.max.concurrent.flow.graph(flow.results$gflow,flow.results$demands)
  flow.results$gamma <- rg.mcf.find.gamma(flow.results$gflow)
  output <- list()
  output$g <- g
  output$sp.results <- sp.results
  output$flow.results <- flow.results
  output$int.results <- int.results
  output$comm <- comm
  network.name <- paste0(dir,network.name,"-")
  cat(network.name,"sp=",sp.results$gamma,",flow=",flow.results$gamma,"int=",int.results$gamma,"\n")
  cat("sp=",sp.results$gamma,",flow=",flow.results$gamma,"int=",int.results$gamma,"\n",
      file=paste0(network.name,"results.txt"))
  file <- paste0(network.name,"integer-flows.csv")
  rg.demands.to.csv(int.results$demands,file)
  file <- paste0(network.name,"sp-flows.csv")
  rg.demands.to.csv(sp.results$demands,file)
  file <- paste0(network.name,"split-flows.csv")
  rg.demands.to.csv(flow.results$demands,file)
  file <- paste0(network.name,"graph.csv")
  write.table(as_adjacency_matrix(igraph.from.graphNEL(g),sparse=FALSE,attr="capacity"),file=file)
  file <- paste0(network.name,"demands.csv")
  rg.demands.to.csv(comm,file=file,paths=FALSE)
  return(output)
}


## test used for PURSUIT TM theory and ICC paper
## N - the number of nodes in the graph
## L - number of demands to generate (NULL will do N(N-1) demands)
## target - the demands is automatically scaled to give at least this much
##          free capacity on the most congested edge for the fractional
##          multicommodity flow (try approx 0.3 for first attempt)
## e - optimisation parameter
## g - a graph to try with, if NULL generate a graph
rg.test.int.versus.nonint.flow <- function(N=8,L=NULL,target=0.0,e=0.1,g=NULL) {
  if(is.null(g)) {
    g <- rg.generate.random.graph(N,2,25)
    
  } else {
    N <- length(nodes(g))
  }
  # number of edges
  M <- length(edgeMatrix(g))/2
  if(is.null(L)) {
    ## generate using lognormal traffic as in
    ## Nucci, Antonio and Sridharan, Ashwin and Taft, Nina,
    ## The problem of synthetically generating IP traffic matrices: initial recommendations,
    ## SIGCOMM Comput. Commun. Rev., vol 35(3) July 2005 pp19-32
    ## this generates a mean of 1 with a variance that matches the paper
    comm <- rg.gen.demands(g,val=rlnorm(N*N,16,1)/exp(16.5))
  } else {
    comm <- rg.gen.demands(g,L,val=rlnorm(L,16,1)/exp(16.5))
  }

  ## set some random capacities (this should be an input parameter)
  g <- rg.set.capacity(g,runif(M,5,10))
  g <- rg.set.weight(g,0.0)
  ##results <- rg.sp.max.concurrent.flow(g,comm)
  results <- rg.minimum.congestion.flow(g,comm,e=e,progress=TRUE)
  ## now work out how much to scale demand to meet a max-flow target
  multiplier <- (1-target)/(1-results$gamma)
  ## and now scale the results
  comm <- rg.rescale.demands(comm,multiplier)
  ## calculate the flow when using strict shortest path routing
  sp.results <- rg.sp.max.concurrent.flow(g,comm)
  ## calculate the flow when using minimum congestion flow (fractional)
  ## but this time with the scaled demand
  flow.results <- rg.minimum.congestion.flow(g,comm,e=e,progress=TRUE)
  ## count how many flows we have, 
  gcount <- rg.count.edge.flows(g,flow.results$demands)
  flow.results$max.mpls <- max(as.double(edgeData(gcount,att="weight")))
  flow.results$max.bf <- max(as.double(lapply(graph::edges(g),length)))
  ## calculate the integer results
  int.results <- rg.max.concurrent.flow.int.c(g,flow.results$demands,
                                              e=e,progress=TRUE)

  output <- list()
  output$g <- g
  output$sp.results <- sp.results
  output$flow.results <- flow.results
  output$int.results <- int.results
  output$comm <- comm
  cat("sp=",sp.results$gamma,",flow=",flow.results$gamma,"int=",int.results$bestgamma,"\n")
  return(output)
}

rg.gen.all.pairs.demands <- function(g,val=1.0) {
  nodes <- nodes(g)
  comm <- list()
  n=1;

  k <- 1
  vallen <- length(val)
  for(i in nodes) {
    for(j in nodes) {
      if ( i != j ) {
        comm[[as.character(n)]] <- list(source=i,sink=j,demand=val[k])
        k <- k + 1
        if(k > vallen) k <- 1
        n <- n+1
      }
    }
  }
  return(comm)
    
}


## Generate a random graph and augmented graph representing the
## wavelength switched paths through the network.
## n - number of nodes
## mindeg - minimum node degree (default = 3)
## maxdeg - maximum node degree (default = 7)
## grid - wavelength numbers to be used on each span (ITU grid 1-73), special
##        value of zero means the augmented graph is the same as the original
##        graph
## lamdacap - the bandwidth of the wavelengths in Gb/s (default 1.0)
##
## output - graph$g the original graph
##          graph$lg the augmented lambda graph
##          graph$accessnodes the list of access nodes (all of the original graph)

rg.generate.random.lambda.graph <- function(g,n=10,mindeg=3,maxdeg=7,grid=c(0),lambdacap=1.0) {

  #g <- rg.generate.random.graph(n,mindeg,maxdeg)

  nodelist <- c()

  accessnodes <- c()
  for (n in nodes(g)) {
    accessnodes <- c(accessnodes,paste("A",n,sep="."))
  }
  nodelist <- accessnodes
  
  for (n in nodes(g)) {
    for(l in grid) {
      nodelist <- c(nodelist,paste("A",n,l,sep="."))
    }
  }

  
  for (e in rg.edgeL(g)) {
    for(l in grid) {
      s <- e[1]
      t <- e[2]
      nodelist <- c(nodelist,paste(s,t,l,sep="."))
    }
  }
  print(nodelist)

  ag <- new("graphNEL",nodes=nodelist,edgemode="directed")


  for(l in grid) {
    for(n in nodes(g)) {
      ag <- addEdge(paste("A",n,sep="."),
                    paste("A",n,l,sep="."),
                    ag,
                    c(0))
      for(u in graph::edges(g)[[n]] ) {
        if( u == n) next
        
        ag <- addEdge(paste("A",n,l,sep="."),
                      paste(n,u,l,sep="."),
                      ag,
                      c(0))
      }
    }
  }
  for(l in grid) {
    for (e in rg.edgeL(g)) {
      s <- e[1]
      t <- e[2]

      ag <- addEdge(paste(s,t,l,sep="."),
                    paste("A",t,sep="."),
                    ag,
                    c(0))
      
      for(u in graph::edges(g)[[t]] ) {
        if( u == t) next
        ## s.t to t.access
        
        ## s.t to t.u

        ag <- addEdge(paste(s,t,l,sep="."),
                      paste(t,u,l,sep="."),
                      ag,
                      c(0))
      }
    }
  }

  ag <- rg.set.capacity(ag,lambdacap)
  graph <- list()
  graph$g <- g
  graph$ag <- ag
  graph$anodes <- accessnodes
  return(graph)
}


rg.sp.lambda.max.concurrent.flow <- function(g,demands) {
  for(c in names(demands)) {
    demands[[c]]$flow <- 0
  }

  nodes <- nodes(g)
  gsol <- g
  punults.done <- rep(FALSE,length(nodes))

  rg.set.all.graph.edge.weights(gsol)
  
  results <- list()
  results$g <- g
  results$demands <- demands

  ccount <- 1
  for (i in demands) {
    startind=which(nodes==i$source)
    finind=which(nodes==i$sink)
    if(!punults.done[startind]) {
      g.sp[[startind]] <- dijkstra.sp(g,start=i$source)$penult
      punults.done[startind]=TRUE
    }
    path <- extractPath(startind,finind,g.sp[[startind]])
    print(path)
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
  
  results$demands <- demands
  results$gflow <- gsol
  results$lambda <- lambda
  results$gamma <- gamma
    
  return(results)
}

rg.dfs <- function(g) {
  data <- list()
  data$color <- rep(0,length(nodes(g))+1)
  data$pred <- rep(NULL,length(nodes(g)))
  for (u in nodes(g)) {
    cat("I am visiting", u,"\n")
    if(data$color[as.integer(u)] == 0 ) {
      data <- mydfs.visit(g,u,data)
    }
  }
}

rg.dfs.visit <- function(g,u,data) {
  data$color[as.integer(u)] <- 1
  print(data$color[as.integer(u)])
  for (v in graph::edges(g,u)[[1]]) {
    cat("testing edge:",v,":\n")
    if(data$color[as.integer(v)] == 0 ) {
      data$pred[as.integer(v)] <- as.integer(u)
      data <- mydfs.visit(g,v,data)
    }
  }
  
  data$color[as.integer(u)] <- 2
  print(data$color[as.integer(u)])
  cat("I have finished at ",u,"\n")
  return(data)
}


analyse.runs <- function(nums=NULL,maxattempts=1,progress=FALSE,e=0.005) {

  if(is.null(nums)) {
    files <- Sys.glob("run[0-9]*")
  } else {
    files <- c()
    for(num in nums) {
      files <- c(files,Sys.glob(paste("run",num,".*",sep="")))
    }
  }

  files <- files[order(as.numeric(gsub("^.*run([0-9]*).*$","\\1",files,files)))]
  foundratio <- as.double(c())
  countzero <- 0
  neverfound <- 0
  faillist <- c()
  numpaths <- 0
  countruns <- 0
  for(file in files) {
##  for(i in seq(1:50)) {
##    file <- paste("run",i,".robj",sep="")
    cat("----------------- Doing ",file,"\n")
    load(file)
    if(is.null(res$error)) {
      scenario <- res$scenario
      if(scenario$intgammalist[[1]] > 0.0) {

        attempts <- 0
        countgamma <- 0
        while(countgamma == 0 && attempts < maxattempts) {
          retlist <- analyse.inner(scenario,foundratio,file,progress=progress,e=e)
          countgamma <- retlist$res$countgamma
          attempts <- attempts + 1
          if(retlist$res$countgamma == 0)
            countzero <- countzero +1
        }
        if(retlist$res$countgamma == 0) {
          neverfound <- neverfound +1
          faillist <- c(file,faillist)
        }
        countruns <- countruns + 1
        numpaths <- numpaths + retlist$numpaths

        foundratio <- retlist$foundratio
      }
    } else {
      ##cat("was null\n")
    }
  }
  cat("numzero=",countzero,
      "min=",min(foundratio)," mean=",mean(foundratio)," max=",max(foundratio),
      "\n")
  cat("av numpaths=",numpaths/countruns,"\n")
  if(neverfound >0) {
    cat("Did not find",neverfound," for:\n")
    cat(gsub("[^0-9]+([0-9]+)[^0-9]+","\\1",faillist),"\n")
  }
  
  

}

analyse.inner.part1 <- function(nums,e=0.005) {

  if(is.null(nums)) {
    files <- Sys.glob("run[0-9]*")
  } else {
    files <- c()
    for(num in nums) {
      files <- c(files,Sys.glob(paste("run",num,".*",sep="")))
    }
  }
  
  files <- files[order(as.numeric(gsub("^.*run([0-9]*).*$","\\1",files,files)))]
  for(file in files) {
##  for(i in seq(1:50)) {
##    file <- paste("run",i,".robj",sep="")
    cat("----------------- Doing ",file,"\n")
    load(file)
    if(is.null(res$error)) {
      scenario <- res$scenario
      if(scenario$intgammalist[[1]] > 0.0) {
        results <- rg.minimum.congestion.flow(scenario$g,scenario$commr,e=e,progress=TRUE,permutation="lowest")

      }
    }
  }
  results$scenario <- scenario
  return(results)
}

analyse.inner <- function(scenario,foundratio,file,progress,e) {
  
  results <- rg.minimum.congestion.flow(scenario$g,scenario$commr,e=e,progress=progress,permutation="lowest")
  #res <- rg.try.single.demands(scenario$g,scenario$commr,e=0.005,permutation="lowest",
  #                             scenario=scenario,progress=progress)
  #numpaths <- rg.count.paths(res$demands)

  numpaths <- rg.count.paths(results$demands)

  res <- rg.max.concurrent.flow.int.c(
                                      scenario$g,
                                      results$demands,
                                      e=e,
                                      #eInternal=0.001,
                                      scenario=scenario,
                                      progress=progress,
                                      permutation="lowest"
                                      )
  res$scenario <- scenario
  foundratio <- c(foundratio,as.double(res$countgamma)/as.double(res$phases))
  cat(file,":",res$countgamma,"/",res$phases,
      " gamma=",res$scenario$intgammalist[[1]],
      "best found=",res$bestgamma,"\n")
  cat("pathdifcount=",res$pathdiffcount,"\n")
  cat("phasepathdiffcount=",res$phasepathdiffcount,"\n")
  retlist <- list()
  retlist$res <- res
  retlist$foundratio <- foundratio
  retlist$numpaths <- numpaths
  return(retlist)
}

rg.try.single.demands <- function(g,demands,e=0.1,progress=FALSE,permutation="fixed",scenario) {

  newdemands <- list()
  for(i in 1:length(demands)) {
    tmp <- list()
    tmp[["1"]] <- demands[[i]]
    results <- rg.minimum.congestion.flow(g,tmp,e=0.01,progress=progress,permutation=permutation)
    newdemands[[i]] <- results$demands[[1]]
  }
  res <- rg.max.concurrent.flow.int.c(g,newdemands,e=e,progress=progress,scenario=scenario,permutation=permutation)
  return(res)
}





rg.test.idea <- function(maxtime=1.0,N=9,L=10,filebase=NULL) {
  ## maximum time in hours for a run
  ## this is based upon 24e6 patch combinations taking 106 seconds
  hoursperpath <- 1470 / ( 3600 * 253e6  )
  maxpathcomb <- maxtime / hoursperpath

  filename <- paste(filebase,".robj",sep="")
  scenario <- rg.gen.random.scenario(N=N,L=L)

  starttime <- Sys.time()
  cat("starting at ",format(starttime),"\n")

  
  
  if(scenario$pathcomb > maxpathcomb) {
    cat("scenario =",scenario$pathcomb," too big\n")
#        as.double(maxpathcomb)/as.double(scenario$pathcomb) * maxtime,
    
    cat("it would have taken ",
        scenario$pathcomb * hoursperpath,
        "hours\n")
    res <- list()
    res$scenario <- scenario
    res$error <- "Scenario too big"
    if(!is.null(filebase)) {
      save(res,file=filename)
    }
    
    return(res)
  }
  estfinished <- starttime + scenario$pathcomb * hoursperpath * 3600
  cat("Estimate finish at ",format(estfinished),"\n")
  scenario <- rg.testeverypath.c(scenario)
  res <-
    rg.fleischer.max.concurrent.flow.restricted.c(
                                                  scenario$g,
                                                  scenario$results$demands,
                                                  e=0.01,
                                                  bestgamma=
                                                  scenario$intgammalist[[1]])
  res$scenario <- scenario
  cat("N=",N," L=",L," pathcomb=",scenario$pathcomb,
      "gammahits=",res$countgamma,"out of ",res$phases," phases\n")
  
  if(!is.null(filebase)) {
    save(res,file=filename)
  }
  
  return(res)
}


### Generate a random graphNEL and scaled commodities so that gamma is 0.9
### N    Number of nodes in graph
### L    Number of commodities
### returns output$pathcomb - number of combinations of paths
###         output$g  - graphNEL with capacity
###         output$results - from running rg.minimum.congestion.flow()
###         output$commr - the commodities (scaled to give gamma a reasonable value)
###         output$comm - the original commodities
rg.gen.random.scenario <- function(N=6,L=8,target=0.9,progress=FALSE,e=0.02) {
  g <- rg.generate.random.graph(N,3,6)
  M <- length(edgeMatrix(g))/2
  comm <- rg.gen.demands(g,L,runif(L,1,5))
  g <- rg.set.capacity(g,runif(M,5,10))
  g <- rg.set.weight(g,0.0)
  results <- rg.minimum.congestion.flow(g,comm,e=e,progress=progress)
  ## need to customize this
  multiplier <- (1-target)/(1-results$gamma)
  commr <- rg.rescale.demands(comm,multiplier)
  results <- rg.minimum.congestion.flow(g,commr,e=e,progress=progress)
  runtotal <- 1
  for(i in results$demands) {
    val <- length(names(i$paths))
    runtotal <- val * runtotal
  }
  output <- list()
  output$pathcomb <- runtotal
  output$g <- g
  output$results <- results
  output$comm <- comm
  output$commr <- commr
  return(output)
}

### Brute force search for the best path. It returns data with new values
### input data$g - graph graphNEL with capacity edge attribute set
###       data$results data$g data$comm data$commr data$pathcomb
###       this is as from rg.gen.random.scenario()
### returns data as above plus
###         data$intdemands
###         data$gintflow
###         WARNING: it will overboook
###         WARNING: it seems to keep the edges attribute from the graph
###                  it copied, this should be set to null.

rg.testeverypath.c <- function(data,progress=FALSE) {
  recordlen <- 100
  demands <- data$results$demands
  intdemands <-  demands

  g <- data$g
  vdemands <- as.double(lapply(demands,"[[","demand"))
  vcapacity <- as.double(edgeData(g,attr="capacity"))
  vdemandflow <- rep(0.0,length(vcapacity))
  L <- length(demands)
  
  edgeMap <- adjMatrix(g)
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

  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = data$pathcomb, style=3)
  } else {
    pb <- NULL
  }

  retlist <- .Call("rg_test_every_path_inner",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   progress,
                   pb,
                   recordlen,
                   environment()
                   );

  if( progress != FALSE) {
    close(pb)
    cat("Wrapping up and creating ",recordlen, "best graphs\n")
  }

  data$gintflowlist <- list()
  data$intdemandslist <- list()
  data$intgammalist <- retlist$gammas
  recorddemands <- retlist$pathptrs
  for(j in seq(1:recordlen) ) {
    g <- data$g
    g <- rg.set.weight(g,0)
    bestPathptr <- recorddemands[[j]] + 1
#    cat("test1\n")
#    print(bestPathptr[[i]])
#    print(demands[[i]])
#    cat("test2\n")
    for(i in seq(1:L)) {
      intdemands[[i]]$paths <- list()
      intdemands[[i]]$flow <- demands[[i]]$flow
      intdemands[[i]]$paths[[names(demands[[i]]$paths[bestPathptr[[i]]])]] <-
        demands[[i]]$flow
      path <- names(demands[[i]]$paths)[ bestPathptr[[i]] ]
      pv <- as.vector(strsplit(path,"|",fixed=TRUE)[[1]])
      
      val <- demands[[i]]$demand
      
      if(length(pv) == 1) {
        edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") <-
          edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") +val
      } else {
        
        fromlist <- pv[1:{length(pv)-1}]
        tolist <- pv[2:{length(pv)}]
        newvals <- as.double(edgeData(g,as.character(fromlist),as.character(tolist),"weight")) + val
        edgeData(g,as.character(fromlist),as.character(tolist),"weight") <- newvals
      }
      
    }

    data$gintflowlist[[j]] <- g

    data$intdemandslist[[j]] <- intdemands
  }
  data$gintflow <- data$gintflowlist[[1]]
  data$intdemands <- data$intdemandslist[[1]]
  return(data)


}


rg.testeverypath <- function(data,progress=FALSE) {

  ## record this many of the best, sorted by best gamma
  recordlen <- 100
  recordgamma <- rep(-Inf,recordlen)
  recorddemands <- rep(list(NULL),recordlen)

  L=length(data$results$demands)
  pathptr <- rep(1,times=L)
  pathsz <- c()
  demands <-  data$results$demands
  count <- 1
  for(i in demands) {
    pathsz[count] <- length(i$paths)
    count <- count +1
  }
  finished=FALSE
  count=1
  intdemands <-  demands
  for(d in names(intdemands)) {
    intdemands[[d]]$paths <- list()
    intdemands[[d]]$flow <- 0.0
  }
  intgraph <- data$g
  bestgamma <- -Inf
  bestPathptr <- c()
  count <- 0

  if(progress != FALSE) {
    pb <- txtProgressBar(min = 0,
                         max = data$pathcomb, style=3)
  } else {
    pb <- NULL
  }
  updatepb <- as.integer(ceiling( data$pathcomb / 100.0))

  while(!finished) {

    if(progress != FALSE && count %% updatepb == 0) {
      setTxtProgressBar(pb,count)
    }

    #print(pathptr)

    ## pathptr is counter to path selectors

    ## add load to each path
    ## for each demand (index)
    g <- data$g
    
    for(i in seq(1:L)) {
      path <- names(demands[[i]]$paths)[ pathptr[i] ]
      pv <- as.vector(strsplit(path,"|",fixed=TRUE)[[1]])

#     print(path)
      val <- demands[[i]]$demand
      
      if(length(pv) == 1) {
        edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") <-
          edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") +val
      } else {
        
        fromlist <- pv[1:{length(pv)-1}]
        tolist <- pv[2:{length(pv)}]
        newvals <- as.double(edgeData(g,as.character(fromlist),as.character(tolist),"weight")) + val
        edgeData(g,as.character(fromlist),as.character(tolist),"weight") <- newvals
      }
      
    }
    
    ## calc max load
    weights <- as.double(edgeData(g,attr="weight"))
    capacities <- as.double(edgeData(g,attr="capacity"))
    gamma <- min((capacities - weights)/capacities)
    ## record best gamma and pathptr so far
    if( gamma > bestgamma ) {
      bestgamma <- gamma
      bestPathptr <- pathptr
      cat("best gamma so far =",bestgamma,"\n")
    }
    if( gamma > recordgamma[recordlen] ) {
      recordgamma[recordlen] <- gamma
      recorddemands[[recordlen]] <- pathptr
      o <- order(recordgamma,decreasing=TRUE)
      recordgamma <- recordgamma[o]
      recorddemands <- recorddemands[o]
      
    }
    #cat("found=",gamma," best=",bestgamma,"\n")

    #cat(recordgamma,"\n")
    increment=TRUE
    for(i in seq(1:L)) {
      if(pathsz[i]>1) {
        if(pathptr[i]%%pathsz[i] == 0)
          incrementnext <- TRUE
        else
          incrementnext <- FALSE
        if(increment)
          pathptr[i] <- pathptr[i]%%pathsz[i]+1
        if(incrementnext && increment)
          increment <- TRUE
        else
          increment <- FALSE
      }
    }

    
    if(increment)
      finished <- TRUE
    count <-  count + 1
    # just for testing
    #if(count > 20)
     # finished=TRUE
  }

  if(progress != FALSE)
    close(pb)

  ## have found best, now fill integer flow graph and demands

  data$gintflowlist <- list()
  data$intdemandslist <- list()
  data$intgammalist <- recordgamma
#  cat("bestpathptr=",recorddemands[[1]],"\n")
  for(j in seq(1:recordlen) ) {
    g <- data$g

    bestPathptr <- recorddemands[[j]]
    for(i in seq(1:L)) {
      intdemands[[i]]$flow <- demands[[i]]$flow
      intdemands[[i]]$paths <- NULL
      intdemands[[i]]$paths[[names(demands[[i]]$paths[bestPathptr[i]])]] <-
        demands[[i]]$flow
      path <- names(demands[[i]]$paths)[ bestPathptr[i] ]
      pv <- as.vector(strsplit(path,"|",fixed=TRUE)[[1]])
      
      val <- demands[[i]]$demand
      
      if(length(pv) == 1) {
        edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") <-
          edgeData(g,as.character(pv[1]),as.character(pv[2]),"weight") +val
      } else {
        
        fromlist <- pv[1:{length(pv)-1}]
        tolist <- pv[2:{length(pv)}]
        newvals <- as.double(edgeData(g,as.character(fromlist),as.character(tolist),"weight")) + val
        edgeData(g,as.character(fromlist),as.character(tolist),"weight") <- newvals
      }
      
    }

    data$gintflowlist[[j]] <- g

    data$intdemandslist[[j]] <- intdemands
  }
  data$gintflow <- data$gintflowlist[[1]]
  data$intdemands <- data$intdemandslist[[1]]
  return(data)
}

calcDexplicit <- function(gdual) {
    return(sum(as.double(edgeData(gdual,attr="capacity"))*
        as.double(edgeData(gdual,attr="weight"))))
  }



###system.time(test1(gdual))
###   user  system elapsed 
### 44.084   0.283  44.589 
test1 <- function(g) {

  var <- 1.0000001
  val <- 1.0
  em <- edgeMatrix(g)


  
  for(i in seq(1:10000)) {
    w <- as.double(edgeData(g,attr="weight"))
    w <- w*var
    g <- rg.set.weight(g,w)
  }

  return(g)


}

###system.time(test2(gdual))
###   user  system elapsed 
### 43.785   0.214  44.427 
test2 <- function(g) {

  var <- 1.0000001
  val <- 1.0
  em <- edgeMatrix(g)


  
  for(i in seq(1:10000)) {
    w <- as.double(edgeData(g,attr="weight"))
    w <- w*var
    em <- edgeMatrix(g)
    if (match("weight",names(edgeDataDefaults(g)),nomatch=1)) {
      edgeDataDefaults(g,"weight") <- 1.0
    }
    edgeData(g,
             from=as.character(em[1,]),
             to=as.character(em[2,]),
             attr="weight") <- val
  }
  
  return(g)


}


rg.test2 <- function() {
  pb <- txtProgressBar(min=0,max=5,style=3)
  for(i in 1:5) {
    Sys.sleep(1)
    setTxtProgressBar(pb,i)
    Sys.sleep(1)
    close(pb)
  }
}

rg.explore.start <- function(g,demands,e=0.1,progress=FALSE,intdemands=NULL) {

  initresults <- rg.fleischer.max.concurrent.flow.c(g,demands,e=e)
  cat("initresult$lambda=",initresults$lambda,
      "initresults$beta=",initresults$beta,"\n")
  initdemands <- initresults$demands
#  res <- rg.explore(g,initdemands,e=e,progress=progress,intdemands=intdemands)

  res <- rg.explore.outer(g,initdemands,e=e,progress=progress,intdemands=intdemands)

  ## fix flow to equal demand
  for(d in 1:length(res)) {
    res[[d]]$paths[1] <- res[[d]]$demand
    res[[d]]$flow <- res[[d]]$demand
  }

  gint <- rg.max.concurrent.flow.graph(g,intdemands)
  gexplore <- rg.max.concurrent.flow.graph(g,res)

  gammaint <- rg.mcf.find.gamma(gint)
  gammaexplore <- rg.mcf.find.gamma(gexplore)


  for(d in 1:length(res)) {
    cat(names(res[[d]]$paths[1]), ", ",
        names(intdemands[[d]]$paths[1]))
    if( identical (names(res[[d]]$paths[1]),
                   names(intdemands[[d]]$paths[1]))) {
      cat("  OK \n")
    } else {
      cat(" !!!!!! Different \n")
    }      
        

  }
  
  cat("Int gamma=",gammaint," gexplore gamma=",gammaexplore,"\n")
  return(res)

}
rg.explore.outer <- function(g,demands,e=0.1,progress=FALSE,intdemands=NULL) {

  res <- rg.explore(g,demands,e=e,progress=progress,intdemands=intdemands)

  finished <- TRUE
  for(d in 1:length(res)) {
    if (length(res[[d]]$paths) > 1 )
      finished <- FALSE
    cat("num paths in ",d," = ",length(res[[d]]$paths),"\n")
  }

  if( ! finished ) {
    res <- rg.explore.outer(g,res,e=e,progress=progress,intdemands=intdemands)
  }
  
  return(res)
}


rg.explore <- function(g,initdemands,e=0.1,progress=FALSE,intdemands=NULL) {
  ## calc split flows

  builtdemands <- initdemands
  
  ## this line not strictly needed
  #res <- rg.fleischer.max.concurrent.flow.restricted.c(g,initdemands,
  #                                                   updateflow=FALSE,e=e)

  #cat("res$beta=",res$beta," res$betar=",res$betar,
  #    " res$lambda=",res$lambda,"\n")

  ## try one path at a time in each demand and calculate lambda, record maximum

  overallMaxMetric <- -Inf
  overallMaxMetricD <- 0

  demandpaths <- list()
  demandpathsbest <- list()
  j <- 1
  for(d in initdemands) {
    demandpaths[[j]] <- list()
    demandpathsbest[[j]] <- 0
    k <- 1
    for(p in names(d$paths)) {
      demandpaths[[j]][[k]] <- list()
      demandpaths[[j]][[k]]$name <- p
      demandpaths[[j]][[k]]$lambda <- 0.0
      demandpaths[[j]][[k]]$beta <- 0.0
      k <- k + 1
    }
    j <- j + 1
  }
  
  for(j in 1:length(initdemands)) {
    maxmetric <- -Inf
    maxmetrici <- i
    cat("demand number",j," num paths ",length(initdemands[[j]]$paths),"---------------\n")
    if (length(initdemands[[j]]$paths) > 1 ) {
      for(i in 1:length(initdemands[[j]]$paths)) {
        d <- initdemands

        ## remove the current path
        cat("removing path ",names(d[[j]]$paths[i]))
        #cat(" flow= ",d[[j]]$paths[[i]])
        d[[j]]$paths[i] <- NULL
        ##cat(names(d[[j]]$paths))
        res <- rg.max.concurrent.flow.capacity.restricted.c(g,d,e=e)

        gamma <- rg.mcf.find.gamma(res$gflow,res$lambda)
        
        #metric <- res$lambda * res$betar / initdemands[[j]]$paths[[i]]
        metric <- gamma
        cat(" metric=",metric," res$betar=",res$betar,
            " res$lambda=",res$lambda,"\n")
        demandpaths[[j]][[i]]$lambda <- res$lambda
        demandpaths[[j]][[i]]$beta <- res$betar

        if ( metric > maxmetric ) {
          maxmetric <- metric
          maxmetrici <- i
          demandpathsbest[[j]] <- i
        }

      }
      if(maxmetric > overallMaxMetric ) {
        cat("better metric by ",overallMaxMetric - maxmetric,"\n")
        overallMaxMetric <- maxmetric
        overallMaxMetricD <- j
      } 
    } else {
      cat("  skipping, only one path",names(initdemands[[j]]$paths),"\n")
    }


    ## summarise results
    cat("max metric=",maxmetric," in path",
        names(initdemands[[j]]$paths[maxmetrici]),
        " path number ",maxmetrici)
    if( identical(names(initdemands[[j]]$paths[maxmetrici]),
                  names(intdemands[[j]]$paths[1]))) {
      cat(" !! in intdemands\n") }
    else cat("\n")
    
  }

  bestpathtoremove <- demandpathsbest[[overallMaxMetricD]]
  cat("----- summary-----------------------------------------------\n")
  cat("overall best metric is for demand ",overallMaxMetricD,", path ",
      bestpathtoremove,"\n")

  nametoremove <- names(builtdemands[[overallMaxMetricD]]$paths[bestpathtoremove])

  builtdemands[[overallMaxMetricD]]$paths[bestpathtoremove] <- NULL
  
  cat(" removing path ",nametoremove,"\n")
  cat("  integer demand is ",names(intdemands[[overallMaxMetricD]]$paths[1]),"\n")
  if( identical(names(intdemands[[overallMaxMetricD]]$paths[1]),
                nametoremove)) {
    cat("ERROR removing INTEGER DEMAND\n")
  }
  return(builtdemands)
}

## note demands must be the result of running
##     initresults <- rg.fleischer.max.concurrent.flow.c(g,demands,e=e)
## demands <- initresults$demands
rg.explore2 <- function(g,initdemands,e=0.1,progress=FALSE,intdemands=NULL) {
  ## calc split flows

  builtdemands <- initdemands
  
  ## this line not strictly needed
  #res <- rg.fleischer.max.concurrent.flow.restricted.c(g,initdemands,
  #                                                   updateflow=FALSE,e=e)

  #cat("res$beta=",res$beta," res$betar=",res$betar,
  #    " res$lambda=",res$lambda,"\n")

  ## try one path at a time in each demand and calculate lambda, record maximum

  overallMaxLambda <- -Inf
  overallMaxLambdaD <- 0
  overallMaxLambdaBeta <- -Inf
  bestSumBetaLambda <- -Inf
  bestSumBetaLambdaD <- 0
  worstFlow <- Inf
  worstFlowD <- 0
  worstFlowi <- 0

  demandpaths <- list()
  demandpathsbest <- list()
  j <- 1
  for(d in initdemands) {
    demandpaths[[j]] <- list()
    demandpathsbest[[j]] <- 0
    k <- 1
    for(p in names(d$paths)) {
      demandpaths[[j]][[k]] <- list()
      demandpaths[[j]][[k]]$name <- p
      demandpaths[[j]][[k]]$lambda <- 0.0
      demandpaths[[j]][[k]]$beta <- 0.0
      k <- k + 1
    }
    j <- j + 1
  }

  res <- rg.max.concurrent.flow.capacity.restricted.c(g,initdemands,e=e)
  cat(" res$betar=",res$betar,
      " res$lambda=",res$lambda,"\n")

  d <- res$demands
  
  for(j in 1:length(initdemands)) {
    minlambda <- Inf
    maxlambda <- -Inf
    betaAtMaxLambda <- -Inf
    maxbeta <- -Inf
    minbeta <- Inf
    maxlambdai <- 0
    maxbetai <- 0
    cat("demand number",j," num paths ",length(initdemands[[j]]$paths),"---------------\n")
    if (length(d[[j]]$paths) > 1 ) {
      for(i in 1:length(d[[j]]$paths)) {
        cat(names(d[[j]]$paths[i])," flow=",
            d[[j]]$paths[[i]],"\n")
        if (worstFlow > d[[j]]$paths[[i]] / d[[j]]$flow ) {
          worstFlow <- d[[j]]$paths[[i]] / d[[j]]$flow
          worstFlowD <- j
          worstFlowi <- i
          cat("NEW worst flow\n")
        }
      }
    } else {
      cat("  skipping, only one path",names(d[[j]]$paths),"\n")
    }
    
  }
  
  bestpathtoremove <- names(d[[worstFlowD]]$paths[worstFlowi])
  cat("----- summary-----------------------------------------------\n")
  cat("overall worst flow is for demand ",worstFlowD,", path ",
      bestpathtoremove,"\n")



  builtdemands[[worstFlowD]]$paths[bestpathtoremove] <- NULL
  
  cat(" removing path ",bestpathtoremove,"\n")
  cat("  integer demand is ",names(intdemands[[worstFlowD]]$paths[1]),"\n")
  if( identical(names(intdemands[[worstFlowD]]$paths[1]),
                bestpathtoremove)) {
    cat("ERROR removing INTEGER DEMAND\n")
  }
  return(builtdemands)
}

### Same as rg.fleischer.max.concurrent.flow.restricted()
### except that it tests to see if the integerdemands are in
### the demand paths used in a single phase.
### When running this on a small example it gave the suprising result
### that the integer solution paths were never used in the phase.
### Out of the eight solution paths only seven were ever the same! This
### was for seven out of the many (40000?) phases.
rg.fleischer.max.concurrent.flow.restricted.test <- function(g,demands,e=0.1,updateflow=TRUE, progress=FALSE,integerdemands=NULL,bestgamma=NULL) {

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
  M <- length(vdemands)
  edgeMap <- adjMatrix(g)

  demandpaths <- list()
  demandpathflows <- list()
  integerdemandsindex <- list()
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
      if( identical(p, names(integerdemands[[j]]$paths[1]))) {
        integerdemandsindex[[j]] <- k
      }
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

  
  doubreq <- 2 * ceiling(1/e * log(m/(1-e),base=(1+e)))
  ## updatepb used to measure progress as percent
  updatepb <- as.integer(ceiling(doubreq / 100.0))

  if(progress != FALSE)
    pb <- txtProgressBar(min = 0,
                        max = doubreq, style=3)
  phases <- 0
  totalphases <- 0
  countallint <- 0
  countgamma <- 0
  D <- calcD()

  pdemands <- c(length(vdemands),1:(length(vdemands) -1 ))
  demseq <- 1:length(vdemands)

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

    allinteger <- 0
    demseq <- demseq[pdemands]
    vcaptmp <- vcapacity
    vweights <- rep(0.0,L)
    underdemands <- TRUE

    tmppaths <- list()
    
    for(i in demseq) {
      demand <- vdemands[i]
      D <- calcD()
      while( D < 1 && demand > 0 ) {
        p <- findshortestpath(demandpaths[[i]],vlength)
        ##cat(integerdemandsindex[[i]]," p=",p,"\n")
        tmppaths[[i]] <- as.integer(p)
        if ( identical(as.integer(p),as.integer(integerdemandsindex[[i]])) &&
            allinteger >= 0 ) {
          allinteger <- allinteger + 1
        }
        vweights[ demandpaths[[i]][[p]] ] <- vweights[ demandpaths[[i]][[p]] ] +
          demand
        caponpath <- vcaptmp[ demandpaths[[i]][[p]] ]
        mincap <- min(demand,caponpath)
        if(mincap < demand)
         underdemands <- FALSE
        demand <- demand - mincap
        lengths <- vlength[ demandpaths[[i]][[p]] ]
        lengths <- lengths * (1 + (e*mincap) / vcapacity[ demandpaths[[i]][[p]] ])
        vlength[ demandpaths[[i]][[p]] ] <- lengths
        demandpathflows[[i]][[p]] <- demandpathflows[[i]][[p]] + mincap
        vdemandflow[i] <- vdemandflow[i] + mincap
        #vcaptmp[ demandpaths[[i]][[p]] ] <-
        #  vcaptmp[ demandpaths[[i]][[p]] ] - mincap
        
        
        D <-  calcD()
      }
      
    }
    gamma = min( 1 - vweights/vcapacity )
#    cat("gamma=",gamma,"\n")
    if(underdemands && D < 1) {
      if(allinteger == M) {
        countallint <- countallint +1
      }
      if(gamma >= bestgamma * 0.99999) {
        countgamma <- countgamma + 1

      }
    }
#    cat("NUM INTEGER PATHS",allinteger,"\n")

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
                 lambda=lambda,phases=totalphases,e=e,vlength=vlength,
                 countallint=countallint,countgamma=countgamma)

  if(progress != FALSE)
    close(pb)

  retval
}


rg.max.concurrent.flow.capacity.restricted.c <- function(g,demands,e=0.1,updateflow=TRUE,progress=FALSE) {

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

  doubreq <- 2 * ceiling(1/e * log(m/(1-e),base=(1+e)))

  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = doubreq, style=3)
  } else {
    pb <- NULL
  }

  retlist <- .Call("rg_max_concurrent_flow_capacity_restricted_c",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   e,
                   progress,
                   pb,
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
                 lambda=lambda,phases=retlist$totalphases,e=e,vlength=retlist$vlengths)

  
  return(retval)
}


### Compute the number of combinations of selecting paths from two path sets
### input
### Given c ways of chosing alternative paths from another path set size m
### there are ( n )
###           ( m )    ways (Binomial coeficient)
### so for m = 0 ... n numbers chosen (choose none, then two etc) we have
### number of comb = Sum_0^n
rg.num.path.comb <- function(n) {
  sum <- 0
  for(m in 0:n)
    sum <- sum + choose(n,m)
  return(sum)
}

rg.dual.primal.ratio <- function(numedges,e,eInternal=NULL)  {
  if(is.null(eInternal))
    eInternal <- e
  
  delta <- (numedges / (1-e)) ^ (-1/e)

  ratio <- e * log(1/delta,1+e) / ( (1-e)*log((1-e)/(numedges*delta)))
  return(ratio)

}

rg.run.reed <- function(scenario) {
  g <- scenario$g
  demands <- scenario$commr
  results <- rg.minimum.congestion.flow(g,demands,e=0.02,progress=TRUE,permutation="lowest")
  res <- rg.max.concurrent.flow.int.c(g,results$demands,e=0.02,progress=TRUE,permutation="lowest")
  cat("Gamma = ",res$bestgamma,"\n")
}

rg.run.sp <- function(scenario) {
  g <- scenario$g
  demands <- scenario$commr
  res <- rg.sp.max.concurrent.flow(g,demands)
  cat("Gamma = ",res$gamma,"\n")
}

rg.run.alice <- function() {
  for(i in 10:20) {
    s <- rg.gen.random.scenario(N=i,L=2*i,target=0.2)
    rg.run.reed(s)
    rg.run.sp(s)
  }

}

rg.max.int.flow <- function(g,demands,e=0.1,progress=FALSE,
                            permutation="fixed",deltaf=1.0) {
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
  
  res.2opt <- rg.fleischer.max.concurrent.flow.c(g,demands,e=e2,
                                                 updateflow=FALSE,
                                                 progress=progress,permutation,
                                                 deltaf)
  
  ## now have 2-opt solution as beta_est so  beta < beta_est < 2*beta
  ## scaling by beta_est/2 give best demand for reduced running time
  scalef <- res.2opt$beta /2
  opt2scale <- scalef
  
  demands <- rg.rescale.demands(demands,scalef)


  if(progress != FALSE)
    progress <- "Main calculation"

  res <- rg.fleischer.max.concurrent.flow.with.int.c(g,demands,e=e,
                                                     progress=progress,updateflow=FALSE,
                                                     permutation,
                                                     deltaf)
  
  for(n in names(demands)) {
    res$demands[[n]]$demand <- savedemands[[n]]$demand
  }
  
  res$lambda <- calcLambda(res$demands)
    
  res$beta <- calcBeta(res$demands,res$gdual)
    
  res$estlambdascale <- estlambdascale
  res$opt2scale <- opt2scale
  res
}


rg.fleischer.max.concurrent.flow.with.int.c <- function(g,demands,e=0.1,updateflow=TRUE,progress=FALSE,permutation="lowest",deltaf=1.0) {


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
  retlist <- .Call("rg_fleischer_max_concurrent_flow_with_int_c",
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
                   permutation,
                   deltaf
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
  delta <- deltaf * (ne / (1-e)) ^ (-1/e) 

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

### testing prescaled rg.max.concurrent.flow.int Ran this 20150920 -
### it showed that the prescaled is much faster (mainly for high
### lambda) for very similar quality. Gamma was the same. They both
### had slightly different number of blocked demands when lambda was
### low (lambda<1.0). On average (over 50 results with blocking) the
### ratio of (prescaled-nonscaled)/prescaled was 0.014 (ie nonscaled
### performed very slightly better but there was only 1% in it,
### sometimes nonscaled was better), this was not a statistically
### significant result (p-value 0.08)
rg.max.concurrent.flow.int.test <- function() {
    Nrange <- seq(10,100,10)
    results <<- data.frame(N=numeric(),
                          Lambda=numeric(),
                          Blocked=numeric(),
                          Gamma=numeric(),
                          Runtime=numeric(),
                          Type=character())
    for(N in Nrange) {
        g <- rg.generate.random.graph(N)
        demands <- rg.gen.demands(g,val=round(runif(5,1,4)))
        numDemands <- length(demands)
        g <- rg.set.capacity(g,1)
        linres <- rg.minimum.congestion.flow(g,demands)
        ##linres <- rg.minimum.congestion.flow(g,demands)
        scales <- seq(0.1,2.0,0.2)*linres$lambda
        for(scale in scales) {
            scaled <- rg.rescale.demands(demands,scale)
            timing <- system.time(res <- rg.max.concurrent.flow.int(g,scaled))
            blocked <- numDemands - rg.count.accepted.demands(res$demands)
            interim <- data.frame(N=N,
                                  Lambda=res$lambda,
                                  Gamma=res$gamma,
                                  Blocked=blocked,
                                  Runtime=timing[["user.self"]],
                                  Type="Prescaled")
            print(interim)
            results <<- rbind(results,interim)
            timing <- system.time(res <- rg.max.concurrent.flow.int(g,scaled,
                                                                    prescaled=FALSE))
            blocked <- numDemands - rg.count.accepted.demands(res$demands)
            interim <- data.frame(N=N,
                                  Lambda=res$lambda,
                                  Gamma=res$gamma,
                                  Blocked=blocked,
                                  Runtime=timing[["user.self"]],
                                  Type="NoScaling")
            print(interim)
            results <<- rbind(results,interim)
        }
    }
    return(results)
}
