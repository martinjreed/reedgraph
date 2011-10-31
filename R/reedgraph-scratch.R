### combinations possible with
### k hash functions, m length
rg.bloom.comb <- function(k,m) {

  prod <- 1.0
  for(i in m:(m-k+1)) {
    print(prod)
    prod <- i * prod
  }
  return(prod)
}

### estimate according to Bose 2007, not quite correct
### k hash functions, m length, n items in filter
### this is a strict lower bound for k>=2
rg.bloom.false <- function(k,m,n) {
  return((1-(1-1/m)^(k*n))^k)
}

### optimium k for Bloom filter
rg.bloom.optimum.k <- function(m,n) {
  return(log(2) * m /n)
}

rg.analyse.tests <- function(res) {

  var <- data.frame()
  for(N in unique(res$N)) {
    var <- rbind(var,lapply(res,subset,res$N==N & res$L==res$N * res$N))
  }
  return(var)
}

### Returns a graph with number of flows in each edge
rg.count.edge.flows <- function(g,demands) {
  for(i in 1:length(demands)) {
    for(j in names(demands[[i]]$paths)) {
      ##results$demands[[i]]$paths <- 1.0
      demands[[i]]$paths[[j]] <- 1.0
    }
  }
  gcount <- rg.max.concurrent.flow.graph(g,demands)
  return(gcount)
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



## test for PURSUIT TM theory
##
rg.test.int.versus.nonint.flow <- function(N=8,L=NULL,target=0.0,e=0.1,g=NULL) {
  if(is.null(g)) {
    g <- rg.generate.random.graph(N,2,25)
    
  } else {
    N <- length(nodes(g))
  }
  M <- length(edgeMatrix(g))/2
  if(is.null(L)) {
    ## generate using lognormal traffic as in
    ## Nucci, Antonio and Sridharan, Ashwin and Taft, Nina},
    ## The problem of synthetically generating IP traffic matrices: initial recommendations},
    ## SIGCOMM} Comput. Commun. Rev., vol 35(3) July 2005 pp19-32
    ## this generates a mean of 1 with a variance that matches the paper
    comm <- rg.gen.demands(g,val=rlnorm(N*N,16,1)/exp(16.5))
  } else {
    comm <- rg.gen.demands(g,L,val=rlnorm(L,16,1)/exp(16.5))
  }

  ## set some random capacities
  g <- rg.set.capacity(g,runif(M,5,10))
  g <- rg.set.weight(g,0.0)
  ##results <- rg.sp.max.concurrent.flow(g,comm)
  results <- rg.minimum.congestion.flow(g,comm,e=e,progress=TRUE)
  ## need to customize this
  multiplier <- (1-target)/(1-results$gamma)
  comm <- rg.rescale.demands(comm,multiplier)
  sp.results <- rg.sp.max.concurrent.flow(g,comm)
  runtime <- as.double(system.time(flow.results <- rg.minimum.congestion.flow(g,comm,e=e,progress=TRUE))[1])
  flow.results$runtime <- runtime
  gcount <- rg.count.edge.flows(g,flow.results$demands)
  flow.results$max.mpls <- max(as.double(edgeData(gcount,att="weight")))
  flow.results$max.bf <- max(as.double(lapply(edges(g),length)))
  
  int.results <- rg.max.concurrent.flow.int.c(g,flow.results$demands,e=e,progress=TRUE)

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
      for(u in edges(g)[[n]] ) {
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
      
      for(u in edges(g)[[t]] ) {
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
  for (v in edges(g,u)[[1]]) {
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


### permutation: how the demands are chosen can be either
###              c(....) integers specifying demand order
###              "fixed" done in fixed order
###              "random" done in random order
###              "lowest" done in lowest cost (lowest dual path) order
###

rg.max.concurrent.flow.int.c <- function(g,demands,e=0.1,eInternal=NULL,updateflow=TRUE,progress=FALSE,scenario=NULL,permutation="fixed",deltaf=1.0) {

  ## note this is not the dual value
  calcD <- function() {
    sum(vcapacity*vlength)
  }

  if(is.null(eInternal))
    eInternal <- e

  if(is.numeric(permutation))
    permutation <- permutation - 1
  else if(permutation == "fixed")
    permutation <- 0:(length(demands)-1)
  else if(permutation == "random")
    permutation <- -1
  else if(permutation == "lowest")
    permutation <- -2

  savedemands <- demands
  vdemands <- as.double(lapply(demands,"[[","demand"))
  vcapacity <- as.double(edgeData(g,attr="capacity"))
  vlength <- rep(0.0,length(vcapacity))
  vdemandflow <- rep(0.0,length(vcapacity))
  L <- length(vlength)
  delta <- deltaf * (L / (1-e)) ^ (-1/e)
  vlength <- delta / vcapacity
  edgeMap <- adjMatrix(g)

  demands.paths <- as.integer(c())
  bestpaths <- as.integer(c())
  rdemandpaths <- list();

  if(!is.null(scenario)) {
    bestgamma=scenario$intgammalist[[1]]
    for(d in 1:length(demands)) {
      rdemandpaths[[d]] <- list()
      demands.paths <- append(demands.paths,d-1)
      demands.paths <- append(demands.paths,length(demands[[d]]$paths))
      missing <- TRUE
      
      for(p in 1:length(demands[[d]]$paths)) {
        pathi <-
          as.integer(strsplit
                     (names(demands[[d]]$paths[p]),"|",fixed=TRUE)[[1]])
        pathi <- pathi -1
        if(identical(names(demands[[d]]$paths[p]),
                     names(scenario$intdemands[[d]]$paths))) {
          cat("best path in",d,names(scenario$intdemands[[d]]$paths),", ",p,"\n")
          missing <- FALSE
          bestpaths[d] <- p
        }
        rdemandpaths[[d]][[p]] <- pathi
        demands.paths <- append(demands.paths,length(pathi))
        demands.paths <- append(demands.paths,pathi)
      }
      if(missing)
        cat("best path in",d,"is MISSING!\n")
      
    }
    print(bestpaths)
  } else {
    bestgamma <- Inf
  }
  ## code the paths as integers
  ## [demandno(e.g=1) numdempaths lengthpath1 path1[1] path1[2] ...
  ##  lengthpath2 path2[1] path2[2]....demandno(e.g=2)....]

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

  # old for ori
  # doubreq <- 2 * ceiling(1/e * log(m/(1-e),base=(1+e)))
  doubreq <- 2* ceiling(log(1.0/delta,1+eInternal))
  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = doubreq, style=3)
  } else {
    pb <- NULL
  }

  retlist <- .Call("rg_max_concurrent_flow_int_c",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   e,
                   eInternal,
                   progress,
                   pb,
                   bestgamma,
                   environment(),
                   permutation,
                   deltaf
                   );
  
  if( progress != FALSE) {
    close(pb)
  }

  foundbestpaths <- retlist$bestpaths + 1
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

  if(FALSE) {
    for(d in 1:length(demandpaths)) {
      cat("demand ",d,": ")
      for(p in 1:length(demandpaths[[d]])) {
        mincap = Inf
        for(ed in demandpaths[[d]][[p]]) {
          if(vcapacity[ed+1] < mincap) {
            mincap = vcapacity[ed+1]
          }
        }
        cat(" ",retlist$pathcount[[d]][[p]])
        if(mincap < savedemands[[d]]$demand) {
          cat("u")
        } else {
          cat("o")
        }
        
        if(foundbestpaths[[d]] == p)
          cat("F")
        
        if(identical(names(demands[[d]]$paths[p]),
                     names(scenario$intdemands[[d]]$paths))) {
                                        #if(bestpaths[[d]] == p)
          cat("X")
        }
      }
      cat("\n")
    }
  }

  intdemands <- list()
  bestpaths=retlist$bestpaths+1

  for(i in names(demands)) {
    intdemands[[i]]$source <- demands[[i]]$source
    intdemands[[i]]$sink <- demands[[i]]$sink
    intdemands[[i]]$demand <- demands[[i]]$demand
    intdemands[[i]]$flow <- demands[[i]]$demand
    intdemands[[i]]$paths <- list()
    intdemands[[i]]$paths[[names(demands[[i]]$paths[bestpaths[as.integer(i)]])]] <-
      demands[[i]]$demand

  }

  gflowint <- rg.max.concurrent.flow.graph(gdual,intdemands)
  
  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,betar=betar,
                 lambda=lambda,phases=retlist$totalphases,e=e,vlength=retlist$vlengths,
                 countgamma=retlist$countgamma,
                 bestgamma=retlist$bestgamma,
                 bestpaths=bestpaths,
                 pathdiffcount=retlist$pathdiffcount,
                 phasepathdiffcount=retlist$phasepathdiffcount,
                 gammavals=retlist$gammavals,
                 betavals=retlist$betavals,
                 lambdavals=retlist$lambdavals,
                 intdemands=intdemands,
                 gflowint <- gflowint)

  
  return(retval)
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
