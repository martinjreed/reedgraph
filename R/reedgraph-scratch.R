analyse.runs <- function(nums=NULL,maxattempts=3) {

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
          retlist <- analyse.inner(scenario,foundratio,file)
          countgamma <- retlist$res$countgamma
          attempts <- attempts + 1
          if(retlist$res$countgamma == 0)
            countzero <- countzero +1
        }
        if(retlist$res$countgamma == 0) {
          neverfound <- neverfound +1
          faillist <- c(file,faillist)
        }

        foundratio <- retlist$foundratio
      }
    } else {
      ##cat("was null\n")
    }
  }
  cat("numzero=",countzero,
      "min=",min(foundratio)," mean=",mean(foundratio)," max=",max(foundratio),
      "\n")
  if(neverfound >0) {
    cat("Did not find",neverfound," for:\n")
    cat(gsub("[^0-9]+([0-9]+)[^0-9]+","\\1",faillist),"\n")
  }
  
  

}

analyse.inner <- function(scenario,foundratio,file) {
  
  results <- rg.minimum.congestion.flow(scenario$g,scenario$commr,e=0.005,progress=TRUE)
  res <- rg.max.concurrent.flow.int.c(
                                      scenario$g,
                                      results$demands,
                                      e=0.005,
                                      scenario=scenario,
                                      progress=TRUE
                                      )
  res$scenario <- scenario
  foundratio <- c(foundratio,as.double(res$countgamma)/as.double(res$phases))
  cat(file,":",res$countgamma,"/",res$phases,
      " gamma=",res$scenario$intgammalist[[1]],
      "best found=",res$bestgamma,"\n")
  retlist <- list()
  retlist$res <- res
  retlist$foundratio <- foundratio
  return(retlist)
}


rg.max.concurrent.flow.int.c <- function(g,demands,e=0.1,updateflow=TRUE,progress=FALSE,scenario) {

  ## note this is not the dual value
  calcD <- function() {
    sum(vcapacity*vlength)
  }
  bestgamma=scenario$intgammalist[[1]]
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
  bestpaths <- as.integer(c())
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
      if(identical(names(demands[[d]]$paths[p]),
                   names(scenario$intdemands[[d]]$paths))) {
        cat("best path in",d,names(scenario$intdemands[[d]]$paths),", ",p,"\n")

        bestpaths[d] <- p
      }
      rdemandpaths[[d]][[p]] <- pathi
      demands.paths <- append(demands.paths,length(pathi))
      demands.paths <- append(demands.paths,pathi)
    }

  }
  print(bestpaths)
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

  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))

  if(progress != FALSE) {
    pb <- txtProgressBar(title = "progress bar", min = 0,
                         max = doubreq, style=3)
  } else {
    pb <- NULL
  }

retlist <- .Call("rg_max_concurrent_flow_int_c",
#  retlist <- .Call("rg_fleischer_max_concurrent_flow_restricted_c_test",
                   demandpaths,
                   vdemands,
                   vcapacity,
                   e,
                   progress,
                   pb,
                   bestgamma,
                   bestpaths-1,
                   environment()
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

  if(TRUE) {
    for(d in 1:length(demandpaths)) {
      cat("demand ",d,": ")
      for(p in 1:length(demandpaths[[d]])) {
        mincap = Inf
        for(e in demandpaths[[d]][[p]]) {
          if(vcapacity[e+1] < mincap) {
            mincap = vcapacity[e+1]
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
        
        #if(bestpaths[[d]] == p)
         # cat("O")
      }
      cat("\n")
    }
  }

  retval <- list(demands=demands,gflow=gflow,gdual=gdual,beta=beta,betar=betar,
                 lambda=lambda,phases=retlist$totalphases,e=e,vlength=retlist$vlengths,
                 countgamma=retlist$countgamma,
                 bestgamma=retlist$bestgamma,
                 bestpaths=retlist$bestpaths+1)

  
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
rg.gen.random.scenario <- function(N=6,L=8) {
  g <- rg.generate.random.graph(N,3,4)
  M <- length(edgeMatrix(g))/2
  comm <- rg.gen.demands(g,L,runif(L,1,5))
  g <- rg.set.capacity(g,runif(M,5,10))
  g <- rg.set.weight(g,0.0)
  results <- rg.minimum.congestion.flow(g,comm,e=0.01)
  ## need to customize this
  commr <- rg.rescale.demands(comm,(1+results$gamma)*0.9)
  results <- rg.minimum.congestion.flow(g,commr,e=0.01)
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

  
  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))
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

  doubreq <- ceiling(2/e * log(m/(1-e),base=(1+e)))

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


