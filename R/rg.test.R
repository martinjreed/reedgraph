rg.test <- function (g,demands) {
  res <- rg.max.concurrent.flow.prescaled(g,demands,e=0.1,progress=TRUE)
  res2 <- rg.max.concurrent.flow.prescaled(g,demands,e=0.1,progress=TRUE,ccode=FALSE)

  print("using c")
  rg.fleischer.max.concurrent.flow.stats(res)
  print("using r")
  rg.fleischer.max.concurrent.flow.stats(res2)
  
}



