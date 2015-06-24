## Simulates the null distribution of max LOD score across genome ##
## ... are parameters to scanone() ##
## n.sim/n.batch ~ 250

nullsims <- function(cross, n.sim=1000, n.batch=4, verbose=TRUE, ...){
  if(n.sim %% n.batch != 0)
    stop("Make batchsize a factor of n.sim.")
  batchsize <- n.sim/n.batch
  np <- qtl::nphe(cross)
  results <- NULL
  for(i in 1:n.batch) {
    if(verbose) print(paste("batch",i))
    cross$pheno[,(1+np):(batchsize+np)] <- rnorm(batchsize*qtl::nind(cross))
    out <- scanone(cross, pheno.col=(1+np):(batchsize+np), ...)
  results.batch <- apply(out[,-(1:2)],2,max,na.rm=TRUE)
  results <- c(results,results.batch)
  }
  results
}
