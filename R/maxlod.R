## scan.out = output of scanone() ##

maxlod <- function(scan.out){
  chrs <- levels(scan.out[,1])
# the maximum lod scores on each chromosome (output is dim nchr x ntraits)
  maxlod <-  apply(scan.out[,-(1:2)], 2, tapply, scan.out[,1], max)
# the locations (take average of those locations with the maximum LOD score)
# [list by chromosome, each being a vector of peak positions
  maxlod.pos <- by(scan.out,scan.out$chr,function(x)
       apply(x[,-(1:2)],2,function(a,b)
       mean(b[!is.na(a) & a==max(a, na.rm=TRUE)]),x[,2]))
  attributes(maxlod.pos) <- list(names=
                 attributes(maxlod.pos)$dimnames[[1]])
  list(maxlod=maxlod,maxlod.pos=maxlod.pos)
}

## List of positions by chromosome with max LOD scores above
## a certain threshold
maxlod.sigpos <- function(maxlod.out,sig.lod=5)
   mapply(function(pos,m) pos[m>=sig.lod],
       maxlod.out$maxlod.pos,as.data.frame(t(maxlod.out$maxlod)))


