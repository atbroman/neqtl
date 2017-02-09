
neqtl <- function(sigpos.out, chr, pos, win=5, inc=0.2)
     smoothall(sigpos.out, chr, pos, win, inc)

smoothall <-
function(themax, thechr, thepos, window, inc)
{
  thesmooth <- vector("list", length(themax))
  names(thesmooth) <- names(themax)
  for(i in names(themax))
      thesmooth[[i]] <- smoothchr(themax[[i]], thepos[thechr==i],
                                  window, inc)
   out <- NULL
  for(i in 1:length(thesmooth))
    out <- rbind(out, data.frame(chr=rep(names(themax)[i], nrow(thesmooth[[i]])),
                    pos=thesmooth[[i]][,1], nqtl=thesmooth[[i]][,2]))
  class(out) <- c("scanone", "data.frame")

  rownames(out) <- paste("c", out[,1], ".loc", 1:nrow(out), sep="")

  out
}

## Uses postions from thepos for smoothing: ATB 9/10/09 ##
## In theloc by=0.2 was outside the seq() function--moved it inside  ATB 12/15/09 ##
## removed  u <- unique(theloc), theloc already unique ATB 5/24/12 ##
## added user-defined increment (inc) ATB 2/09/16 ##

smoothchr <-
function(themax, thepos, window, inc)
{
  theloc <- sort(unique(c(thepos, seq(0, max(thepos), by=inc))))
  temploc <- c(themax, theloc)
  tempval <- c(rep(1, length(themax)), rep(0, length(theloc)))
  o <- order(temploc)
  temploc <- temploc[o]
  tempval <- tempval[o]
  smoothed <- runningmean(temploc, tempval,
                 at=theloc,window=window, what="sum")
  return(cbind(theloc, smoothed))
}
