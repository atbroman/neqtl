
## annot is a matrix of chr (factor),
## and pos (in same units as maxlod$maxlod.pos)
## rownames(annot) =  names of expression traits

## maxlod.out is output of maxlod() (formerly maxlod and maxlod.pos)
## sigpos.out is output of sigpos() (formerly maxlod5.pos)

maxlod.cistrans <- function(maxlod.out,sigpos.out,
                          annot,win=5){
  annot.chr <- annot[,1]
  annot.pos <- annot[,2]
  names(annot.chr) <- names(annot.pos) <- rownames(annot)
  trait.pos <- split(annot.pos,annot.chr)
## omit traits with unknown location
  maxlod <- lapply(1:nrow(maxlod.out$maxlod),
                   function(i,x) x[i,], maxlod.out$maxlod)
  wh <- lapply(sigpos.out, function(x,a)
        names(x)[names(x) %in% a], names(annot.chr)[!is.na(annot.chr)])
  maxlod.knownloc <- mapply(function(w,m) m[w],wh,maxlod)
  maxlod.pos.knownloc <- mapply(function(w,s) s[w],wh,sigpos.out)

  ## Matches names of positions on same chr if NA then may be on othr chr ##
  maxlod.pos.trait <- mapply(function(mxpos,trpos){
     x <- trpos[match(names(mxpos),names(trpos))]
     names(x) <- names(mxpos)
     x},maxlod.pos.knownloc,trait.pos[names(maxlod.pos.knownloc)])

  dpos <- mapply(function(mxpos,trpos) abs(mxpos-trpos),
                 maxlod.pos.knownloc,maxlod.pos.trait)
  mxpos.cis <- mapply(function(mxpos,d) mxpos[!is.na(d) & d<=win/2],
                 maxlod.pos.knownloc,dpos)
  mxpos.vaguelycis <- mapply(function(mxpos,d) mxpos[!is.na(d)],
                 maxlod.pos.knownloc,dpos)
  mxpos.forsuretrans <- mapply(function(mxpos,d) mxpos[is.na(d)],
                 maxlod.pos.knownloc,dpos)
  mxpos.trans <- mapply(function(mxpos,d) mxpos[is.na(d) | d>win/2],
                 maxlod.pos.knownloc,dpos)
  mxlod.cis <- mapply(function(mxlod,d) mxlod[!is.na(d) & d<=win/2],
                 maxlod.knownloc,dpos)
  mxlod.trans <- mapply(function(mxlod,d) mxlod[is.na(d) | d>win/2],
                 maxlod.knownloc,dpos)

  maxpos.unlist <- function(mxpos,mxlod)
    data.frame(id=unlist(lapply(mxpos,names)),
       peaks.chr=factor(rep(names(mxpos),
         unlist(lapply(mxpos,length))),levels=names(mxpos)),
       do.call(rbind,mapply(function(pos,lod)
               cbind(peaks.pos=pos,peaks.lod=lod),mxpos,mxlod)),
       row.names=NULL)
  peaks <- rbind(cbind(maxpos.unlist(mxpos.cis,mxlod.cis),cis=1),
                    cbind(maxpos.unlist(mxpos.trans,mxlod.trans),cis=0))
  peaks[,paste("trait",c("chr","pos"),sep=".")] <-
          annot[match(peaks$id,rownames(annot)),]
  
  list(maxpos=maxlod.pos.knownloc,
        maxlods=maxlod.knownloc,
        maxpos.cis=mxpos.cis,
        maxpos.vaguelycis=mxpos.vaguelycis,
        maxpos.forsuretrans=mxpos.forsuretrans,
        maxpos.trans=mxpos.trans,
        maxlods.cis=mxlod.cis,
        maxlods.trans=mxlod.trans,
        maxpos.trait=maxlod.pos.trait,
        peaks=peaks)
}

