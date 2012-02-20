
## annot is a matrix of chr (factor),
## and pos (in same units as maxlod$maxlod.pos)
## rownames(annot) =  names of expression traits

## maxlod.out is output of maxlod() (formerly maxlod and maxlod.pos)
## sigpos.out is output of sigpos() (formerly maxlod5.pos)

peaks <- function(maxlod.out,sigpos.out,
                          annot,win=5){
  annot.chr <- annot[,1]
  annot.pos <- annot[,2]
  names(annot.chr) <- names(annot.pos) <- rownames(annot)
  probe.pos <- split(annot.pos,annot.chr)
## omit probes with unknown location
  maxlod <- lapply(1:nrow(maxlod.out$maxlod),
                   function(i,x) x[i,], maxlod.out$maxlod)
  wh <- lapply(sigpos.out, function(x,a)
        names(x)[names(x) %in% a], names(annot.chr)[!is.na(annot.chr)])
  maxlod.knownloc <- mapply(function(w,m) m[w],wh,maxlod)
  maxlod.pos.knownloc <- mapply(function(w,s) s[w],wh,sigpos.out)

  ## Matches names of positions on same chr if NA then may be on othr chr ##
  maxlod.pos.probe <- mapply(function(pkpos,prpos){
     x <- prpos[match(names(pkpos),names(prpos))]
     names(x) <- names(pkpos)
     x},maxlod.pos.knownloc,probe.pos[names(maxlod.pos.knownloc)])

  dpos <- mapply(function(pkpos,prpos) abs(pkpos-prpos),
                 maxlod.pos.knownloc,maxlod.pos.probe)
  pkpos.cis <- mapply(function(pkpos,d) pkpos[!is.na(d) & d<=win/2],
                 maxlod.pos.knownloc,dpos)
  pkpos.vaguelycis <- mapply(function(pkpos,d) pkpos[!is.na(d)],
                 maxlod.pos.knownloc,dpos)
  pkpos.forsuretrans <- mapply(function(pkpos,d) pkpos[is.na(d)],
                 maxlod.pos.knownloc,dpos)
  pkpos.trans <- mapply(function(pkpos,d) pkpos[is.na(d) | d>win/2],
                 maxlod.pos.knownloc,dpos)
  pklod.cis <- mapply(function(pklod,d) pklod[!is.na(d) & d<=win/2],
                 maxlod.knownloc,dpos)
  pklod.trans <- mapply(function(pklod,d) pklod[is.na(d) | d>win/2],
                 maxlod.knownloc,dpos)

  peaks.unlist <- function(pkpos,pklod)
    data.frame(id=unlist(lapply(pkpos,names)),
       peak.chr=factor(rep(names(pkpos),
         unlist(lapply(pkpos,length))),levels=names(pkpos)),
       do.call(rbind,mapply(function(pos,lod)
               cbind(peak.pos=pos,peak.lod=lod),pkpos,pklod)),
       row.names=NULL)
  peaks.data.frame <- rbind(cbind(peaks.unlist(pkpos.cis,pklod.cis),cis=1),
                    cbind(peaks.unlist(pkpos.trans,pklod.trans),cis=0))
  peaks.data.frame[,paste("probe",c("chr","pos"),sep=".")] <-
                   annot[match(peaks.data.frame$id,rownames(annot)),]
  
  list(peakpos=maxlod.pos.knownloc,
        peaklods=maxlod.knownloc,
        peakpos.cis=pkpos.cis,
        peakpos.vaguelycis=pkpos.vaguelycis,
        peakpos.forsuretrans=pkpos.forsuretrans,
        peakpos.trans=pkpos.trans,
        peaklods.cis=pklod.cis,
        peaklods.trans=pklod.trans,
        peakprobe=maxlod.pos.probe,
        peaks=peaks.data.frame)
}

