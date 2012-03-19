## Adjusted from qtlview:::plot.cistrans() author Brian Yandell ##
## pos.peaks = peaks element from peaks() function ##
## ... are additional parameters to plot() and points() ##
## Colors are divided by equal quantiles of the LOD scores ##
## col.legend adds a color legend values are TRUE, FALSE, or a
## vector of the form c(xlim1,xlim2) as a proportion of 
## the plotting area

plotpeaks <- function(pos.peaks, map=NULL, n.col=256,
     cbreaks=quantile(pos.peaks$peaks.lod,probs=seq(0,1,length=n.col+1)),
     chr.peaks=levels(factor(pos.peaks$peaks.chr,exclude=NULL)),
     chr.trait=levels(factor(pos.peaks$trait.chr,exclude=NULL)),
     clab.peaks=chr.peaks,clab.trait=chr.trait,
     xlab="Chromosome of Peak Score",
     ylab="Chromosome of Transcript",
     col.legend=TRUE, lims.legend=c(-0.08,0.3),
     q.legend=c(0.025,0.25,0.5,0.75,0.975),
     ...){

  ## Checks ##
  if(length(cbreaks) != n.col+1) stop("cbreaks must have length n.col+1")
  for(i in c("peaks","trait")){
    if(!all(get(paste("chr",i,sep=".")) %in%
       levels(pos.peaks[[paste(i,"chr",sep=".")]]))){
       warning(paste("chr.",i," not in range",sep=""))
       assign(paste("chr",i,sep="."),
       get(paste("chr",i,sep="."))[
          !is.na(match(levels(pos.peaks[[paste(i,"chr",sep=".")]]),
                get(paste("chr",i,sep="."))))])
     }
    if(is.numeric(get(paste("chr",i,sep="."))))
      assign(paste("chr",i,sep="."),as.character(get(paste("chr",i,sep="."))))
    if(any(is.na(match(get(paste("clab",i,sep=".")),get(paste("chr",i,sep="."))))))
      stop(paste("clab.",i," must be a subset of chr.",i,sep=""))
  }
  if(any(is.na(pos.peaks$trait.pos))){
    warning("Setting trait.pos==NA to trait.pos=0; position in alternate units may be used for plotting purposes")
    pos.peaks$trait.pos[is.na(pos.peaks$trait.pos)] <- 0
  }

  ## Chr lengths equal for peak and trait pos ##
  ## Trait may have chromosomes not in genetic map ##
  len.pos <- function(axis.name="peaks"){
    if(is.null(map) | axis.name != "peaks")
      tapply(pos.peaks[[paste(axis.name,"pos",sep=".")]],
             pos.peaks[[paste(axis.name,"chr",sep=".")]],
             function(x) max(x,na.rm=TRUE))
    else
      sapply(map,function(x) max(x))
  } 
  pklen <- len.pos(axis.name="peaks")
  prlen <- len.pos(axis.name="trait")
  len <- apply(cbind(pklen,prlen[names(prlen) %in% names(pklen)]),1,
         function(x){if(all(is.na(x))) NA else max(x,na.rm=TRUE)})

  ## Subset of Chromosomes ##
  pos.peaks <- pos.peaks[pos.peaks$peaks.chr %in% chr.peaks &
                   pos.peaks$trait.chr %in% chr.trait &
                   !is.na(pos.peaks$trait.chr),]
  if(nrow(pos.peaks)==0) stop("No peaks in chr.peaks & trait.chr")

  ## Allow for different number of chromosomes ##
  axis.pos <- function(axis.name="peaks",len){
    lchr <- get(paste("chr",axis.name,sep="."))
    llen <- ifelse(len==0,5,len)[lchr]
    cumpos <- c(0,cumsum(llen))
    names(cumpos) <- c(names(cumpos)[-1],"max")
    pchr <- as.character(pos.peaks[[paste(axis.name,"chr",sep=".")]])
    ppos <- (pos.peaks[[paste(axis.name,"pos",sep=".")]]-
             llen[pchr]/2)*0.9+llen[pchr]/2+cumpos[pchr]
    list(len=llen,cumpos=cumpos,chr=pchr,pos=ppos)
  }
  ## Positions for Peaks ##
  pk <- axis.pos(axis.name="peaks",len=len)
  ## Positions for Transcript ##
  pr <- axis.pos(axis.name="trait",
     len=c(len,prlen[!(names(prlen) %in% names(pklen))]))

  ## CIS / Trans colors ##
  cols.trans <- rev(rainbow(n.col, start = 0, end = 2/3))
  cols.cis <- rev(gray(seq(0, 0.8, len = n.col)))
  lods.cat <- cut(pos.peaks$peaks.lod,breaks=cbreaks,
        labels=cbreaks[-length(cbreaks)]+0.5*diff(cbreaks),
        include.lowest=TRUE)
  lods.ord <- order(pos.peaks$cis,pos.peaks$peaks.lod)

  ## Plot ##
  par(xpd=TRUE,mgp=c(1.2,0.5,0),las=1)
  plot(pk$cumpos[c(1,length(pk$cumpos))],
       pr$cumpos[c(1,length(pr$cumpos))],
       xlab=xlab,ylab=ylab,
       xaxt="n",type="n",xaxs="i",yaxs="i",yaxt="n",
       ...)
  ## Grid lines ##
  axis(side=1,at=pk$cumpos[-c(1,length(pk$cumpos))],
       labels=FALSE,tck=1,col.ticks="lightgray")
  axis(side=2,at=pr$cumpos[-c(1,length(pr$cumpos))],
       labels=FALSE,tck=1,col.ticks="lightgray")
  ## Axis Labels ##
  axis(side=1,at=(pk$len/2+pk$cumpos[-length(pk$cumpos)]),
       labels=clab.peaks[match(names(pk$len),clab.peaks)], ...)
  axis(side=2,at=(pr$len/2+pr$cumpos[-length(pr$cumpos)]),
       labels=clab.trait[match(names(pr$len),clab.trait)], ...)
  points(pk$pos[lods.ord],pr$pos[lods.ord],
         col=ifelse(pos.peaks$cis[lods.ord]==1,
           cols.cis[as.numeric(lods.cat)][lods.ord],
           cols.trans[as.numeric(lods.cat)][lods.ord]),
           ...)
  if(col.legend==TRUE){
  ## lims is expressed as a fraction of the side ##
    a <- par("usr")
    lims <- c(lims.legend,-0.12,-0.14)
    dx <- a[2]-a[1]
    dy <- a[4]-a[3]
    b <- c(a[1]+lims[1:2]*dx,a[3]+lims[3:4]*dy)
    rect(b[1]+((1:n.col)-1)/n.col*(b[2]-b[1]), b[3],
         b[1]+(1:n.col)/n.col*(b[2]-b[1]), b[4],
         col=cols.trans,density=NA)
    rect(b[1]+((1:n.col)-1)/n.col*(b[2]-b[1]),b[4],
         b[1]+(1:n.col)/n.col*(b[2]-b[1]),b[4]+(b[4]-b[3]),
         col=cols.cis,density=NA)
    segments(b[1],b[3],b[2],b[3])
    segments(b[1]+q.legend*(b[2]-b[1]),
             rep(b[3],length(q.legend)),
             b[1]+q.legend*(b[2]-b[1]),
             rep(b[3]-0.5*(b[4]-b[3]),length(q.legend)))
    qi <- round(quantile(as.numeric(levels(lods.cat)),probs=q.legend),2)
    text(b[1]+c(q.legend[1],q.legend)*(b[2]-b[1]),
         b[3]-c(1.5,rep(0.5,length(q.legend)))*(b[4]-b[3]),
         c("LOD",qi),adj=c(0.5,0),
         cex=ifelse(!is.na(match("cex",names(list(...)))),
           list(...)$cex,1)*0.7)
   }
}

