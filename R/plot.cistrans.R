## Adjusted from qtlview:::plot.cistrans() author Brian Yandell ##
## pos.peaks = peaks element from maxlod.cistrans() function ##
## ... are additional parameters to plot() and points() ##
## Colors are divided by equal quantiles of the LOD scores ##
## col.legend adds a color legend values are TRUE, FALSE, or a
## vector of the form c(xlim1,xlim2) as a proportion of 
## the plotting area

plot.cistrans <- function(pos.peaks, n.col=256,
     cbreaks=quantile(pos.peaks$peaks.lod,probs=seq(0,1,length=n.col+1)),
     chr.peaks=NULL,chr.trait=chr.peaks,cis.gray=TRUE,
     xlab="Chromosome of Peak Score",
     ylab="Chromosome of Transcript",
     col.legend=TRUE, lims.legend=c(-0.08,0.3),
     q.legend=c(0.025,0.25,0.5,0.75,0.975),...){
  if(length(cbreaks) != n.col+1) stop("cbreaks must have length n.col+1")
  ## Plot subset of Chromosomes ##
  if(!is.null(chr.peaks))
     pos.peaks <- pos.peaks[pos.peaks$peaks.chr %in% chr.peaks &
                  pos.peaks$trait.chr %in% chr.trait &
                  !is.na(pos.peaks$trait.chr),]
  ## Allow for different number of chromosomes ##
  axis.pos <- function(axis.name="peaks"){
    len <- tapply(pos.peaks[[paste(axis.name,"pos",sep=".")]],
           pos.peaks[[paste(axis.name,"chr",sep=".")]],
           function(x) ifelse(length(x)>0,max(x,na.rm=TRUE),0))
    len[is.na(len)] <- 0
    cumpos <- c(0,cumsum(len)[-length(len)])
    pchr <- as.numeric(pos.peaks[[paste(axis.name,"chr",sep=".")]])
    ppos <- (pos.peaks[[paste(axis.name,"pos",sep=".")]]-len[pchr]/2)*0.9+
      len[pchr]/2+cumpos[pchr]
    list(len=len,cumpos=cumpos,chr=pchr,pos=ppos)
  }
  ## Positions for Peaks ##
  pk <- axis.pos(axis.name="peaks")
  ## Positions for Transcript ##
  pr <- axis.pos(axis.name="trait")

  ## CIS / Trans colors ##
  cols.trans <- rev(rainbow(n.col, start = 0, end = 2/3))
  cols.cis <- rev(gray(seq(0, 0.8, len = n.col)))
  lods.cat <- cut(pos.peaks$peaks.lod,breaks=cbreaks,
        labels=cbreaks[-length(cbreaks)]+0.5*diff(cbreaks),
        include.lowest=TRUE)
  lods.ord <- order(pos.peaks$peaks.lod)

  ## Plot ##
  par(xpd=TRUE)
  plot(pk$pos,pr$pos,xaxt="n",mgp=c(1.5,0.5,0),
       xlab="",ylab="",type="n",xaxs="i",yaxs="i",yaxt="n",...)
  mtext(xlab,side=1,line=1.4)
  mtext(ylab,side=2,line=1.4)
  a <- par("usr")
  segments(c(pk$cumpos[pk$len>0][-1],
             rep(a[1],length(pr$len[pr$len>0])-1)),
           c(rep(a[3],length(pk$len[pk$len>0])-1),
             pr$cumpos[pr$len>0][-1]),
           c(pk$cumpos[pk$len>0][-1],
             rep(a[2],length(pr$len[pr$len>0])-1)),
           c(rep(a[4],length(pk$len[pk$len>0])-1),
             pr$cumpos[pr$len>0][-1]),
           col="lightgray")
  segments(c(pk$len[pk$len>0]/2+pk$cumpos[pk$len>0],
             rep(a[1],length(pr$len[pr$len>0]))),
           c(rep(a[3],length(pk$len[pk$len>0])),
             pr$len[pr$len>0]/2+pr$cumpos[pr$len>0]),
           c(pk$len[pk$len>0]/2+pk$cumpos[pk$len>0],
             rep(a[1]-0.005*(a[2]-a[1]),length(pr$len[pr$len>0]))),
           c(rep(a[3]-0.01*(a[4]-a[3]),length(pk$len[pk$len>0])),
             pr$len[pr$len>0]/2+pr$cumpos[pr$len>0]))
  mtext(names(pk$len)[pk$len>0],at=(pk$len[pk$len>0]/2+pk$cumpos[pk$len>0]),
        line=0.3,side=1)
  mtext(names(pr$len)[pr$len>0],at=(pr$len[pr$len>0]/2+pr$cumpos[pr$len>0]),
        line=0.3,side=2,las=2)
  points(pk$pos[lods.ord],pr$pos[lods.ord],
         col=ifelse(pos.peaks$cis[lods.ord]==1 & cis.gray==TRUE,
           cols.cis[as.numeric(lods.cat)][lods.ord],
           cols.trans[as.numeric(lods.cat)][lods.ord]),
         pch=19,cex=0.5, ...)
  if(col.legend==TRUE){
  ## lims is expressed as a fraction of the side ##
    lims <- c(lims.legend,-0.12,-0.14)
    dx <- a[2]-a[1]
    dy <- a[4]-a[3]
    b <- c(a[1]+lims[1:2]*dx,a[3]+lims[3:4]*dy)
    rect(b[1]+((1:n.col)-1)/n.col*(b[2]-b[1]), b[3],
         b[1]+(1:n.col)/n.col*(b[2]-b[1]), b[4],
         col=cols.trans,density=NA)
    if(cis.gray==TRUE)
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
         c("LOD",qi),cex=0.7,adj=c(0.5,0))
   }
}

