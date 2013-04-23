## from library(qtlview) ##

## From qtlview:::calc.hc ##
## orders traits by correlation of LOD & position ##
## Remove chr & pos first ##
hc.ord <- function(lods,cluster=TRUE){
   if(cluster & ncol(lods) >= 4) {
     ## Sort by hierarchical cluster ordering.
     rlod <- apply(as.matrix(lods), 2,
             function(x) if(all(x == 0)) x else (x / max(x, na.rm = TRUE)))
     hord <- hclust(dist(t(rlod)))$order
   } else hord <- 1:ncol(lods)
 }

## x is a scanone object ##
## bwd=borderwidth in same units as pos ##
## annot is chr & pos of traits with rownames of traits ##
## cisnam is a vector of trait names that are CIS traits, use
## unlist(sapply(peaks.out$maxpos.cis,names))
image.scanone <- function (x,
    chr = levels(x$chr),
    lodcolumn = 1:min(ncol(x)-2,100),
    n.col = 256,
    allow.neg = FALSE,
    threshold.lod=0,
    cluster=TRUE,
    bwd=0.001,
    annot,
    cisnam, ...){ 
  if(!is.factor(x$chr)) x$chr <- factor(x$chr,levels=unique(x$chr))
  # pull out desired chromosomes
  chr <- qtl:::matchchr(chr, unique(x[,1]))
  
  ## Colors ##
  cols <- rev(rainbow(n.col, start = 0, end = 2/3))
  
  ## Subset first ##
  chrpos <- droplevels(x[x$chr %in% chr,1:2])
  lod <- x[x$chr %in% chr, lodcolumn+2, drop = FALSE]

  ## Check LOD ##
  if (!allow.neg && any(!is.na(lod) & lod < -1e-06)) {
    u <- !is.na(lod) & lod < 0
    n <- sum(u)
    warning(n, " LOD scores <0, set to 0")
    lod[u] <- 0
  }
  if (any(!is.na(lod) & lod == Inf)) {
    u <- !is.na(lod) & lod == Inf
    n <- sum(u)
    warning(n, " LOD scores =Inf, set to NA")
    lod[u] <- NA
  }

  ## Drop traits that don't meet threshold LOD ##
  lod <- lod[,apply(lod,2,function(x,tl)
    ifelse(length(x)>0,ifelse(all(is.na(x)),FALSE,
           max(x,na.rm=TRUE)>=tl),FALSE),threshold.lod),drop=FALSE]
  
  if(ncol(lod) == 0){
    print(paste("No traits with LOD above",threshold.lod))
    return()
  }

  ## Cluster ##
  hc <- hc.ord(lod,cluster)
  ## LOD wrt maxlod ##
  lod <- lod[,hc,drop=FALSE]
  plod <- apply(lod,2,function(x) x/max(x,na.rm=TRUE))

  if(!missing(annot))
    annot <- annot[annot[,1] %in% chr,]
  
  label2 <- colnames(lod)
# label4 <- sprintf("%.1f",apply(lod,2,
#           function(x) round(max(x,na.rm=TRUE),1)))
  label4 <- format(apply(lod,2,
            function(x) round(max(x,na.rm=TRUE),1)))
  
  mai.orig <- mai <- par("mai")
  mai[2] <- mai[4]+max(strwidth(label2,units="inch"))
  mai[4] <- mai[4]+max(strwidth(label4,units="inch"))
  par(mai=mai)

  if(length(chr)>1){
    clen <- tapply(chrpos$pos,chrpos$chr,function(x)
                   diff(range(x,na.rm=TRUE)))
    ccum <- c(0,cumsum(clen+c(rep(bwd,length(clen)-1),0)))
    names(ccum) <- c(names(clen),"end")

    ## Add space between chromosomes ##
    ppos <- do.call(rbind,by(chrpos,chrpos$chr,function(x,d){
                x <- rbind(x,
                data.frame(chr=x$chr[1],pos=x$pos[nrow(x)]+d/2))
                x$pos0 <- x$pos-x$pos[1]
                x},bwd))
    ppos$cumpos <- ppos$pos0+ccum[ppos$chr]
    pplod <- as.matrix(do.call(rbind,by(plod,chrpos$chr,
              function(x) rbind(x,rep(NA,ncol(x))))))
    
    image(ppos$cumpos[-nrow(ppos)], 1:ncol(pplod),
          pplod[-nrow(pplod),], ylab = "", xlab = "Chromosome",
        col = cols, xaxt="n",yaxt = "n") 
    axis(side = 1, at = ccum[2:(length(ccum)-1)]-bwd/2, tck=1,
         labels=FALSE, ...)
    axis(side = 1, at = ccum[-length(ccum)]+(clen+bwd)/2,
       labels=names(clen), ...)
    if(!missing(annot)){
      trnames <- colnames(pplod)[tin <- colnames(pplod) %in% rownames(annot)]
      points(ccum[annot[trnames,1]]+annot[trnames,2]-
          tapply(chrpos$pos,chrpos$chr,function(x) x[1])[annot[trnames,1]],
          (1:ncol(pplod))[tin],pch=1,...)
    }
  } else {
    image(x$pos[x$chr==chr], 1:ncol(lod), plod, ylab = "",
        xlab = paste("Chromosome",chr),
        col = cols, yaxt = "n")
    if(!missing(annot)){
      trnames <- colnames(plod)[tin <- colnames(plod) %in% rownames(annot)]
      points(annot[trnames,2],(1:ncol(lod))[tin],pch=1,...)
    }
  }
  axis(side = 2, at = 1:ncol(lod),
    labels = label2,las = 1, ...)
  if(!missing(cisnam) & !missing(annot)){
      cis <- (colnames(lod) %in% rownames(annot)) &
             (colnames(lod) %in% cisnam)
    axis(side=2, at = (1:ncol(lod))[cis],
         labels=label2[cis],col.axis="red",
         tick=FALSE,lty=0,las=1, ...)
    }
  rmgp <- par("mgp")+
    c(0,max(strwidth(label4,units="inch"))/par("csi"),0)
  axis(side = 4, at= 1:ncol(lod),
    labels = label4 ,las=1, mgp=rmgp,hadj=1, ...)
  axis(side = 4, at=1.04*(ncol(lod)-1),tick=FALSE,
       labels="LOD",las=1,mgp=rmgp,hadj=1,padj=0,xpd=TRUE, ...)

  par(mai=mai.orig)
}

