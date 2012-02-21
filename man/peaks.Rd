\name{peaks}
\alias{peaks}
\title{
Proximity of Peaks to Transcript Location
}
\description{
Determines distance between position of maximum LOD score and trait
position, based on annotation.  Assigns CIS/Trans status based on a
window size (default=5)
}
\usage{
peaks(maxlod.out, sigpos.out, annot, win = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxlod.out}{
Output from \code{\link{maxlod}} function
}
  \item{sigpos.out}{
Output from \code{\link{maxlod.sigpos}} function
}
  \item{annot}{
    A data.frame consisting of two columns:
    \code{chr} (Chromosome name, as a factor, with levels in the same
      order as the original \code{\link{scanone}} object.  Additional
      levels may be added at the end) and 
    \code{pos} (Position of trait in the same units as
      \code{\link{scanone}} object).
     The rownames should be the trait IDs, same as pheno names in original
       cross object.}
  \item{win}{Window size, in same units as \code{scanone} object
    (e.g. cM or Mb)}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list corresponding to the number of chromosomes:
 \item{peakpos}{List of maximum LOD positions, named by trait, with
   existing location in annotation}
 \item{peaklods}{List of maximum LOD scores, named by trait, with
   existing location in annotation}
 \item{peakpos.cis}{List of maximum LOD positions, named by trait,
   within +/- 1/2 window size of annotation location}
 \item{peakpos.vaguelycis}{List of maximum LOD positions, named by
   trait, on the same chromosome as annotation location}
 \item{peakpos.forsuretrans}{List of maximum LOD positions, named by
   trait, on different chromosome as annotation location}
 \item{peakpos.trans}{List of maximum LOD positions, named by trait,
   outside of +/- 1/2 window size, or on different chromosome of
   annotation location} 
 \item{peaklods.cis}{List of maximum LOD scores, named by trait,
   within +/- 1/2 window size of annotation location}
 \item{peaklods.trans}{List of maximum LOD scores, named by trait,
   outside of +/- 1/2 window size, or on different chromosome of
   annotation location} 
 \item{peakprobe}{List of annotation positions, named by trait}
 \item{peaks}{\code{data.frame} of traits with maximum lod score (by
   chromosome) higher than a significance threshold (as defined in
   \code{\link{maxlod.sigpos}}), columns named "peaks" refer
   to linkage, columns named "probe" refer to annotation.  This object
   can be used for plotting.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Karl W. Broman and Aimee Teo Broman
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link[qtlview]{cistrans}}
}
\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

## Fake annotation ##
annot.bc <- data.frame(chr=factor(c("5","1"),levels=names(fake.bc$geno)),pos=c(20,10))
row.names(annot.bc) <- c("pheno1","pheno2")

scanone.bc <- scanone(fake.bc,pheno.col=1:2)
maxlod.bc <- maxlod(scanone.bc)
sig.bc <- maxlod.sigpos(maxlod.bc,sig.lod=3)

peaks.bc <- peaks(maxlod.bc,sig.bc,annot.bc)

## To plot this in library(qtlview) ##
library(qtlview)
peaks.view <- peaks.bc$peaks

names(peaks.view) <- c("id","peak.chr","peak.pos.cM","peak.score",
                       "cis","trait.chr","trait.pos.cM")
## cM map--if Mb map is available, use this in map list parameter ##
gmap <- pull.map(fake.bc)
plot(cistrans(peaks.view,cross.name="fake.bc",use.cM=TRUE,peak.chr=1:6,
trait.chr=1:6,map=list(Mb.map=gmap,cM.map=gmap)),show.cumscore=FALSE,pch=16,cex=5)
}