\name{plotpeaks}
\alias{plotpeaks}

\title{
Plot position of trait vs. position of QTL
}
\description{
Plots trait position, based on annotation against position of QTL, based
on position of maximum LOD score.  Uses output from
\code{maxlod.cistrans} function, CIS traits are shown in grayscale.
}
\usage{
plotpeaks(pos.peaks, map=NULL,n.col = 256, cbreaks = quantile(pos.peaks$peaks.lod, probs = seq(0, 1, length = n.col + 1)), chr.peaks = levels(factor(pos.peaks$peaks.chr,exclude=NULL)), chr.trait = levels(factor(pos.peaks$trait.chr,exclude=NULL)), xlab = "Chromosome of Peak Score", ylab = "Chromosome of Transcript", col.legend = TRUE, lims.legend = c(-0.08, 0.3), q.legend = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos.peaks}{\code{peaks} element from \code{maxlod.cistrans}
    function. \code{peaks} is a \code{data.frame}}
  \item{map}{\code{map} object from cross object, used to determine
    chromosome lengths; if not available, estimated from pos.peaks}
  \item{n.col}{Number of colors to use, default is 256}
  \item{cbreaks}{Breaks to categorize LOD score, length(cbreaks) = n.col+1}
  \item{chr.peaks}{Character vector of subset of peak chromosomes,
    defaults to all chromosomes} 
  \item{chr.trait}{Character vector of trait chromosomes,
    defaults to all chromosomes}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{col.legend}{logical, TRUE to add a color legend to the plot,
    added along the x-axis}
  \item{lims.legend}{numeric vector of the form c(xlim1,xlim2) as a
    proportion of the plotting area, where legend should be placed}
  \item{q.legend}{numeric vector between 0 and 1 indicating LOD categories to
    label in the color legend}
  \item{\dots}{additional paramters for \code{plot} and \code{points}}
}
\details{
  The default color categories are divided by equally spaced
  quantiles of the peak LOD scores for maximum color range. If
  cis.gray=FALSE CIS traits are plotted with the same color as the trans
  traits. 
}
\value{Outputs a plot}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Brian Yandell and Aimee Teo Broman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[qtlview]{cistrans}}
}
\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

## Fake annotation ##
annot.bc <- data.frame(chr=factor(c("5","1"),
            levels=names(fake.bc$geno)),pos=c(20,10))
row.names(annot.bc) <- c("pheno1","pheno2")

scanone.bc <- scanone(fake.bc,pheno.col=1:2)
maxlod.bc <- maxlod(scanone.bc)
sig.bc <- maxlod.sigpos(maxlod.bc,sig.lod=3)

peaks.bc <- maxlod.cistrans(maxlod.bc,sig.bc,annot.bc)

plotpeaks(peaks.bc$peaks)

## Using qtlview ##
#peaks <-peaks.bc$peaks
#names(peaks)[c(2:4,7)] <- c("peak.chr","peak.pos.cM","peak.score","trait.pos.cM")
#plot(cistrans(peaks,cross.name="fake.bc",
#   maps=list(Mb.map=pull.map(fake.bc),cM.map=pull.map(fake.bc))))

}