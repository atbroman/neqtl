\name{plotpeaks}
\alias{plotpeaks}

\title{
Plot position of trait vs. position of QTL
}
\description{
Plots trait position, based on annotation against position of QTL, based
on position of maximum LOD score.  Uses output from
\code{peaks} function, CIS traits are shown in grayscale.
}
\usage{
plotpeaks(pos.peaks,
map=NULL, 
n.col = 256,
cbreaks = quantile(pos.peaks$peaks.lod,
          probs = seq(0, 1, length = n.col+ 1)),
chr.peaks = levels(factor(pos.peaks$peaks.chr,exclude=NULL)),
chr.trait = levels(factor(pos.peaks$trait.chr,exclude=NULL)),
clab.peaks = chr.peaks,
clab.trait = chr.trait,
xlab = "Chromosome of Peak Score",
ylab = "Chromosome of Transcript",
col.legend = TRUE,
llims = c(-0.08, 0.3, -0.12, -0.14),
lticks = c(0.025, 0.25, 0.5, 0.75, 0.975),
lcex = 0.8,
 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos.peaks}{\code{peaks} element from \code{peaks}
    function. \code{peaks} is a data.frame.  See Details}
  \item{map}{\code{map} object from cross object, used to determine
    chromosome lengths; if not available, estimated from pos.peaks}
  \item{n.col}{Number of colors to use, default is 256}
  \item{cbreaks}{Breaks to categorize LOD score, length(cbreaks) =
    n.col+1, see Details}
  \item{chr.peaks}{Character vector of subset of peak chromosomes,
    defaults to all chromosomes} 
  \item{chr.trait}{Character vector of trait chromosomes,
    defaults to all chromosomes}
  \item{clab.peaks}{Character vector for x-axis labels, must be a
    subset of chr.peaks}
  \item{clab.trait}{Character vector for y-axis labels, must be a subset
    of chr.trait}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{col.legend}{logical, TRUE to add a color legend along the
    x-axis; colors correspond to LOD scores}
  \item{llims}{numeric vector of the form c(x0,xl,y0,y1) as a
    fraction of the plotting area, where legend should be placed}
  \item{lticks}{numeric vector between 0 and 1 indicating quantiles of
    LOD categories to label in the legend}
  \item{lcex}{magnification of legend text size relative to the
    currrent setting of cex}
  \item{\dots}{additional paramters for \code{points} and \code{text}}
}
\details{
  For plotting purposes, the pos.peaks should have non-missing trait.pos
  for all traits; occasionally, traits will reside on chromosomes that
  are not part of the marker set, and will have unknown cM positions.
  In this case, Mb positions may be substituted to show trans
  relationships.
    
  The default color categories are divided by equally spaced
  quantiles of the peak LOD scores for maximum color range.
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
  \code{\link{peaks}}
  \code{\link[qtlview]{cistrans}} (from qtlview library)
}
\examples{
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")
maxlod.f2 <- maxlod(scan.f2)
sigpos.f2 <- maxlod.sigpos(maxlod.f2,sig.lod=3)

## Fake annotation ##
data(fake.annot)

peaks.f2 <- peaks(maxlod.f2,sigpos.f2,fake.annot)

plotpeaks(peaks.f2$peaks,map=pull.map(fake.f2expr))

}