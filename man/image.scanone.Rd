\name{image.scanone}
\alias{image.scanone}

\title{
Heatmap of scanone object
}
\description{
Heatmap of scanone object; useful for a subset of expression traits
(<100).
}
\usage{
image.scanone(x,chr = levels(x$chr), lodcolumn = 1:min(ncol(x)-2,100), n.col = 256, allow.neg = FALSE, threshold.lod = 0, cluster = TRUE, bwd = 0.001, annot, cisnam, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"scanone"}, as output by \code{\link[qtl]{scanone}}
}
  \item{chr}{
Optional vector indicating chromosomes to plot; this should be a
character vector referring to chromosomes by name
}
  \item{lodcolumn}{
A vector of integers indicating which of the LOD score columns to plot.
Although the number of columns is suggested to be <100, larger numbers
are allowed.
}
  \item{n.col}{
Number of colors / bins to use for LOD scores
}
  \item{allow.neg}{
Logical, if FALSE, sets negative LOD scores to 0
}
  \item{threshold.lod}{
Logical, default is 0.  If traits do not have LOD >= threshold.lod, they
are excluded from plot.
}
  \item{cluster}{
Logical; if TRUE, hierarchical clustering used to order traits by
similar LOD patterns.  If FALSE, same order as specified in lod.column
}
  \item{bwd}{
Numeric value for space between chromosomes in plot; in same units as
"pos" values in x.
}
  \item{annot}{
    Optional data.frame of trait annotation consisting of two columns:
    \code{chr} (Chromosome name, as a factor, with levels in the same
    order as x$chr) and
    \code{pos} (Position of trait in the same units as x$pos)
    rownames(annot) should be the trait IDs, the same as colnames(x)[-(1:2)].
}
  \item{cisnam}{
Optional character vector of trait names that are defined as CIS; these
will be labelled in red on the plot.  See Details.
}
  \item{\dots}{
Optional parameters to pass to axis(), points() or mtext() functions
}
}
\details{
The \code{cisname} vector can be taken from the output of \code{\link{peaks}}:
Let peaks.out be output from the peaks() function.  Then use either
unlist(sapply(peaks.out$maxpos.cis, names)) or
peaks.out$peaks$id[peaks.out$peaks$cis==1]. 
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Brian Yandell and Aimee Teo Broman
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{peaks}}
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

cisnam <- peaks.f2$peaks$id[peaks.f2$peaks$cis==1]
# unlist(sapply(peaks.f2$maxpos.cis,names))

image.scanone(scan.f2,lodcolumn=1:100,threshold.lod=3,
   cisnam=cisnam,annot=fake.annot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
