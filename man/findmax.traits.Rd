\name{findmax.traits}
\alias{findmax.traits}

\title{
  List Traits with QTLs at genomic positions
}
\description{
  For a given vector of positions, list traits that map to those positions
}
\usage{
findmax.traits(chr, pos, sigpos.out, maxlod.out, peaks.out, win = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
Character vector of chromosomes, length equal to \code{pos}
}
  \item{pos}{
Numeric vector of positions, on chromosomes listed in \code{chr}
}
  \item{sigpos.out}{
List of positions with significant LOD scores from \code{\link{maxlod.sigpos}} function
}
  \item{maxlod.out}{
Optional list output from \code{\link{maxlod}} function
}
  \item{peaks.out}{
Optional list output from \code{\link{peaks}} function
}
  \item{win}{
Numeric value for window size; trait positions with significant LOD
  scores within +/- 1/2 \code{win} of \code{chr}, \code{pos} are output
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\note{
  Use in conjunction with \code{findmax}
}

\value{
  Output is a \code{data.frame} with values:
  \item{index}{Index of chr/pos vector}
  \item{chr}{Chromosome}
  \item{pos}{Position on Chromosome}
  \item{id}{Trait ID, from \code{cross} phenotypes}
  \item{peak.pos}{Position of significant maximum LOD score, based on
  output from \code{maxlod.sigpos}}
  \item{lod}{if \code{maxlod} dataset is provided, maximum LOD score for
  traits is included}
  \item{cis}{if \code{peaks} dataset is provided, identifier for CIS
  trait is provided: cis=1 is CIS trait, based on \code{peaks}
  definition, 0 if trans, and NA if unknown}

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Karl W. Broman and Aimee Teo Broman
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(multitrait)
multitrait <- calc.genoprob(multitrait, step=2)

scan.multi <- scanone(multitrait,pheno.col=1:nphe(multitrait))
maxlod.multi <- maxlod(scan.multi)
sigpos.multi <- maxlod.sigpos(maxlod.multi,sig.lod=4)

## Random position ##
findmax.traits(chr=c("5","1"),pos=c(36,87),sigpos.multi)

## Use findmax() ##
n.multi <- neqtl(sigpos.multi,chr=scan.multi[,1],pos=scan.multi[,2],win=5)
x <- findmax(n.multi)

findmax.traits(sigpos.multi,chr=x$chr,pos=x$pos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
