\name{findpeaks.traits}
\alias{findpeaks.traits}

\title{
  List Traits with QTLs at genomic positions
}
\description{
  For a given vector of positions, list traits that map to those positions
}
\usage{
findpeaks.traits(chr, pos, sigpos.out, maxlod.out, peaks.out, win = 5)
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
  Use in conjunction with \code{findpeaks}
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
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")
maxlod.f2 <- maxlod(scan.f2)
sigpos.f2 <- maxlod.sigpos(maxlod.f2,sig.lod=3)

## Random position ##
findpeaks.traits(chr=c("5","1"),pos=c(36,87),sigpos.f2)

## Use findpeaks() ##
n.f2 <- neqtl(sigpos.f2,chr=scan.f2[,1],pos=scan.f2[,2],win=5)
x <- findpeaks(n.f2)

findpeaks.traits(sigpos.f2,chr=x$chr,pos=x$pos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
