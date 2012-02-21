\name{probes.bypos}
\alias{probes.bypos}

\title{
  List Transcripts with QTLs at genomic positions
}
\description{
  For a given vector of positions, list transcripts that map to those positions
}
\usage{
probes.bypos(chr, pos, sigpos.out, maxlod.out, peaks.out, win = 5)
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
Numeric value for window size; transcript positions with significant LOD
  scores within +/- 1/2 \code{win} of \code{chr}, \code{pos} are output
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Output is a \code{data.frame} with values:
  \item{index}{Index of chr/pos vector}
  \item{chr}{Chromosome}
  \item{pos}{Position on Chromosome}
  \item{id}{Transcript ID, from \code{cross} phenotypes}
  \item{peak.pos}{Position of significant maximum LOD score, based on
  output from \code{peaks}}
  \item{lod}{if \code{maxlod} dataset is provided, maximum LOD score for
  transcripts is included}
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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

scanone.bc <- scanone(fake.bc,pheno.col=1:2)
maxlod.bc <- maxlod(scanone.bc)
sig.bc <- maxlod.sigpos(maxlod.bc,sig.lod=3)

probes.bypos(chr=c("5","2"),pos=c(18,34),sig.bc)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
