\name{neqtl}
\alias{neqtl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothed count of traits mapping to positions across genome
}
\description{
Counts the number of traits with significant maximum LOD scores at each
position along the genome
}
\usage{
neqtl(sigpos.out, chr, pos, win = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigpos.out}{
Output from \code{\link{maxlod.sigpos}} function, a list of positions of
significant maximum LOD scores
}
  \item{chr}{
Vector of chromosomes (as a factor) for count, corresponding to \code{pos},
e.g. first column of \code{scanone} object
}
  \item{pos}{
Vector of positions for count, corresponding to \code{chr}, e.g. second
column of \code{scanone} object
}
  \item{win}{
Numeric value of window size, count will be number of traits with
significant maximum LOD score within the window, centered at
\code{c(chr,pos)}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{scanone} object with \code{cbind(chr,pos)} as the first two
columns, and number of transcripts with significant maximum LOD score in
the pheno column
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Karl W. Broman and Aimee Teo Broman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[qtl]{scanone}}
  \code{\link[qtl]{plot.scanone}}
}
\examples{
data(multitrait)
multitrait <- calc.genoprob(multitrait, step=2)

scan.multi <- scanone(multitrait,pheno.col=1:nphe(multitrait))
maxlod.multi <- maxlod(scan.multi)
sigpos.multi <- maxlod.sigpos(maxlod.multi,sig.lod=4)

n.multi <-
neqtl(sigpos.multi,chr=scan.multi[,1],pos=scan.multi[,2],win=5)
plot(n.multi)

}
