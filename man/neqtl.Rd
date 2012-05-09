\name{neqtl}
\alias{neqtl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothed count of traits mapping to positions across genome
}
\description{
  Counts the number of positions within a window centered at a position
  on a chromosome.
  
  Given a list of positions of significant LOD scores, counts the number
  of traits with significant maximum LOD scores at a given set of
  positions along the genome
}
\usage{
neqtl(sigpos.out, chr, pos, win = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigpos.out}{
Output from \code{\link{maxlod.sigpos}} function, a list of positions of
significant maximum LOD scores

Any list made up of a vector of positions for each chromosome.
}
  \item{chr}{
Vector of chromosomes (as a factor) for count, corresponding to pos,
e.g. first column of scanone object
}
  \item{pos}{
Vector of positions for count, corresponding to chr, e.g. second
column of scanone object
}
  \item{win}{
Numeric value of window size, count will be number of traits with
significant maximum LOD score within the window, centered at
c(chr,pos)}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{scanone} object with \code{cbind(chr,pos)} as the first two
columns, and number of transcripts with significant maximum LOD score in
the pheno column
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Karl W. Broman and Aimee Teo Broman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[qtl]{scanone}}
  \code{\link[qtl]{plot.scanone}}
}
\examples{
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

## pheno.col=1:nphe(fake.f2expr) takes awhile ##
scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")
maxlod.f2 <- maxlod(scan.f2)
sigpos.f2 <- maxlod.sigpos(maxlod.f2,sig.lod=3)

n.f2 <-neqtl(sigpos.f2,chr=scan.f2[,1],pos=scan.f2[,2],win=5)
plot(n.f2)

}
