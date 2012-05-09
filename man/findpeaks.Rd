\name{findpeaks}
\alias{findpeaks}

\title{
Locates positions of high numbers of eQTL
}
\description{
Given output from \code{neqtl}, finds positions for maximum numbers of
eQTL.  Will find next highest maximum outside of \code{window} from
peak, for \code{n.peaks} number of maximums
}
\usage{
findpeaks(neqtl.out, lodcolumn = 1, window = 5, n.peaks = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{neqtl.out}{
Output from \code{neqtl} function, a \code{scanone} object of number of
eQTL by position within a window
}
\item{lodcolumn}{
    Numeric value for which LOD column to use, default is first column
  }
  \item{window}{
Numeric value of window size, greater than or equal to \code{win}
parameter used in \code{neqtl} function, so neighboring maximums will not
share transcripts
}
  \item{n.peaks}{
Number of peaks to find (e.g. for n.peaks=10, the top 10 positions of
maximum values)
}
}

\value{data.frame of chr and pos of peak positions}

%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{
Karl W. Broman and Aimee Teo Broman
}
%%\note{

%%}

\seealso{
\code{\link{findpeaks.traits}}
}

\examples{
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")
maxlod.f2 <- maxlod(scan.f2)
sigpos.f2 <- maxlod.sigpos(maxlod.f2,sig.lod=3)

n.f2 <-neqtl(sigpos.f2,chr=scan.f2[,1],pos=scan.f2[,2],win=5)
findpeaks(n.f2)

}
