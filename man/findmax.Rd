\name{findmax}
\alias{findmax}

\title{
Locates positions of high numbers of eQTL
}
\description{
Given output from \code{neqtl}, finds positions for maximum numbers of
eQTL.  Will find next highest maximum outside of \code{window} from
peak, for \code{n.max} number of maximums
}
\usage{
findmax(neqtl.out, lodcolumn = 1, window = 5, n.max = 10)
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
  \item{n.max}{
Number of maximums to find (e.g. for n.max=10, the top 10 positions of
maximum values)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Karl W. Broman and Aimee Teo Broman
}
\note{

}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(multitrait)
multitrait <- calc.genoprob(multitrait, step=2)

scan.multi <- scanone(multitrait,pheno.col=1:nphe(multitrait))
maxlod.multi <- maxlod(scan.multi)
sigpos.multi <- maxlod.sigpos(maxlod.multi,sig.lod=4)

n.multi <-
neqtl(sigpos.multi,chr=scan.multi[,1],pos=scan.multi[,2],win=5)

findmax(n.multi)

}
