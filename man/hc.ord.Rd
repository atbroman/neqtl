\name{hc.ord}
\alias{hc.ord}

\title{
Order traits by similarity of LOD scores
}
\description{
Traits are ordered using Hierarchical clustering of scanone output by similar
positions of QTL. Can be used outside of \code{imagelods} function to
order traits, if several imagelod() plots need to be made (eg. for
number of traits > 100) 

}
\usage{
hc.ord(lods, cluster = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lods}{
Scanone output with c("chr","pos") columns removed
}
  \item{cluster}{
Logical; if TRUE, order of columns is clustered, otherwise, order is
same as input order.
}
}
%\details{
%}
\value{
  Integer vector of lodcolumn order
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
\code{\link[qtlview]{plot.aug.scanone}}
}
\examples{
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")

hc.ord(scan.f2[-(1:2)],cluster=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
