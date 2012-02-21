\name{nullsims}
\alias{nullsims}

\title{
Null distribution for Maximum LOD score
}
\description{
Simulates a null distribution for the maximum LOD score for a genome
scan of a single QTL model
}
\usage{
nullsims(cross, n.sim = 1000, n.batch = 4, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{n.sim}{Number of simulations to run, default is 1000}
  \item{n.batch}{Number of batches, n.sim / n.batch should be around ~250}
  \item{verbose}{TRUE or FALSE, prints batch number}
  \item{\dots}{Additional parameters for \code{\link[qtl]{scanone}}}}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric vector of maximum LOD scores, of length \code{n.sim}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Karl W. Broman and Aimee Teo Broman
}
\note{
Use \code{\link[stats]{quantile}}(,probs=0.95) to determine LOD threshold
for 5\% false discovery rate

Recommend for \code{n.sim}>1000, code be run in BATCH mode
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[qtl]{scanone}}
  \code{\link[stats]{quantile}}
    
}
\examples{
data(hyper)
# Genotype probabilities
hyper <- calc.genoprob(hyper, step=2.5)

nd <- nullsims(hyper,n.sim=100,n.batch=1,method="hk")

quantile(nd,probs=0.95)
}

\keyword{models}

