\name{maxlod}
\alias{maxlod}

\title{
Maximum LOD score and Position
}
\description{
Calculates maximum LOD score for each chromosome of a
\code{\link[qtl]{scanone}} object, and position of that max on the chromosome
}
\usage{
maxlod(scan.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{scan.out}{An object of class \code{scanone}. See
   \code{\link[qtl]{scanone}} for details.}
}
%\details{
%}
\value{
Returns list of 
\item{maxlod}{matrix of maximum LOD, made up of n rows and m columns,
  where n is the number of chromosomes in \code{scanone} object and m
  is the number of traits}
\item{maxlod.pos}{List with n vectors of m positions, where n is the
  number of chromosmes in \code{scanone} and m is the number of traits}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Karl W. Broman and Aimee Teo Broman
}
\note{
 Output from is mainly used in subsequent functions
 \code{maxlod.sigpos},\code{peaks}, and \code{findpeaks.traits}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{maxlod.sigpos}}
\code{\link{peaks}}
\code{\link{findpeaks.traits}}
}

\examples{
data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)

scan.f2 <- scanone(fake.f2expr,pheno.col=1:100, method="hk")
maxlod.f2 <- maxlod(scan.f2)
}
