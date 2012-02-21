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
  where m is the number of chromosomes in \code{scanone} object and m
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
 Output from \code{maxlod} is mainly used in subsequent functions
 \code{neqtl} and \code{peaks}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

scanone.bc <- scanone(fake.bc,pheno.col=1:2)
maxlod.bc <- maxlod(scanone.bc)

}
