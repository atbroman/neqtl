\name{maxlod.sigpos}
\alias{maxlod.sigpos}

\title{
Positions of Maximum LOD scores above a Threshold
}
\description{
A list of positions with maximum LOD scores above a threshold; a subset
of the \code{maxlod.pos} output from the \code{maxlod} function.  Mainly
used in subsequent functions
}
\usage{
maxlod.sigpos(maxlod.out, sig.lod = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxlod.out}{
List object, output from the \code{maxlod} function
}
  \item{sig.lod}{
Threshold used to determine significant LOD score, see
\code{\link{nullsims}} function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List with n vectors of m1 positions, where n is the
  number of chromosmes in \code{scanone} and m1 is the number of traits
  with maximum LOD score above a given threshold}
}
\references{
Jansen, Risart ??
}
\author{
  Karl W. Broman and Aimee Teo Broman
  }
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{maxlod}}
}
\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc, step=2.5)

scanone.bc <- scanone(fake.bc,pheno.col=1:2)
maxlod.bc <- maxlod(scanone.bc)
maxlod.sigpos(maxlod.bc,sig.lod=3)

}
