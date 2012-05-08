\name{fake.annot}
\alias{fake.annot}
\docType{data}
\title{
 Simulated annotation for fake.f2expr cross
}
\description{
 Annotation data for fake.f2expr expression traits denoting position on
 the genome.
}
\usage{data(fake.annot)}
\format{
  A data.frame with a column for chromosome (chr, as factor), and
  position (pos, as numeric)
}
\details{
  There are 1000 simulated positions, one corresponding to each
  expression trait
}
\examples{
data(fake.annot)
\dontrun{
library(qtl)

# random seed
set.seed(5297830)

# load map
data(map10)

# drop X chr
map10 <- map10[-20]

# chromosome lengths
L <- summary(map10)[1:19,"length"]
names(L) <- as.character(1:19)

# simulate genes in different positions
n.genes <- 1000
chr <- factor(sample(names(L), n.genes, prob=L/sum(L), repl=TRUE), levels=1:19)
pos <- runif(n.genes, 0, L[chr])
genenam <- paste("gene", 1:n.genes, sep="")
genes <- data.frame(chr=chr, pos=pos)
rownames(genes) <- genenam


# which ones have cis?
cis <- sample(c(FALSE,TRUE), n.genes, prob=c(9,1), repl=TRUE)

# trans eqtl on chr 1 at 50 cM
trans <- sample(c(FALSE,TRUE), n.genes, prob=c(8,2), repl=TRUE)
trans[cis] <- FALSE # make cis ones be not trans

# cis effects
cis.eff <- runif(n.genes, 0.2, 4.9)
cis.eff[!cis] <- 0

# trans effects
trans.eff <- runif(n.genes, 0.2, 0.9)
trans.eff[!trans] <- 0

# simulate annotation
chr <- factor(sample(names(L), n.genes, prob=L/sum(L), repl=TRUE), levels=1:19)
pos <- runif(n.genes, 0, L[chr])
fake.annot <- cbind(chr,pos)
rownames(fake.annot) <- genenam
fake.annot[cis,] <- c(genes[cis,"chr"],genes[cis,"pos"]+
                 runif(sum(cis),min=-5,max=5))

save(fake.annot,file="neqtl/data/fake.annot.RData")
}
}
\keyword{datasets}
