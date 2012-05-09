\name{fake.f2expr}
\alias{fake.f2expr}
\docType{data}
\title{
 Simulated expression traits for F2 intercross
}
\description{
  Simulated expression traits for an F2 intercross, obtained using
  \code{\link[qtl]{sim.cross}}, with \code{\link[qtl]{map10}} as the
  map object.
}
\usage{data(fake.f2expr)}
\format{
  An object of class \code{cross}.  See\code{\link[qtl]{read.cross}} for
  details.
}
\details{
  There are 500 F2 individuals with 1000 simulated expression traits. A
  random 10\% are considered proximal (or CIS) traits, and a random 20\%
  are considered distal or (or trans) traits.
}
\examples{
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
fake.annot <- data.frame(chr=factor(sample(names(L),
           n.genes, prob=L/sum(L), repl=TRUE), levels=1:19),
           pos=runif(n.genes, 0, L[chr]),row.names=genenam)
fake.annot[cis,"chr"] <- genes[cis,"chr"]
fake.annot[cis,"pos"] <-
  apply(cbind(genes[cis,"pos"]+runif(sum(cis),min=-2.5,max=2.5),
        L[genes[cis,"chr"]]),1,function(x) max(0,min(x)))

#save(fake.annot,file="neqtl/data/fake.annot.RData")

# simulate cross; qtl for all cis locations plus the one trans QTL
n.ind <- 500
x <- sim.cross(map10,
               model=rbind(c(1, 50, 0, 0), cbind(as.numeric(chr), pos, 0, 0)),
               type="f2", n.ind=n.ind)
origqtlchr <- c(1, chr)
origqtlpos <- c(50, pos)
o <- order(origqtlchr, origqtlpos)
# qtl get sorted by position so we have to go back
qtlgeno <- x$qtlgeno
qtlgeno[,o] <- qtlgeno # restore to original order
x <- clean(x) # drop qtl geno from object

# simulate phenotypes
pheno <- matrix(rnorm(n.ind*n.genes), ncol=n.genes)

# make them correlated
pheno <- pheno + rnorm(n.ind)

# add effects for cis traits
for(i in which(cis)) 
  pheno[,i] <- pheno[,i] + (qtlgeno[,i+1]-2)*cis.eff[i]

# add effects for trans traits
for(i in which(trans))
  pheno[,i] <- pheno[,i] + (qtlgeno[,1]-2)*trans.eff[i]

# add column names; make data frame; add to cross
colnames(pheno) <- genenam
x$pheno <- as.data.frame(pheno)

fake.f2expr <- x
#save(fake.f2expr,file="neqtl/data/fake.f2expr.RData")
}

data(fake.f2expr)
fake.f2expr <- calc.genoprob(fake.f2expr)
out <- scanone(fake.f2expr, method="hk", phe=1:nphe(fake.f2expr))

}
\seealso{
  \code{\link[qtl]{sim.cross}}
  \code{\link[qtl]{map10}}
  }
\keyword{datasets}
