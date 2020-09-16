\name{genotypeD}
\alias{genotypeD}
\alias{is.genotypeD}
\title{genotypeD Object}
\description{Extends the genotype class for dominant marker data}
\usage{
genotypeD(a1, locus=NULL)
}
\arguments{
  \item{a1}{vector of scored genotypes (0 or 1) for dominant markers}
  \item{locus}{object of class locus, gene, or marker, holding information about the source of this genotype.}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{genotype}, \code{summary.genotypeD}}
\examples{
\dontrun{
l1<-rbinom(100,1,0.5)
l1<-genotypeD(l1)
}
}
\keyword{classes}

