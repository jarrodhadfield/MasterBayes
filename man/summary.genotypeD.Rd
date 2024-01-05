\name{summary.genotypeD}
\alias{summary.genotypeD}
\alias{summary}
\title{genotypeD Object}
\description{creates and object containing allele and genotype frequency for genotypeD objects}
\usage{
\method{summary}{genotypeD}(object, \dots)
}
\arguments{
  \item{object}{genotypeD object}
  \item{\dots}{other arguments to be passed}
}
\value{
  \item{locus}{locus information field (if present)} 
  \item{allele.names}{vector of allele names: 0 and 1}
  \item{allele.freq}{estimated allele frequencies with finite sample size correction (Lynch and Milligan 1994)}
  \item{genotype.freq}{frequencies of observed genotypes (phenotypes)}
}
\references{Lynch M. and Milligan B.G. (1994) Molecular Ecology 3 91-99}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{genotype}, \code{summary.genotypeD}}
\examples{
\dontrun{
l1<-rbinom(100,1,0.5)
l1<-genotypeD(l1)
summary(l1)
}
}
\keyword{misc}

