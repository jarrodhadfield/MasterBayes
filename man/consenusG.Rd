\name{consensusG}
\alias{consensusG}
\title{Obtains a consensus genotype from duplicate samples}
\description{
A function for obtaining a consensus genotype from duplicate samples. The amount of missing data is minimised, and preference is given to samples with lower genotyping error  
}
\usage{
consensusG(GdP, cat.levels=NULL, gmax=FALSE, het=FALSE)
}
\arguments{
  \item{GdP}{a \code{GdataPed} object}
  \item{cat.levels}{order of genotyping error rate categories, with most reliable category first}
  \item{gmax}{logical; if a most represented genotype exists should it be saved}
  \item{het}{logical; should heterozygotes be saved over homozygotes - overrides \code{cat.levels}}

}
\value{
  \item{GdP}{a \code{GdataPed} object}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{GdataPed}}}
\keyword{manip}
\keyword{models}
