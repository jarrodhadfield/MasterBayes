\name{simpedigree}
\alias{simpedigree}
\title{Simulates a Pedigree given a Log-Linear Model}
\description{
  Given a \code{\link{PdataPed}} object simulates a pedigree according to the linear model defined by \code{formula} and user specified parameter values for the given model.
}
\usage{
simpedigree(PdP, beta=NULL, nUS=NULL)
}
\arguments{
 \item{PdP}{a \code{\link{PdataPed}} object}
 \item{beta}{parameter vector for the model defined by the \code{formula} argument in \code{\link{PdataPed}}} 
 \item{nUS}{vector for the size of the unsampled population(s) defined in the \code{USdam} and \code{USsire} arguments passed to \code{\link{PdataPed}}.  Parmeters for the unsampled female population come before the male population.}
}
\value{
   \item{ped}{pedigree in 3 columns: id, dam, sire. Base individuals have NA as parents}
   \item{USsire.data}{binary vector indicating unsampled sire records (1)}
   \item{USsire.formula}{variable of the form \code{expression(\link{varPed}(...))} that can be included in the \code{formula} argument of \code{\link{PdataPed}} so that unobserved male records are effectively hidden}
   \item{USdam.data}{binary vector indicating unsampled dam records (1)}
   \item{USdam.formula}{variable of the form \code{expression(\link{varPed}(...))} that can be included in the \code{formula} argument of \code{\link{PdataPed}} so that unobserved male records are effectively hidden}
}
\references{Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{MCMCped}
\keyword{datagen}
