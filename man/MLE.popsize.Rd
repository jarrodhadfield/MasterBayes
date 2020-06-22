\name{MLE.popsize}
\alias{MLE.popsize}
\title{Maximum Likelihood Estimation of the Unsampled Population Size}
\description{  
Finds the MLE for the number of unsampled males and/or females following Nielsen \emph{et al.} (2001). The size of the unsampled population can vary over time and space, and genotyping error is accomodated using the CERVUS model of genotyping error (Kalinwoski \emph{et al.} 2006). 
}
\usage{
MLE.popsize(X.list, USdam=FALSE, USsire=FALSE, nUS=NULL, 
   ped=NULL, shrink=NULL)
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
 \item{USdam}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled females exists. If \code{USdam} is a character vector it must have the same length as the number of offspring (\code{length(X.list$X)}) with factor levels representing sub-populations (in time or space) over which the number of unsampled females vary.}  
  \item{USsire}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled males exists. if    \code{USsire} is a character vector it must either have the same length as the number of offspring (\code{length(X.list$X)}) with factor levels representing sub-populations (in time or space) over which the number of unsampled males vary, or alternatively \code{"USdam"}, in which case the unsampled male and female populations are constrained to be equal.}  
  \item{nUS}{optional starting vector for the size of the unsampled population.  Parmeters for the unsampled female population come before the male population.}
  \item{ped}{optional pedigree with id, dam and sire in ech column}
  \item{shrink}{optional scalar for the variance defining the ridge-regression likelihood penalisation.}
}
\value{
 \item{nUS}{vector of MLE's for the size of the unsampled population. Lower bound is 1e-5 for numerical stability.}
 \item{C}{large sample variance-covariance matrix of \code{nUS} MLE's}
}
\references{Nielsen. R. \emph{et.al} Genetics (2001) 157 4 1673-1682}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{
Nielsen's original model does not account for genotyping error, and estimation of the unsampled population size is VERY sensitive to the level of genotyping error.  This function implements a commonly used approxiamtion for genotyping error that ignores pedigree information. For many problems this approximation seems valid, but appears to break down when estimating the size of the unsampled population size.  Bayesian estimation of the unsampled population size (see \code{\link{MCMCped}}) that uses an exact solution for genotyping error is more robust.
}
\seealso{\code{\link{MCMCped}}, \code{\link{popsize.loglik}}}
\examples{
\dontrun{
data(WarblerP)
data(WarblerG)

GdP<-GdataPed(WarblerG)
res1<-expression(varPed("offspring", restrict=0))

PdP<-PdataPed(formula=list(res1), data=WarblerP, USsire=TRUE, USdam=TRUE)

X.list<-getXlist(PdP=PdP, GdP=GdP, E2=0.02)

nUS<-MLE.popsize(X.list, USsire=TRUE, USdam=TRUE)
nUS
}
}
\keyword{models}
\keyword{optimize}
