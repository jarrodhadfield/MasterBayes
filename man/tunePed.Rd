\name{tunePed}
\alias{tunePed}
\alias{is.tunePed}
\title{tunePed Object}
\description{
An object containing scaling constants for the tuning parameters used in the Metropolis-Hastings updates.  The tuning parameters should be set so that the Metropolis-Hastings acceptance rates lie between 0.2 and 0.5.  Initial tuning parameters for beta and the unsampled population size are obtained from the large sample variance-covariances of the Maximum Likelihood estimates.
}
\usage{
tunePed(E1 = NULL, E2 = NULL, beta = NULL, USdam = NULL, 
   USsire = NULL)
}
\arguments{
  \item{E1}{vector of scaling parameters for E1}
  \item{E2}{vector of scaling parameters for E2}
  \item{beta}{vector which is multiplied by \code{sqrt(10)} to get scaling parameters for beta}
  \item{USdam}{vector which is multiplied by \code{10} to get scaling parameters for the number of unsampled females}
  \item{USsire}{vector which is multiplied by \code{10} to get scaling parameters for the number of unsampled males}
}
\details{
  The proposal distribution for all parameters is the multivariate normal, the variances of which are the large sample  variance covariances of the Maximum Likelihood estimates multiplied by the scaling constants.  For all parameters except \code{beta}, the covariance matrix for the proposal distribution has all off-diagonal elements set to zero.  These parameters must be positive and so the proposal distribution is reflected at zero. A diagonal covariance matrices ensures that the proposal distribution remains symetric.  For beta the covariances are not constrained at zero, and so the matrices are multiplied by the scaling constants in a way that preserves the correlational structure.  The tuning parameters for the error rates are the scaling constants multiplied by 3e-5. 
}
\value{list containing the arguments passed}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
data(WarblerG)
A<-extractA(WarblerG)

ped<-matrix(NA, 100,3)
ped[,1]<-1:100

G<-simgenotypes(A, ped=ped, E1=0.1, E2=0.001, no_dup=2)
GdP<-GdataPed(G=G$Gobs, id=G$id)

model1<-MCMCped(GdP=GdP, nitt=1500, thin=1, burnin=500)

# The proposal distribution is to conservative for E1
# and the update is accepted about 70\% of the time

plot(model1$E1)
autocorr(model1$E1)

# Succesive samples from the posterior distribution are 
# strongly autocorrelated.  Should of course run the chain
# for longer with a larger thinning interval, but a greater
# tuning parameter helps (now 3e-4, rather than 3e-5):

model2<-MCMCped(GdP=GdP, tP=tunePed(E1=10), nitt=1500, 
  thin=1, burnin=500)

plot(model2$E1)
autocorr(model2$E1)
}
}
\keyword{classes}
  
