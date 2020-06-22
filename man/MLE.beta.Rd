\name{MLE.beta}
\alias{MLE.beta}
\title{Maximum Likelihood Estimation of Beta}
\description{  
Finds MLE for beta given a pedigree, via a call to optim. Beta is the paramater vector of a multinomial log-linear model.}
\usage{
MLE.beta(X.list, ped, beta=NULL, nUSdam=NULL, nUSsire=NULL, shrink=NULL)
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
  \item{ped}{pedigree with id, dam and sire in ech column}
  \item{beta}{optional starting vector for beta}
  \item{nUSdam}{optional number of unsampled females. Only required if unsampled females have known phenotype.}
  \item{nUSsire}{optional number of unsampled males. Only required if unsampled males have known phenotype.}
  \item{shrink}{optional scalar for the variance defining the ridge-regression likelihood penalisation.}
}
\value{
 \item{beta}{vector of MLE's for beta}
 \item{C}{large sample variance-covariance matrix of beta MLE's}
}
\references{
Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31
Smouse P.E. \emph{et al} (1999) Journal of Evolutionary Biology 12 1069-1077
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}, \code{\link{beta.loglik}}}
\examples{
\dontrun{
data(WarblerP)
data(WarblerG)

GdP<-GdataPed(WarblerG)

res1<-expression(varPed("offspring", restrict=0))
var1<-expression(varPed(c("lat", "long"), gender="Male", 
  relational="OFFSPRING"))
res2<-expression(varPed("terr", gender="Female", relational="OFFSPRING",
  restrict="=="))

PdP<-PdataPed(formula=list(var1,res1,res2), data=WarblerP, USsire=FALSE)

X.list<-getXlist(PdP=PdP, GdP=GdP, E2=0.005)

ped<-MLE.ped(X.list)$P
beta<-MLE.beta(X.list, ped) 
beta
}
}
\keyword{models}
\keyword{optimize}
