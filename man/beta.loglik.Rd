\name{beta.loglik}
\alias{beta.loglik}
\title{Log-Likelihood of Beta}
\description{Log-likelihood of beta given a pedigree and phenotypic data. Beta is the parameter vector for the multinomial log-linear model. Intended to be used within the function \code{\link{MLE.beta}}}
\usage{
beta.loglik(X, dam_pos=NULL, sire_pos=NULL, par_pos=NULL, beta=NULL,
   beta_map=NULL, merge=NULL, mergeN=NULL, nUS=c(0,0), shrink=NULL)
}
\arguments{
  \item{X}{list of design matrices for each offspring.  Each element should either have dam (D) and/or sire (S) matrices, or a composite Dam/Sire (DS) matrix. See \code{\link{varPed}} for model types}
  \item{dam_pos}{position of each offspring's mother in the dam design matrix}
  \item{sire_pos}{position of each offspring's mother in the sire design matrix}
  \item{par_pos}{position of each offspring's parents in the composite dam/sire matrix}
  \item{beta}{parameter vector}
  \item{beta_map}{vector that maps \code{beta} onto the design matrices (see \code{\link{getXlist}})}
  \item{merge}{optional vector that indicates columns of for which the parameter is transformed using the argument \code{merge} in varPed}
  \item{mergeN}{optional list of matrices for each offspring the columns of which refer to merged variables and the rows to the number of individuals that fall into each category defined by \code{merge})}
  \item{nUS}{vector of the number of unsampled females and males, respectively. Only required if unsampled individuals have known phenotype.}
  \item{shrink}{optional scalar for the variance defining the ridge-regression likelihood penalisation.}
}
\value{
  log-likelihood of beta given the pedigree and X.
}
\references{
Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31
Smouse P.E. \emph{et al} (1999) Journal of Evolutionary Biology 12 1069-1077
              }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{Intended to be used within \code{\link{MLE.beta}}}
\seealso{\code{\link{MLE.beta}}, \code{\link{MCMCped}}, \code{\link{varPed}}, \code{\link{getXlist}}}
\examples{
\dontrun{
data(WarblerP)
data(WarblerG)

GdP<-GdataPed(WarblerG)

res1<-expression(varPed("offspring", relational=FALSE, restrict=0))
var1<-expression(varPed(c("lat", "long"), gender="Male", 
  relational="OFFSPRING"))
res2<-expression(varPed("terr", gender="Female", relational="OFFSPRING", 
  restrict="=="))

PdP<-PdataPed(formula=list(var1,res1,res2), data=WarblerP)

# probability of paternity is modelled as a function of  distance 

X.list<-getXlist(PdP=PdP, GdP=GdP)

ped<-MLE.ped(X.list)$P

# get ML pedigree from genetic data alone

X<-lapply(X.list$X, function(x){list(S=x$XSs)})

# Extract Design matrices for Sires 

sire_pos<-match(ped[,3][as.numeric(names(X))], X.list$id)
sire_pos<-mapply(function(x,y){match(x, y$sire.id)}, sire_pos, X.list$X)

# row number of each design matrix corresponding to the ML sire. 

beta<-seq(-0.065,-0.0325, length=100)
beta_Loglik<-1:100
  for(i in 1:100){
     beta_Loglik[i]<-beta.loglik(X, sire_pos=sire_pos, beta=beta[i], 
     beta_map=X.list$beta_map) 
  }

plot(beta_Loglik~beta, type="l", main="Profile Log-likelihood for beta")
}
}
\keyword{models}
