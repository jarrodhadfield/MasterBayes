\name{startPed}
\alias{startPed}
\alias{is.startPed}
\title{startPed Object}
\description{An object containing the starting parameterisation of a model, and logical variables indicating wether parameters should be estimated or fixed at the starting parameterisation.  By default the starting parameterisation is obtained through a mixture of Maximum Likelihood and heuristic techniques.}
\usage{
startPed(G=NULL, id=NULL, estG=TRUE, A=NULL, estA=TRUE, E1=NULL, 
   estE1=TRUE, E2=NULL, estE2=TRUE, ped=NULL, estP=TRUE, 
   beta=NULL, estbeta=TRUE,  USdam=NULL, estUSdam=TRUE, 
   USsire=NULL, estUSsire=TRUE, shrink=NULL)
}
\arguments{
  \item{G}{list of genotype objects}
  \item{id}{vector of indivual id's for \code{G}}
  \item{estG}{logical; should genotypes be estimated?}
  \item{A}{list of allele frequencies}
  \item{estA}{logical; should base-population allele frequencies be estimated? }
  \item{E1}{if Wang's (2004) model of genotyping error for co-dominant markers is used this is a vector of probabilities of an allele dropping out. If CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers is used this parameter is not used. If Hadfield's (2009) model of genotyping error for dominant markers is used this is a vector of probabilities of a dominant allele being scored as a recessive allele. Default=0.005.}
  \item{estE1}{logical; should E1 estimated?}
 \item{E2}{if Wang's (2004) or CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers are used this is a vector of probabilities of an allele being miss-scored. In the CERVUS model errors are not independent for the two alleles within a genotype and so if a genotyping error has occurred at one allele then a genotyping error occurs at the other allele with probability one. Accordingly, \code{E2}(2-\code{E2}) is the per-genotype rate defined in CERVUS. If Hadfield's (2009) model of genotyping error for dominant markers is used this is a vector of probabilities of a recessive allele being scored as a dominant allele.  Default=0.005.}
  \item{estE2}{logical; should E2 be estimated?}
  \item{ped}{pedigree in 3 columns: id, dam, sire. Base individuals have NA as parents.}
  \item{estP}{logical; should the pedigree be estimated?}
  \item{beta}{vector of population-level parameters}
  \item{estbeta}{logical; should the population-level parameters be estimated?}
  \item{USdam}{vector of unsampled female population sizes}
  \item{estUSdam}{logical; should the female population sizes be estimated?}
  \item{USsire}{vector of unsampled male population sizes}
  \item{estUSsire}{logical or character; if \code{TRUE} the male population size is estimated separately from the female population size, if \code{"USdam"} male and female population sizes are constrained to be the same.}
  \item{shrink}{optional scalar for the variance defining the ridge-regression likelihood penalisation used to obatain starting values for beta and/or unsampled population sizes.}
}
\details{
If \code{estG=FALSE} an approximation is used for genotyping error.  In this case error rates and allele frequencies are not estimated but fixed at the starting parameterisation.  If indivdiuals have been typed more than once, then the approxiamtion only uses the genotype that first appears in the \code{GdP$G} object passed to \code{\link{MCMCped}}.  If \code{A} is not specified estimates are taken directly from \code{GdP$G} using \code{\link{extractA}}.  If \code{E1} and \code{E2} are not specified they are set to 0.005.  Note that if the approximation for genotyping error is used with codominant markers, Wang's (2005) model is not used, and the CEVUS model (Marshall 1998) is adopted.  In this case \code{E2} is the per-allele error rate and \code{E2}(2-\code{E2}) is the per-genotype error rate used by CERVUS.  If \code{dam} and \code{sire} are not specified the most likely set of parents given the genetic data are used (see \code{\link{MLE.ped}}).  The starting value of \code{beta}, if not given, is the MLE of beta given the starting pedigree (see \code{\link{MLE.beta}}).   The starting values of \code{USdam} and \code{USsire}, if not given, are the MLE based on the genotype data (see \code{\link{MLE.popsize}}).     
}
\value{list containing the arguments passed}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
# In this example we simulate a pedigree and then fix the
# pedigree and estimate the population level paarmeters

data(WarblerP)

var1<-expression(varPed(c("lat", "long"), gender="Male", 
  relational="OFFSPRING"))

# paternity is to be modelled as a function of distance 
# between offspring and male territories

res1<-expression(varPed("offspring", restrict=0))

# indivdiuals from the offspring generation are excluded as parents

res2<-expression(varPed("terr", gender="Female", relational="OFFSPRING",
  restrict="=="))

# mothers not from the offspring territory are excluded
 
PdP<-PdataPed(formula=list(var1,res1,res2), data=WarblerP, USsire=FALSE)
simped<-simpedigree(PdP, beta=-0.25)

# simulate a pedigree where paternity drops with distance (beta=-0.25)

sP<-startPed(ped=simped$ped, estP=FALSE)
model1<-MCMCped(PdP=PdP, sP=sP, nitt=3000, thin=2, burnin=1000) 
plot(model1$beta)

# The true underlying value is -0.25
}
}
\keyword{classes}
   
