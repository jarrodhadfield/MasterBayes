\name{MCMCped}
\alias{MCMCped}
\title{Markov chain Monte Carlo Methods for Pedigree Reconstruction and Analysis}
\description{
  Markov chain Monte Carlo methods for estimating the joint posterior distribution of a pedigree and the parameters that predict its structure using genetic and non-genetic data.  These parameters can be associated with covariates of fecundity such as a sexually selected trait or age, or can be associated with spatial or heritable traits that relate parents to specific offspring.  Population size, allele frequencies, allelic dropout rates, and stochastic genotyping error rates can also be simultaneously estimated. 
}
\usage{
MCMCped(PdP=PdataPed(), GdP=GdataPed(), sP=startPed(), tP=tunePed(), 
   pP=priorPed(), mm.tol=999, nitt = 13000, thin = 10,  burnin = 
   3000, write_postG = FALSE, write_postA=FALSE, write_postP = 
   "MARGINAL", checkP = FALSE, jointP = TRUE, DSapprox=FALSE, verbose=TRUE)
}
\arguments{
  \item{PdP}{optional \code{\link{PdataPed}} object containing phenotypic data}
  \item{GdP}{optional \code{\link{GdataPed}} object containing genetic data}
  \item{sP}{optional \code{\link{startPed}} object containing starting parameterisation}
  \item{tP}{optional \code{\link{tunePed}} object containg tuning parameters for Metropolis Hastings updates}
  \item{pP}{optional \code{\link{priorPed}} object containg prior specifications}
  \item{mm.tol}{maximum number of mismatches tollerated}
  \item{nitt}{number of MCMC iterations}
  \item{thin}{thinning interval of the Markov chain}
  \item{burnin}{the number of initial iterations to be discarded}
  \item{write_postG}{if \code{TRUE} the marignal posterior distribution of true genotypes is stored}
  \item{write_postA}{if \code{TRUE} the joint posterior distribution of allele frequencies is stored}
  \item{write_postP}{if \code{"MARGINAL"} the marginal distribution of parents is stored. If \code{"JOINT"} the joint distribution of parents (the pedigree) is stored.}
  \item{checkP}{if \code{TRUE} the pedigree is checked for legality, and illegal pedigrees rejected. If \code{FALSE} it is assumed that any potential parent would produce a legal pedigree, i.e one without circuits, in the terminology of graph theory.}
  \item{jointP}{if \code{TRUE} both parents are sampled simultaneously, if \code{FALSE} each parent is sampled conditional on the other. \code{TRUE} should mix faster, but \code{FALSE} should iterate faster, especially when \code{relational="MATE"} is passed to \code{\link{varPed}}}
  \item{DSapprox}{if \code{TRUE} the likelihood for models in which a \code{relational="MATE"} variable is passed is approximated.  This can be much more efficient because the denominator of the multinomial is the summed linear pedictors for combinations in which i=m or j=m where m referes to the "MATE" at the current iteration.}
  \item{verbose}{if \code{TRUE} posterior samples and the Metropolis Hastings accpetance rates of \code{beta}, \code{USdam}, \code{USsire}, \code{E1}, \code{E2} are printed to the screen every 1000 iterations.}
}
\value{
  \item{beta}{an \code{mcmc} object containing samples from the posterior distribution of the population level parameters}
  \item{USdam}{an \code{mcmc} object containing samples from the posterior distribution of the number of unsampled females}
  \item{USsire}{an \code{mcmc} object containing samples from the posterior distribution of the number of unsampled males}
  \item{E1}{an \code{mcmc} object containing samples from the posterior distribution of allelic dropout rates for codominant markers or the probability of mis-scoring a dominant allele as recessive for dominant markers}
  \item{E2}{an \code{mcmc} object containing samples from the posterior distribution of stochasting genotyping error rates for codominant markers or the probability of mis-scoring a recessive allele as dominant for dominant markers}
  \item{G}{list of marginal distributions of true genotypes at each locus}
  \item{A}{list of \code{mcmc} objects containing samples from the posterior distribution of the base population allele frequencies at each locus}
  \item{P}{either samples from the posterior distribution of the pedigree, or the marginal distribution of parents}
}
\references{Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{getXlist}}}
\examples{

data(WarblerP)
data(WarblerG)

GdP<-GdataPed(WarblerG)

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
tP<-tunePed(beta=30)

model1<-MCMCped(PdP=PdP, GdP=GdP, tP=tP, nitt=300, thin=1, burnin=0) 

plot(model1$beta)

}
\keyword{models}
