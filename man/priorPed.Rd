\name{priorPed}
\alias{priorPed}
\alias{is.priorPed}
\title{priorPed Object}
\description{An object containing the prior specifiactions for a model fitted using \code{\link{MCMCped}}.  If prior distributions are not specified then improper priors are used, and a proper posterior distribution cannot be gauranteed.}
\usage{
priorPed(E1=999, E2=999, beta=list(mu=999, sigma=999),
         USdam=list(mu=999, sigma=999),  
         USsire=list(mu=999, sigma=999))
}
\arguments{
  \item{E1}{matrix of parameters for the beta distribution specifying the prior distribution. If Wang's (2004) model of genotyping error for co-dominant markers is used this is the probability of an allele dropping out. If CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers is used this parameter is not used. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a dominant allele being scored as a recessive allele.  Rows correspond to error rate categories, columns to the beta shape parameters.  The order of rows in E1 are the order in which the error rate categories appear in the \code{categories} argument of \code{\link{GdataPed}} (see \code{dbeta}). If \code{perlocus=TRUE} was passed to \code{GdataPed}, then the error rate categories are replicated across loci}
  \item{E2}{matrix of parameters for the beta distribution specifying the prior distribution. If Wang's (2004) or CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers are used this is the probability of an allele being miss-scored. In the CERVUS model errors are not independent for the two alleles within a genotype and so if a genotyping error has occurred at one allele then a genotyping error occurs at the other allele with probability one. Accordingly, \code{E2}(2-\code{E2}) is the per-genotype rate defined in CERVUS. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a recessive allele being scored as a dominant allele.  Rows correspond to error rate categories, columns to the beta shape parameters.  The order of rows in E1 are the order in which the error rate categories appear in the \code{categories} argument of \code{\link{GdataPed}} (see \code{dbeta}). If \code{perlocus=TRUE} was passed to \code{GdataPed}, then the error rate categories are replicated across loci}
  \item{beta}{list containing a vector for the mean, and a matrix for the variance-covariances of a multivariate normal distribution, that specifies the prior distribution for the population level parameters.  The order of \code{beta} is the order in which the parameters appear in the MCMC ouput.}
  \item{USdam}{list containing vectors of means and standard deviations for log normal distributions that specify the prior distribution for the number of unsampled females.  The order of \code{USdam} is the order in which the unsampled dam populations appear in the \code{USdam} argument of \code{\link{PdataPed}} (see \code{dlnorm})}
  \item{USsire}{list containing vectors of means and standard deviations for log normal distributions that specify the prior distribution for the number of unsampled males. The order of \code{USsire} is the order in which the unsampled sire populations appear in the \code{USsire} argument of \code{\link{PdataPed}} (see \code{dlnorm})}
}
\value{list containing the arguments passed}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
# When each individual has only been genotyped once, and no pedigree 
# information exists, there is virtually no information available
# to estimate error rates.  The tiny amount of information comes 
# (dangerously) from the assumption of Hardy-Weinburg equilibrium.  
# The posterior distribution is similar to the prior:

data(WarblerG)
A<-extractA(WarblerG)

ped<-matrix(NA, 100,3)
ped[,1]<-1:100

G<-simgenotypes(A, E1=0.01, E2=0.01, ped=ped, no_dup=1)
GdP<-GdataPed(G=G$Gobs, id=G$id)
pP<-priorPed(E1=matrix(c(40,1600), nrow=1), E2=matrix(c(40,1600), nrow=1))

model1<-MCMCped(GdP=GdP, pP=pP)

#The posterior distribution recovers the prior distribution 

summary(model1$E1)
quantile(rbeta(1000, 40, 1600), prob=c(0.025, 0.25, 0.5, 0.75, 0.975))
}
}
\keyword{classes}
  
