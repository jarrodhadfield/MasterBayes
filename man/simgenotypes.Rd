\name{simgenotypes}
\alias{simgenotypes}
\title{Genotype and Genotyping Error Simulation}
\description{
Simulates genotypes given a pedigree and allele frequencies.  Option exists to simulate observed genotypes given Wangs's (2004) or CERVUS's model (Marshall 1998) of genotyping error for codominat markers or an asymmetric allele based model for dominant markers (Hadfield, 2009).
}
\usage{
simgenotypes(A, E1 = 0, E2 = 0, ped, no_dup = 1, prop.missing=0, marker.type="MSW")
}
\arguments{
  \item{A}{list of allele frequencies at each locus}
  \item{E1}{if Wang's (2004) model of genotyping error for co-dominant markers is used this is the probability of an allele dropping out. If CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers is used this parameter is not used. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a dominant allele being scored as a recessive allele.}
  \item{E2}{if Wang's (2004) or CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers are used this is the probability of an allele being miss-scored. In the CERVUS model errors are not independent for the two alleles within a genotype and so if a genotyping error has occurred at one allele then a genotyping error occurs at the other allele with probability one. Accordingly, \code{E2}(2-\code{E2}) is the per-genotype rate defined in CERVUS. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a recessive allele being scored as a dominant allele.}
  \item{ped}{pedigree in 3 columns: id, dam, sire. Base individuals have NA as parents. All parents must be in id.}
  \item{no_dup}{integer: number of times genotypes are to be observed}
  \item{prop.missing}{proportion of observed genotypes that are missing}
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
}
\value{
  \item{G}{list of genotype objects; true genotypes for each locus}
  \item{Gid}{vector of id names indexing \code{G}}
  \item{Gobs}{list of genotype objects; observed genotypes for each locus}
  \item{id}{vector of \code{id} names indexing \code{Gobs}}
}
\references{
Marshall, T. C. \emph{et al} (1998) Molecular Ecology 7 5 639-655
Kalinowski S.T. \emph{et al} (2007) Molecular Ecology 16 5 1099-1106
Hadfield J. D.  \emph{et al} (2009) \emph{in prep}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{genotype}}
\examples{
pedigree<-cbind(1:10, rep(NA,10), rep(NA, 10))

gen_data<-simgenotypes(A=list(loc_1=c(0.5, 0.2, 0.1, 0.075, 0.025)), 
 E1=0.1, E2=0.1, ped=pedigree, no_dup=1)

summary(gen_data$G[[1]])
summary(gen_data$Gobs[[1]])
}
\keyword{datagen}

