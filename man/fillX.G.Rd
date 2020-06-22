\name{fillX.G}
\alias{fillX.G}
\title{Mendelian Transition Probabilities}
\description{
   This function is primarily intended for use within \code{\link{getXlist}}, and fills in the design matrices of the model with the genetic likelihoods.
}
\usage{
fillX.G(X.list, A, G, E1=0.005, E2=0.005, marker.type="MSW")
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
  \item{A}{list of allele frequencies}
  \item{G}{list of genotype objects; rows must correspond to individuals in the vector \code{X.list$id}}
  \item{E1}{if Wang's (2004) model of genotyping error for co-dominant markers is used this is the probability of an allele dropping out. If CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers is used this parameter is not used. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a dominant allele being scored as a recessive allele.}
  \item{E2}{if Wang's (2004) or CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers are used this is the probability of an allele being miss-scored. In the CERVUS model errors are not independent for the two alleles within a genotype and so if a genotyping error has occurred at one allele then a genotyping error occurs at the other allele with probability one. Accordingly, \code{E2}(2-\code{E2}) is the per-genotype rate defined in CERVUS. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a recessive allele being scored as a dominant allele.}
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
}
\value{
  list of design matrices of the form \code{X.list} containing genetic likelihoods for each offspring.
}
\references{
Marshall, T. C. \emph{et al} (1998) Molecular Ecology 7 5 639-655
Kalinowski S.T. \emph{et al} (2007) Molecular Ecology 16 5 1099-1106
Hadfield J. D.  \emph{et al} (2009) \emph{in prep}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{If a  \code{\link{GdataPed}} object is passed to \code{\link{getXlist}} then the genetic likelihoods will be calculated by default.
}
\seealso{\code{\link{getXlist}}}
\examples{
\dontrun{
data(WarblerG)
A<-extractA(WarblerG)

ped<-matrix(NA, 5,3)
ped[,1]<-1:5
ped[,2]<-c(rep(NA, 4), 1)
ped[,3]<-c(rep(NA, 4), 2)

genotypes<-simgenotypes(A, ped=ped)

sex<-c("Female", "Male", "Female", "Male","Female")
offspring<-c(0,0,0,0,1)

data<-data.frame(id=ped[,1], sex, offspring)

res1<-expression(varPed(x="offspring", restrict=0))

PdP<-PdataPed(formula=list(res1), data=data)
GdP<-GdataPed(G=genotypes$Gobs, id=genotypes$id)

X.list<-getXlist(PdP)
# creates design matrices for offspring (in this case indivdiual "5")

X.list.G<-fillX.G(X.list, A=A, G=genotypes$Gobs, E2=0.005)
# genetic likelihoods are arranged sires within dams 

X.list.G$X$"5"$dam.id
X.list.G$X$"5"$sire.id

# so for this example we have parental combinations 
# ("1","2"), ("1","4"), ("3","2"), ("2","4"):

X.list.G$X$"5"$G

# The true parents have the highest likelihood in this case
}
}
\keyword{manip}
