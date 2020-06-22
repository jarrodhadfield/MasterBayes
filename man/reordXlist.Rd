\name{reordXlist}
\alias{reordXlist}
\title{Reorders Design Matrices}
\description{Reorders design matrices so excluded parents appear last, and high probability parents appear first, thus increasing computational efficiency.
}
\usage{
  reordXlist(X.list, marker.type="MSW")
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}. Mismatch information must be present (see \code{\link{mismatches}})}
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
}
\details{
  The design matrices are reordered by the number of mismatches between a parent and offspring for codominant markers, and by the probability of the offspring genotype conditional on parent genotype for dominant markers.
}
\value{
  \code{X.list} for which parents are reordered 
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{If a  \code{\link{GdataPed}} object is passed to \code{\link{getXlist}} then the design matrices will be reordered by default.
}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
data(WarblerG)
A<-extractA(WarblerG)

ped<-matrix(NA, 5,3)
ped[,1]<-1:5
ped[,2]<-c(rep(NA, 4), 3)
ped[,3]<-c(rep(NA, 4), 4)

genotypes<-simgenotypes(A, ped=ped)

sex<-c("Female", "Male", "Female", "Male","Female")
offspring<-c(0,0,0,0,1)

data<-data.frame(id=ped[,1], sex, offspring)

var1<-expression(varPed(x="offspring", restrict=0))
PdP<-PdataPed(formula=list(var1), data=data)

X.list<-getXlist(PdP)
# creates design matrices for offspring (in this case indivdiual "5")

X.list<-mismatches(X.list, G=genotypes$Gobs)
X.list<-fillX.G(X.list, A=A, G=genotypes$Gobs)

X.list.reord<-reordXlist(X.list)

# The design matrices for the genetic likelihoods are reordered
# by the number of mismatches.  The true parental combination
# now appears first rather than last.

X.list$X$"5"$G
X.list.reord$X$"5"$G
}
}
\keyword{manip}
