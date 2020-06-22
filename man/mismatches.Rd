\name{mismatches}
\alias{mismatches}
\title{Parent-Offspring Genotype Mismatches}
\description{Calculates the number of mismatches between  parental and offspring genotypes, assuming the genotypes of spouses are unknown. Primarily intended to be used inside the function \code{\link{getXlist}} where potential parents can be excluded based on the number of mismatches. Dominant markers do not produce mismatches.
}
\usage{
mismatches(X.list, G, mm.tol=999)
}
\arguments{
   \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
   \item{G}{list of genotype objects, the rows of which must refer to the id vector \code{X.list$id}}
   \item{mm.tol}{maximum number of mismatches that are tolerated before exclusion}
}
\value{
    list of design matrices of the form \code{X.list}, but containing the number of mismatches between parents and offspring.  Potential parents that exceed the number of mismatches specified by \code{mm.tol} are removed from the vectors of potential parents: \code{restdam.id} and \code{restsire.id}.  
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{If a  \code{\link{GdataPed}} object is passed to \code{\link{getXlist}} then the number of mismatches will be calculated by default.
}
\seealso{\code{\link{MCMCped}}}
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

X.list<-getXlist(PdP)
# creates design matrices for offspring (in this case indivdiual "5")

X.list.MM<-mismatches(X.list, G=genotypes$Gobs, mm.tol=0)
# genetic likelihoods are arranged sires within dams 

X.list.MM$X$"5"$mmD
# number of mismatches between offspring "5" and dams "1" and "3"

X.list.MM$X$"5"$mmS
# number of mismatches between offspring "5" and sires "4" and "5"

X.list.MM$X$"5"$restdam.id
X.list.MM$X$"5"$dam.id
# dams with mismatches are excluded  mismatch (mm.tol=0)

X.list.MM$X$"5"$restsire.id
X.list.MM$X$"5"$sire.id
# sires with mismatches are excluded  mismatch (mm.tol=0)
}
}
\keyword{manip}
