\name{legalG}
\alias{legalG}
\title{Legal Genotype Configurations}
\description{
A function for checking whether a set of genotypes have a positive probability given the pedigree. If not, a legal configuration is found using heuristic methods.  Missing genotypes are also replaced with compatible genotypes.
}
\usage{
legalG(G, A, ped, time_born=NULL, marker.type="MSW")
}
\arguments{
  \item{G}{list of \code{genotype} objects}
  \item{A}{list of allele frequencies}
  \item{ped}{pedigree with id in the first column, dam in the second, and sire in the third. The genotypes must be in the same order as the id column}
  \item{time_born}{an optional vector for ordering a pedigree more efficiently (see \code{\link{orderPed}})} 
  \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
}
\value{
  \item{G}{a list of \code{genotype} objects with positive likelihood given the pedigree}
  \item{legal}{logical; TRUE if the the genotype configuration passed to \code{legalG} had a positive likelihood}
}
\references{
Marshall, T. C. \emph{et al} (1998) Molecular Ecology 7 5 639-655
Kalinowski S.T. \emph{et al} (2007) Molecular Ecology 16 5 1099-1106
Hadfield J. D.  \emph{et al} (2009) \emph{in prep}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
data(WarblerG)
A<-extractA(WarblerG[,16:17])

pedigree<-matrix(NA, 8,3)
pedigree[,1]<-1:8
pedigree[,2][5:8]<-c(1,1,2,2)
pedigree[,3][5:8]<-c(3,3,4,4)

G<-simgenotypes(A, E1=0, E2=0.3, ped=pedigree, no_dup=1)

newG<-legalG(G=G$Gobs,A=A,ped=pedigree)
newG$valid

# The input genotypes had a zero probability given the pedigree 
# (because of genotype error) but the output genotypes have 
# positive probability

legalG(newG$G,A,pedigree)$valid
}
}
\keyword{manip}
\keyword{models}
