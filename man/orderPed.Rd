\name{orderPed}
\alias{orderPed}
\title{Orders a Pedigree}
\description{Orders a pedigree so parents come before offspring}
\usage{
orderPed(ped, time_born=NULL)
}
\arguments{
  \item{ped}{pedigree with id, dam and sire in ech column}
  \item{time_born}{an optional vector of birth dates by which the pedigree can be ordered)} 
}
\value{
  an ordered pedigree pedigree with id, dam and sire in each column
}
\note{This function has changed name from \code{order.ped} in earler versions <2.42. \code{order.ped} did not always (rarely) ordered the pedigree correctly.  This new function uses the \code{kindepth} function from the \code{kinship2} package}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
pedigree<-matrix(NA, 8,3)
pedigree[,1]<-1:8
pedigree[,2][5:8]<-c(1,1,2,2)
pedigree[,3][5:8]<-c(3,3,4,4)

pedigree<-pedigree[sample(1:8),]

pedigree2<-orderPed(pedigree)
}
}
\keyword{manip}
