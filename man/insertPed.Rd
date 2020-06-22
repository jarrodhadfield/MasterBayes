\name{insertPed}
\alias{insertPed}
\title{Inserts Founders into a Pedigree}
\description{Inserts Founders into a Pedigree}
\usage{
insertPed(ped, founders=NULL)
}
\arguments{
  \item{ped}{pedigree with id, dam and sire in ech column}
  \item{founders}{optional vector of founder id's. If not specified, then parents without their own pedigree row are inserted}

}
\value{
  a pedigree pedigree with id, dam and sire in each column
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
pedigree<-matrix(NA, 7,3)
pedigree[,1]<-2:8
pedigree[,2][4:7]<-c(1,1,2,2)
pedigree[,3][4:7]<-c(3,3,4,4)

pedigree2<-insertPed(pedigree)
}
}
\keyword{manip}
