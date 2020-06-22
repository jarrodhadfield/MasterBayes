\name{modeP}
\alias{modeP}
\title{Posterior Mode of Parents}
\description{Finds the mode of the posterior marginal distribution of parents}
\usage{
modeP(postP, threshold=0, marginal=FALSE, USasNA=TRUE)
}
\arguments{
  \item{postP}{posterior distribution of parentage}
  \item{threshold}{threshold probability under which ML parents are replaced by NA}
  \item{marginal}{logical; should the marginal mode be calculated from the joint distribution?}
  \item{USasNA}{logical; should usampled parents be replaced by NA?}
}
\details{
Individuals that do not have a parent assignment with a posterior probability exceeding the \code{threshold}, or whose parents belong to the base or unsampled population (if \code{USasNA=TRUE}), have NA as their parents. Please bear in mind that the mode of the marginal distribution (returned by \code{\link{MCMCped}} if \code{write_postP="MARGINAL"}) may be different from the mode of the joint distribution (\code{write_postP="JOINT"}).  For example the male that has the highest marginal probability (marginal with respect to potential mothers) may not be the male that is in the parental category (i.e. dam/sire combination) with the highest probability.  If \code{write_postP="JOINT"} was sepcified, then the mode of the marginal distribution can be obtained by specifying \code{marginal=TRUE}.  The modes are marginal with respect to other offspring and with multigenerational pedigrees may not coincide with the mode of the distribution of pedigrees. 
}
\value{
   \item{P}{pedigree with id in the first column, and dam and sire in the second and third columns}
   \item{prob}{marginal posterior probability of the most likely parental combination (joint) or the most likely mother (marginal)}
   \item{prob.male}{marginal posterior probability of the most likely father (marginal)}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
\dontrun{
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

model1<-MCMCped(PdP=PdP, GdP=GdP, tP=tP, nitt=3000, thin=2, burnin=1000) 

ped<-modeP(model1$P, threshol=0.9)
ped
}
}
\keyword{distribution}
