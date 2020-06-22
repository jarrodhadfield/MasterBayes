\name{modeG}
\alias{modeG}
\title{Posterior Mode of Genotypes}
\description{Finds the mode of the posterior marginal distribution of genotypes}
\usage{
modeG(postG, threshold=0)
}
\arguments{
  \item{postG}{posterior distribution of genotypes from an \code{\link{MCMCped}} model with argument \code{write_postG=TRUE}}
  \item{threshold}{threshold probability under which ML genotypes are replaced by NA}
}

\value{
   \item{G}{list of \code{genotype} objects}
   \item{id}{id vector}
}
\references{Hadfield J.D. \emph{et al}, Molecular Ecology  }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}, \code{genotype}}
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

model1<-MCMCped(PdP=PdP, GdP=GdP, tP=tP, nitt=3000, thin=2, burnin=1000, write_postG=TRUE) 

G<-modeG(model1$G)$G
summary(G[[1]])
}
}
\keyword{distribution}
