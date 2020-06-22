\name{autocorrP}
\alias{autocorrP}
\title{Autocorrelation Function for Parenatge Assignment}
\description{Function for assessing mixing of the Markov chain with respect to parentage assignment.}
\usage{
autocorrP(postP)
}
\arguments{
  \item{postP}{JOINT posterior distribution of parentage}
}
\details{
For each offspring the proportion of transitions is calculated at lags 1, 2, 5, 10, 50 and 100 (i.e. the proportion of times that the parentage assignment at time t is different from the parentage assignment at time t+lag).  The difference between these proportions and the proportion at lag 1 is then calculated, and the mean over offspring given.  When the parentage assignments in successive MCMC iterations are independent these autocorrelation metrics should be randomly distributed about zero and should not decrease with increasing lag.
}
\value{
   matrix
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

# individuals from the offspring generation are excluded as parents

res2<-expression(varPed("terr", gender="Female", relational="OFFSPRING", 
  restrict="=="))

# mothers not from the offspring territory are excluded
 
PdP<-PdataPed(formula=list(var1,res1,res2), data=WarblerP, USsire=FALSE)
tP<-tunePed(beta=30)

model1<-MCMCped(PdP=PdP, GdP=GdP, tP=tP, nitt=3000, thin=1,  burnin=0, write_postP="JOINT") 
autocorrP(model1$P)
}
}
\keyword{distribution}
