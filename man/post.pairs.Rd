\name{post.pairs}
\alias{post.pairs}
\title{Returns pairs of individuals that fall into specific relatedness categories}
\description{
Computes posterior probabilities of pairs of indiviuals falling into specific relatedness categories (parent-offsping, sibs, full-sibs, half-sibs).  Returns those pairs that have a posterior probability greater than some threshold. 
}
\usage{
post.pairs(postP, threshold=0, rel="PO")
}
\arguments{
  \item{postP}{joint posterior distribution of parentage}
  \item{threshold}{threshold probability over which related pairs are returned}
  \item{rel}{relatedness category. Currently \code{"PO"} (Parent-Offspring), \code{"S"} (Sibs), \code{"FS"} (Full-Sibs) and \code{"HS"} (Half-Sibs) are supported.}
}
\value{
   \item{P}{pairs of indiviuals that fall into the \code{rel} category with posterior probability >  \code{threshold}}
   \item{prob}{posterior probability}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{modeP}}}
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

model1<-MCMCped(PdP=PdP, GdP=GdP, tP=tP, nitt=3000, thin=2, burnin=1000, write_postP="JOINT") 

fsib<-post.pairs(model1$P, threshol=0.9, rel="FS")
fsib$P
}
}
\keyword{manip}
\keyword{models}
