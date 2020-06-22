\name{popsize.loglik}
\alias{popsize.loglik}
\title{Log-Likelihood of Unsampled Population Size}
\description{Log-likelihood of the number of unsampled individuals given the genotypes of offspring and potential parents
}
\usage{
popsize.loglik(X, USdam=FALSE, USsire=FALSE, nUS=NULL, ped=NULL, USsiredam=FALSE,
   shrink=NULL)
}
\arguments{
  \item{X}{list for each offspring with elements \code{N} and \code{G}.  \code{N} is a vector conatining the number of parental combinations in each of 4 classes.  \code{G} is a vector conatining the sum of the Mendelian transition probabilities over parental combinations in each class.  The 4 classes are parental combinations where a) both parents are sampled b) only sires are sampled, c) only dams are sampled d) neither parent is sampled.}
 \item{USdam}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled females exists. if    \code{USdam} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled females vary.}  
  \item{USsire}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled males exists. if    \code{USsire} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled males vary.}  
  \item{nUS}{vector for the size of the unsampled populations.  Parmeters for the unsampled female populations come before the male populations.}
  \item{ped}{optional pedigree with id, dam and sire in ech column}
  \item{USsiredam}{logical; if \code{TRUE} male and female unsampled populations sizes are constrained to be equal}
\item{shrink}{optional scalar for the variance defining the ridge-regression likelihood penalisation.}
}
\value{
  log-likelihood of the number of unsampled individuals given the genotype data.
}
\references{Nielsen. R. \emph{et.al} Genetics (2001) 157 4 1673-1682}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{Intended to be used within \code{\link{MLE.popsize}}}
\seealso{\code{\link{MCMCped}}, \code{\link{MLE.popsize}}}
\examples{
\dontrun{
data(WarblerG)
A<-extractA(WarblerG)

sex<-c(rep("Male", 50), rep("Female", 100))
offspring<-c(rep(0,100), rep(1, 50))
terr<-as.factor(rep(1:50, 3))
id<-1:150

res1<-expression(varPed(x="offspring", restrict=0)) 
res2<-expression(varPed(x="terr", gender="Female", relational="OFFSPRING", 
  restrict="=="))  
 
test.data<-data.frame(id, sex, offspring, terr)

PdP<-PdataPed(formula=list(res1, res2), data=test.data)

simped<-simpedigree(PdP)
G<-simgenotypes(A, E1=0, E2=0, ped=simped$ped, no_dup=1)

# remove 25 males at random, leaving 25

rm.males<-sample(1:50, 25, replace=FALSE) 

data.rm<-test.data[-rm.males,]
GdPrm<-GdataPed(G=lapply(G$Gobs, function(x){x[-rm.males]}), 
  id=G$id[-rm.males])

# delete genotype and phenotype records

PdPrm<-PdataPed(formula=list(res1, res2), data=data.rm, USsire=TRUE)

X.listrm<-getXlist(PdP=PdPrm, GdP=GdPrm, A=A, E2=0)

X<-lapply(X.listrm$X, function(x){list(N=c(25,0,1,0), 
  G=c(sum(x$G[1:25]), 0, x$G[26], 0))})

# each offspring has 1 mother and 25 sampled fathers so the 4 classes are:
# a) 1*25 categories with both parents sampled, 0*25 categries with only 
# sires sampled b) 1*1 categories with only dams sired and 0*0 categories
# with both sexes unsampled. 

nUS<-seq(10,40, length=100)
nUS_Loglik<-1:100
for(i in 1:100){
 nUS_Loglik[i]<-popsize.loglik(X, USsire=TRUE, nUS=nUS[i])
}
plot(nUS_Loglik~nUS, type="l", main="Profile Log-likelihood
   for number of unsampled males")
}
}
\keyword{models}
