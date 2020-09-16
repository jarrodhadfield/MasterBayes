\name{PdataPed}
\alias{PdataPed}
\alias{is.PdataPed}
\title{PdataPed Object}
\description{
 \code{PdataPed} creates an object of class \code{PdataPed}, which typically contains the phenotype data to be passed to \code{\link{MCMCped}} and the formula that defines the model to be fitted.
 \code{is.PdataPed} returns \code{TRUE} if \code{x} is of class \code{PdataPed}
}
\usage{
PdataPed(formula, data=NULL, id=data$id, sex=data$sex,
   offspring=data$offspring, timevar=data$timevar, 
   USdam=FALSE, USsire=FALSE)
}
\arguments{
  \item{formula}{list of model predictors of the form \code{expression(\link{varPed}(...))}}
  \item{data}{data frame containing the predictor variables}
  \item{id}{vector of individual identifiers. If not specified, \code{data} must have an \code{id} column}
  \item{sex}{vector of individual sexes (either 'Male' or 'Female' or \code{NA}). If not specified individuals are assumed to be hermpahroditic unless \code{data} has a \code{sex} column}
  \item{offspring}{binary vector indicating whether records belong to offspring (1) or not (0)}
  \item{timevar}{an optional vector indicating cohorts for multigenerational pedigree reconstruction}
  \item{USdam}{logical or character; if \code{TRUE} a single undiferentaited population of unsampled females exists. If \code{USdam} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled females vary.} 
  \item{USsire}{logical or character; if \code{TRUE} a single undiferentaited population of unsampled males exists. If \code{USsire} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled males vary.}  
}
\value{list containing the arguments passed}
\details{
   If the number of unsampled individuals varies over subpopulations, and the parentage of an offspring is not restricted to ceratin subpopulations then the parameters will not be idenifiable.  This can be resolved by using an informative prior (see \code{\link{priorPed}}) for the number of unsampled individuals in each sub-population, or using the \code{restrict} argument in \code{\link{varPed}}.
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
id<-1:20
sex<-sample(c("Male", "Female"),20, replace=TRUE)
offspring<-c(rep(0,18),1,1)
lat<-rnorm(20)
long<-rnorm(20)
mating_type<-gl(2,10, label=c("+", "-"))

test.data<-data.frame(id, offspring, lat, long, mating_type, sex)

res1<-expression(varPed("offspring", restrict=0))
var1<-expression(varPed(c("lat", "long"), gender="Male", 
  relational="OFFSPRING"))
var2<-expression(varPed(c("mating_type"), gender="Female", 
  relational="MATE"))
var3<-expression(varPed("mating_type", gender="Male"))

PdP<-PdataPed(formula=list(res1, var1, var2, var3), data=test.data)

}
\keyword{classes}
 
