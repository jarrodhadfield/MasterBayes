\name{getXlist}
\alias{getXlist}
\title{Design Matrices for the Multinomial Log-Linear Model}
\description{
  Forms design matrices for each offspring, and stores other relevant information. 
}
\usage{
getXlist(PdP, GdP=NULL, A=NULL, E1=0.005, E2=0.005, mm.tol=999)
}
\arguments{
    \item{PdP}{\code{\link{PdataPed}} object}
    \item{GdP}{optional \code{\link{GdataPed}} object}
    \item{A}{optional list of allele frequencies.  If not specified and \code{GdP} exists, allele frequencies are taken from \code{GdP$G} using \code{extractA}}
    \item{E1}{if Wang's (2004) model of genotyping error for co-dominant markers is used this is the probability of an allele dropping out. If CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers is used this parameter is not used. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a dominant allele being scored as a recessive allele.}
    \item{E2}{if Wang's (2004) or CERVUS's (Kalinowski, 2006; Marshall, 1998) model of genotyping error for co-dominant markers are used this is the probability of an allele being miss-scored. In the CERVUS model errors are not independent for the two alleles within a genotype and so if a genotyping error has occurred at one allele then a genotyping error occurs at the other allele with probability one. Accordingly, \code{E2}(2-\code{E2}) is the per-genotype rate defined in CERVUS. If Hadfield's (2009) model of genotyping error for dominant markers is used this is the probability of a recessive allele being scored as a dominant allele.}
    \item{mm.tol}{maximum number of genotype mismatches tolerated for potential parents}
}
\details{
  This is the main R routine for setting up design matrices for the various models that may be defined in the \code{formula} argument of \code{\link{PdataPed}}.  If a \code{\link{GdataPed}} object is passed to \code{getXlist} design matrices of genetic likelihoods are calculated (see \code{\link{fillX.G}}), and the number of mismatches between offspring and parental genotypes are stored (see \code{\link{mismatches}}).  \code{mm.tol} specifies the maximum number of mismatches that are tolerated between an offspring and a parent.  Parents that exceed this number of mismatches are excluded, and the design matrices for non-excluded parents are reordered by the number of mismatches. This increases the efficiency of sampling from the multinomial distribution of parents, because high probability parents appear first.
}
\value{
  \item{id}{vector of unique identifiers taken from \code{PdP}}
  \item{beta_map}{index relating the vector of unique parameters to the columns of the design matrices}
  \item{X}{list of design matrices and other information.}
}
\references{
Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31
Kalinowski S.T. \emph{et al} (2006) Molecular Ecology \emph{in press}
Hadfield J. D. \emph{et al} (2007) \emph{in prep}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\note{
Each element of \code{X} refers to an offspring (\code{names(X)}) and contains vectors for the set of potential parents (\code{restdam.id} and \code{restsire.id}) of each offspring.  Also included are the set of individuals that may have been parents but have been excluded for certain reasons (\code{dam.id} and \code{sire.id}).  Exclusion may have been based on the number of genotype mismatches, or it may have been on biological grounds (See the \code{keep} argument of \code{\link{varPed}}).  Parental id's are stored as integers which correspond to the actual id's stored in \code{id}.  Parental id's greater than the length of \code{id} refer to unsampled parents.

Six types of design matrix are used (\code{XDus}, \code{XDs}, \code{XSus}, \code{XSs}, \code{XDSus}, \code{XDSs}).  \code{XD..} are the design matrices for dams, and \code{XS..} are the design matrices for sires.  The rows of each design matrix are associated with individuals in dam.id and sire.id, respectively.  When interactions between dam and sire variables are modelled, or a \code{\link{varPed}} variable is created using the argument \code{relational="MATE"}, the design matrices vary over parental combinations.  \code{XDS..} are the design matrices for parental combinations with sire's varying the fastest.  Each of these three types of design matrix have two subclasses: \code{s} and \code{us}. \code{s} are design matrices which are fully observed, either because unsampled parents do not exist or because unsampled parents have known phenotypes (see argument \code{USvar} in \code{\link{varPed}}). \code{us} are for design matrices where the phenotypes of unsampled parents are unknown.  The matrices \code{XDus} and \code{Xsus} have a row of \code{NA}'s which correspond to the unsampled parent category.  The design matrix \code{XDSus} will typically have many rows of \code{NA}'s because each sampled parent may be paired to an unsampled individual.

When the argument \code{gender=NULL} is passed to \code{\link{varPed}} the respective columns in the dam and sire design matrices are associated with a single parameter.  Because of this the number of parameters to be estimated may be less than the total number of columns in the 6 design matrices.  \code{beta_map} relates a parameter vector to the columns of the design matrices.  The columns of the design matrices are numbered in the order they are introduced in the preceding paragraph (i.e XDus through to XDSs). The parameter vector is ordered identically except parameters associated with genderless variables are omitted for males. \code{par_order} is similar to \code{beta_map} but relates the order of the parameters specified in the \code{formula} argument to \code{\link{PdataPed}} to the respective columns of the design matrices.

If the argument \code{relational="OFFSPRING"} is specified in \code{\link{varPed}}, or the set of potential parents varies over offspring, the design matrices will vary across offspring.  For this reason I create a design matrix for each offspring irrespective of whether the matrices vary or not.  The design matrices for the genetic likelihoods will always vary over offspring.
}

\seealso{\code{\link{varPed}}, \code{\link{MCMCped}}}
\examples{
\dontrun{
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

X.list<-getXlist(PdP)
X.list$X$"19"$XSs

# For the first offspring we have the design matrix for sires
# The first column represents the distance between each male 
# and each offspring. The second column indicates the male's 
# mating type. Note that contrasts are set up with the first 
# male so the indicator variables may be negative.

matrix(X.list$X$"19"$XDSs, ncol=length(X.list$X$"19"$dam.id), 
   nrow=length(X.list$X$"19"$sire.id))

# incidence matrix indicating whether Females (columns) and Males (rows)
# are the same mating type. Again this is a contrast with the first 
# parental combination (which is +/+) so 0 actually represents parents
# with the same mating type.
}
}
\keyword{manip}
\keyword{models}
