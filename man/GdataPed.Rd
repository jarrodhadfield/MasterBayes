\name{GdataPed}
\alias{GdataPed}
\alias{is.GdataPed}
\title{GdataPed Object}
\description{An object containing genotype data and the categories over which error rates may vary.}
\usage{
GdataPed(G, id = NULL, categories = NULL, perlocus=FALSE, marker.type="MSW")
}
\arguments{
  \item{G}{a list of \code{genotype} objects for each locus, or a \code{data.frame} to be coerced using \code{\link{genotype.list}}}
  \item{id}{a vector of individual identifiers associated with each genotype, individuals can have more than one observed genotype.  If \code{G} is a \code{data.frame} to be coerced and has a column name \code{id}, this will be used.}
 \item{categories}{an optional vector indicating subsets of genotypes that have different error rates.  If \code{G} is a \code{data.frame} to be coerced and has a column name \code{categories}, this will be used.}
 \item{perlocus}{if \code{TRUE} different error rates are estimated for each locus} 
 \item{marker.type}{\code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).}
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
GdP<-GdataPed(WarblerG)
}
}
\keyword{classes}

