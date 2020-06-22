\name{MasterBayes}
\docType{package}
\alias{MasterBayes}
\title{Maximum Likelihood and Markov chain Monte Carlo methods for Pedigree Reconstruction, Analysis and Simulation.}
\description{The primary aim of MasterBayes is to use MCMC techniques to integrate over uncertainty in pedigree configurations estimated from molecular markers and phenotypic data.  Emphasis is put on the marginal distribution of parameters that relate the phenotypic data to the pedigree.  All simulation is done in compiled C++ using the Scythe Statistical Library. More detailed information can be found in \code{vignette("Tutorial", "MasterBayes")}.}
\details{

The motivation behind the package is to approximate the following probability distribution using Markov chain Monte Carlo techniques:

\deqn{p(\bf{\beta} | \bf{G}, \bf{y})}{p(beta | G, y)}
 
where \eqn{\bf{\beta}}{beta} is the vector of parameters of primary interest, \eqn{\bf{G}}{G} are the genetic data and  \eqn{\bf{y}}{y} are phenotypic data. Generally, it is not possible to simulate from the posterior distribution of \eqn{\bf{\beta}}{beta} when the problem is in this form and so I augment the parameter space with the pedigree, \eqn{\bf{P}}{P}:

\deqn{\int_{\bf{P}} p(\bf{\beta}, \bf{P} | \bf{G}, \bf{y})d\bf{P}}{int_P  p(beta, P | G, y)dP}

This simplifies the problem because the likelihood can be expressed more simply:

\deqn{L(\bf{G}, \bf{y} | \bf{\beta}, \bf{P}) = L(\bf{G} | \bf{P})L(\bf{y} | \bf{P}, \bf{\beta})}{L(G, y | beta, P) = L(G | P)L(y | P beta)}

This simplification rests on the assumption that the genetic and non-genetic data are independent after conditioning on the pedigree.  This will generally be true when markers are not linked to QTL's.  The first likelihood, \eqn{L(\bf{G} | \bf{P})}{L(G | P)}, is easily calculated for arbitrary pedigrees using the Elston-Stewart algorithm (Elston, 1971), and is based around the Mendelian transition probability. The second likelihood is obtained by fitting the multinomial log-linear model:

\deqn{L(\bf{y} | \bf{\beta}, \bf{P}) \propto p(\bf{P} | \bf{\beta}, \bf{y})p(\bf{P}).}{L(y | P, beta) = p(P |y, beta)p(P).}

Assuming that the set of possible pedigrees have equal prior probability, and that offspring are independently distributed after conditioning on the predictor variables: 

\deqn{L(\bf{y} | \bf{\beta}, \bf{P}) \propto  \prod^{n_{o}}_{i=1} \frac{e^{\bf{X}^{i}_{p_{i}}\bf{\beta}}}{\sum^{n_{p}}_{j=1}e^{\bf{X}^{i}_{j}\bf{\beta}}}.}{L(y | P, beta) = prod_i^no e^(X_i,p_i beta) / sum_j^np e^(X_i,j beta).}

where \eqn{\bf{X}^{i}_{j}}{X i,j} denotes the \eqn{j^{th}}{jth} row of offspring \eqn{i}{i}'s design matrix formed from the phenotypic data \eqn{\bf{y}}{y}. Each row of the design matrix corresponds to a parental combination. \eqn{n_{o}}{no} and \eqn{n_{p}}{np} denote the number of offspring and the number of potential parental combinations, respectively. \eqn{p_{i}}{p_i} denotes the actual parents of individual \eqn{i}{i} (Smouse, 1999).

This likelihood is evaluated over the probability distribution of the pedigree, \eqn{\bf P}{P}:

\deqn{p(\bf{P} | \bf{G}, \bf{y}, \bf{\beta}).}{p(P|G,y,beta).}

Most other techniques approximate this distribution as \eqn{p(\bf{P} | \bf{G})}{p(P|G)}, and even then tend to use the mode rather than the complete distribution, leading to inferential problems (See the information boxes in Hadfield et al. 2006).


Unfortunately, genotype data are rarely observed with out error and the parents of some offspring may not be sampled.  If the genetic markers are co-dominant then genotyping errors can be handled following either the model of Wang (2004) or CERVUS (Marshall 1998).  When the markers are dominant I model the probabilities of a dominant allele being miss-scored as a recessive and \emph{vice versa}. Denoting the parameters associated with these two forms of genotyping error as \eqn{\varepsilon_{1}}{E1} and \eqn{\varepsilon_{2}}{E2}, and the vector of parental allele frequencies as \eqn{\omega}{w}, two solutions are implemented.  

An exact solution:

\deqn{\int_{\bf{P}} \int_{\bf{G}} \int_{\varepsilon_{1}} \int_{\varepsilon_{2}} \int_{\omega} p(\bf{\beta}, \bf{P}, \bf{G}, \varepsilon_{1}, \varepsilon_{2}. \omega | \bf{G}^{(obs)}, \bf{y})d\bf{P}d\bf{G}d\varepsilon_{1}d\varepsilon_{2}d\omega}{int_P int_G int_E1 int_E2 int_w, p(beta, P, G, E1, E2, w, | G_obs, y)dPdGdE1dE2dw}

where the posterior probability distribution of the error rates, the allele frequencies and the true unobserved genotypes, \eqn{\bf{G}}{G}, are estimated and integrated out.  The conditional distribution of the true genotypes in the exact form is given by:

\deqn{p(\bf{G}^{obs} | \bf{G}, \varepsilon_{1}, \varepsilon_{2})p(\bf{G} | \bf{P}, \omega)}{p(G_obs | G,E1,E2)p(G | P, w).}

The second solution is an approximation to the above equation, and uses point estimates for \eqn{\omega}{w}, \eqn{\varepsilon_{1}}{E1} and \eqn{\varepsilon_{2}}{E2}. The conditional distribution of \eqn{\bf{G}}{G} is derived ignoring the information present in \eqn{\bf{P}}{P}:

\deqn{p(\bf{G}^{obs} | \bf{G}, \varepsilon_{1}, \varepsilon_{2})p(\bf{G} | \omega)}{p(G_obs | G,E1,E2)p(G | w)}

The approximation can be derived analytically, whereas the exact solution requires the Markov chain to be augmented with the true genotypes of all individuals.  This becomes very computer intensive but the approximation breaks down for dominant markers, or models in which the number of unsampled males and/or females is to be estimated. Unsampled parents are dealt with, and their number estimated using an approximation originally due to Nielsen (2001).  An exact solution to the problem has been proposed by Emery \emph{et.al.} (2001) but becomes impractical as the number of unsampled parents gets large. Nielsen's approximation is based around the Mendelian transition probability when a parental genotype is unknown.  This probability is derived using estimates of the allele frequencies at that locus and the assumption of Hardy-Weinberg equilibrium.  


I deal with the fact that unsampled individuals have missing phenotype data by approximating the distribution of the sum of linear predictors across unsampled parents.  This approximation relies on the assumption that the unsampled individuals come from the same statistical population as sampled individuals, and that population sizes are large enough so that the distribution for the sum tends to a normal distribution under the central limit theorem.

Taking \eqn{n}{n} and \eqn{N}{N} as the number of sampled individuals, and the total number of individuals in the population respectively:

\deqn{p(\sum^{N-n}\bf{\hat{p}}^{(miss)} | \bf{\hat{p}}^{(obs)}) \approx N(\frac{N-n}{n}\sum^{n}\bf{\hat{p}}^{(obs)}, \frac{N(N-n)}{n}S_{obs}^{2})}{p(sum(p_miss) | p_obs) = N(N-n mean(p_obs), (N*(N-n)*S^2)/n)}

where \eqn{\bf{\hat{p}}}{p} are vectors of linear predictors for the unsampled \eqn{^{(miss)}}{_miss} and sampled  \eqn{^{(obs)}}{_obs} individuals, respectively (Gelman \emph{et al.}, 2004). \eqn{S_{obs}^{2}}{S^2} is the sample variance of the observed linear predictors.

}
\references{
Elston, R. C. \& Stewart, J. Human Heredity (1971) 21 523-542
Emery, A. M. \emph{et.al} Molecular Ecology (2001) 10 1265-1278
Gelman, A. \emph{et.al} Bayesian Data Analysis \emph{Edition II} (2004) Chapman and Hall
Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31
Marshall, T. C. \emph{et al} (1998) Molecular Ecology 7 5 639-655
Nielsen. R. \emph{et.al} Genetics (2001) 157 4 1673-1682
Smouse P.E. \emph{et al} (1999) Journal of Evolutionary Biology 12 1069-1077
Wang J.L. Genetics (2004) 166 4 1963-1979
 }
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\keyword{package}
