\name{me.test}
\alias{me.test}

\title{A homogeneity Test under the Presence of Measurement Error}
\description{
	This function provides the test statistic and p-value of a homogeneity test of 
	distributions when the observations are measured with error.  
}
\usage{
me.test(W, V, B = 1000, wt = c("Uniform", "Normal"), wt.bd = NULL, 
	wt.prob = 0.99, nGL = 32)
}

\arguments{
  \item{W}{an \eqn{m_x} (>= 2) by \eqn{n_x} matrix of observations.}
  
  \item{V}{an \eqn{m_y} (>= 2) by \eqn{n_x} matrix of observations.}
  
  \item{B}{the number of bootstrap samples. Default is 1000.}
  
  \item{wt}{type of the weight function.}
  
  \item{wt.bd}{lower and upper bound of the weight function. 
	If \code{wt.bd} is not specified, bounds are computed based on
	the deconvoluted distribution function.}
  
  \item{wt.prob}{probability used to compute lower and upper bound.
	Will be ignored if \code{wt.bd} is provided.}
	
  \item{nGL}{the number of nodes for Gaussian quadrature}
}
\details{
	Based on our extensive simulations, we recommend to use \code{uniform}
	weight function with 0.99 probability. 
}
\value{
The output is an object of the class \code{htest} like in \code{t.test}.
	\item{statistic}{the value of the test statistic.}
	
	\item{p.value}{the p-value for the test.}
	
	\item{method}{the character string indicating the weight function.}
	
	\item{alternative}{a character string describing the alternative hypothesis.}
	
	\item{boundary}{lower and upper bound for the weight function.}
}
\references{
Lee, D., Lahiri, S. N. and Sinha, S. A Test of Homegeneity of Distributions when Observations
are Subject to Measurement Errors. \emph{Under revision}.
}
\author{
DongHyuk Lee, Samiran Sinha
}

\examples{
library(statmod)
set.seed(1234)
n <- 200
mx <- my <- 2
X <- rnorm(n, mean = 0, sd = 1)
Y <- rnorm(n, mean = 0.2, sd = 1)
Ux <- matrix(rnorm(n*mx, mean = 0, sd = 0.5), ncol = mx)
Uy <- matrix(rnorm(n*my, mean = 0, sd = 0.5), ncol = my)

W <- X + Ux
V <- Y + Uy
me.test(W, V, wt = "Uniform")
}
\keyword{homogeneity test}
\keyword{measurement error}
