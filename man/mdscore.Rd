\name{mdscore}
\alias{mdscore}
\title{Modified score test for generalized linear models
}

\description{Computes the modified score test based for the coefficients of a generalized linear model.
}

\usage{
mdscore(model = model, X1 = X1, phi = NULL)
}

\arguments{
  \item{model}{an object that stores the results of \code{glm} fit of the model under the null hypothesis.
}
  \item{X1}{the matrix with the columns of the model matrix X that correspond to the coefficients being specified in the null hypothesis.
}
  \item{phi}{the precision parameter.
}
}

\details{
The object \code{fit.model} is obtained using the usual options passed to the \code{glm} function.

}
\value{
The function \code{mdscore()} returns the following list of values:
\item{Sr}{the value of the score statistic.}
\item{Srcor}{the value of the modified score statistic.}
\item{coef}{a vector with the coefficients A1 , A2 and A3.}
\item{n}{the total sample size.}
\item{df}{the number of degrees of freedom of the chi--squared approximations for the tests.}
\item{phi}{the precision parameter used in the computations}
}

\references{
{Cordeiro GM, Ferrari SLP (1991). A Modified Score Test Statistic Having chi-squared Distribution to Order n--1 . \emph{Biometrika}, 78(3), 573--582.}

{Cordeiro GM, Ferrari SLP, Paula GA (1993). Improved Score Tests for Generalized Linear Models. \emph{Journal of the Royal Statistical Society B}, 55(3), 661--674.}

{Cribari-Neto F, Ferrari SLP (1995). Second Order Asymptotics for Score Tests in Generalised Linear Models. \emph{Biometrika}, 82(2), 426--432.}

{da Silva-Junior AHM, da Silva DN, Ferrari SLP  (2014). mdscore: An R Package to Compute Improved Score Tests in Generalized Linear Models. \emph{Journal of Statistical Software}, 61(2), 1-16., \url{http://www.jstatsoft.org/v61/c02/}}
}
\author{
Antonio Hermes M. da Silva-Junior \email{hermes@ccet.ufrn.br}

Damiao N. da Silva \email{damiao@ccet.ufrn.br}
}

\seealso{
\code{\link{summary.mdscore}}
}
\examples{
data(strength)
fitf <- glm(y ~ cut * lot, data = strength,family = inverse.gaussian("inverse"))
summary(fitf)
X <- model.matrix(fitf, data = strength)
fit0 <- glm(y ~ cut + lot, data = strength, family = inverse.gaussian("inverse"))
mdscore(fit0, X1=X[, 7:10])
}

\keyword{ mdscore }
\keyword{ score }
\keyword{ glm }
