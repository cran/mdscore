\name{lr.test}
\alias{lr.test}
\title{Likelihood ratio test for generalized linear models
}

\description{Computes the likelihood ratio test for the coefficients of a generalized linear model.
}

\usage{
lr.test(fit1, fit2)
}

\arguments{
  \item{fit1}{an object that stores the results of \code{glm} fit of the model under the null hypothesis.
}
  \item{fit2}{an object that stores the results of \code{glm} fit of the model under the alternative hypothesis.
}
}

\details{
The objects \code{fit1} and \code{fit2} are obtained using the usual options passed to the \code{glm} function.
}

\note{
Both \code{fit1} and \code{fit2} must have the same \code{family} and link function.
}

\value{
The function \code{lrt.test()} returns the following list of values:
\item{LR}{the value of the likelihood ratio statistic.}
\item{pvalue}{the p value of test under null hypothesis chi-square distribution.}
}

\references{
McCullagh P, Nelder J (1989). \emph{Generalized Linear Models}. Chapman & Hall/CRC, London.

Da Silva DN, Cordeiro GM (2009). "A Computer Program to Improve LR Tests for Generalized Linear Models." \emph{Communications in Statistics -- Simulation and Computation}, 38(10), 2184--2197.
}

\author{
Damiao N. da Silva \email{damiao@ccet.ufrn.br}

Antonio Hermes M. da Silva-Junior \email{hermes@ccet.ufrn.br}

}

\seealso{
\code{\link{mdscore}}

\code{\link{wald.test}}
}
\examples{
data(strength)
fitf <- glm(y ~ cut * lot, data = strength,family = inverse.gaussian("inverse"))
fit0 <- glm(y ~ cut + lot, data = strength, family = inverse.gaussian("inverse"))
lr.test(fit0,fitf)
}

\keyword{ likelihood }
\keyword{ ratio }
\keyword{ glm }
