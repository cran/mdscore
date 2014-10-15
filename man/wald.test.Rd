\name{wald.test}
\alias{wald.test}
\title{Wald test for generalized linear models
}

\description{Computes the Wald score test for the coefficients of a generalized linear model.
}

\usage{
wald.test(model = model, terms)
}

\arguments{
  \item{model}{an object that stores the results of \code{glm} fit of the model under the null hypothesis.
}
  \item{terms}{number of coefficients to be tested under null hypothesis}
}

\details{
The object \code{model} is obtained using the usual options passed to the \code{glm} function.

}
\value{
The function \code{wald.test()} returns the following list of values:
\item{W}{the value of the Wald statistic.}
\item{pvalue}{the p value of test under null hypothesis chi-square distribution.}
}

\references{
McCullagh P, Nelder J (1989). \emph{Generalized Linear Models}. Chapman & Hall/CRC, London.
}

\author{
Damiao N. da Silva \email{damiao@ccet.ufrn.br}

Antonio Hermes M. da Silva-Junior \email{hermes@ccet.ufrn.br}
}

\seealso{
\code{\link{lr.test}}

\code{\link{mdscore}}
}
\examples{
data(strength)
fitf <- glm(y ~ cut * lot, data = strength,family = inverse.gaussian("inverse"))
wald.test(fitf,term=9)
}

\keyword{ Wald }
\keyword{ glm }
