\name{strength}
\alias{strength}
\docType{data}
\title{Impact Strength an Insulating Material}

\description{
The dataset is a subsample of the 5 x 2 factorial experiment given by Ostle and Mensing (1963).}

\usage{data(strength)}
\format{
  A data frame with 30 observations on the following 3 variables.
  \describe{
    \item{\code{cut}}{type of specimen cut.}
    \item{\code{lot}}{lof of the material -- \code{I}, \code{II}, \code{III}, \code{IV} and \code{V}.}
    \item{\code{y}}{observations of the impact strength.}
  }
}

\source{
Ostle B, Mensing RW (1963). \emph{Statistics in Research: Basic Concepts and Techniques for Research Workers}. Iowa State University.
}

\examples{
data(strength)
fitf <- glm(y ~ cut * lot, data = strength,family = inverse.gaussian("inverse"))
summary(fitf)
X <- model.matrix(fitf, data = strength)
fit0 <- glm(y ~ cut + lot, data = strength, family = inverse.gaussian("inverse"))
test <- mdscore(fit0, X1=X[, 7:10])
summary(test)
}
\keyword{datasets}
