\name{summary.mdscore}
\alias{summary.mdscore}

\title{
Summary methods for mdscore objects
}

\description{
  \code{summary} methods for the \code{mdscore} objects
}
  
\usage{
\method{summary}{mdscore}(object, \dots)
}

\arguments{
    \item{object}{object resulting from a run of the \code{\link{mdscore}} function.
        }
    \item{\dots}{not currently used}
}
\references{
{da Silva-Junior AHM, da Silva DN, Ferrari SLP  (2014). mdscore: An R Package to Compute Improved Score Tests in Generalized Linear Models. \emph{Journal of Statistical Software}, 61(2), 1-16., \url{http://www.jstatsoft.org/v61/c02/}}
}
\author{
Damiao N. da Silva \email{damiao@ccet.ufrn.br}
}

\seealso{
\code{\link{mdscore}}
}
\examples{
library(Sleuth3)
d <- transform(case1102, TLrat = Brain/Liver, Ltime = log(Time),
               Lwrat = log((Weight + Loss)/Weight),
               Treat = factor(Treatment == "BD", 
               labels=c("NS", "BD"))
               )
fitf <- glm(TLrat ~ Ltime * Treat + Days + Sex + Lwrat
            + Tumor + Treat*Lwrat, data = d,
            family = Gamma("log")
            )
X <- model.matrix(fitf)
fit0 <- glm(TLrat ~ Ltime * Treat + Sex + Lwrat + Tumor + Days,
            data=d, family=Gamma("log"))
test <- mdscore(fit0, X1=X[,9], phi=NULL)
summary(test)
}

\keyword{ mdscore }
\keyword{ score }
\keyword{ glm }
