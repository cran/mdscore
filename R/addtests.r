# file mdscore/R/addtests.R
# by  Damiao N. da Silva and Antonio Hermes M. da Silva Jr.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

wald.default <- function(model = model, terms){
  
  X   <- model.matrix(model)
  X1  <- as.matrix(X[,terms])
  X2  <- as.matrix(X[,-terms])
  w   <- as.numeric(model$weights)
  n   <- nrow(X)
  pp  <- switch(model$family$family,
                gaussian         = as.numeric(1/summary(model)$dispersion)*n/summary(model)$df.residual,
                Gamma            = as.numeric(gamma.shape(model)$alpha),
                inverse.gaussian = as.numeric(1/summary(model)$dispersion)*n/summary(model)$df.residual,
                poisson          = 1,
                binomial         = 1)
  
  IX2WX2 <- invmat(crossprod(X2, w*X2))  
  X2WX1  <- crossprod(X2, w*X1)  
  C      <- crossprod(IX2WX2, X2WX1) 
  R      <- X1 - crossprod(t(X2), C) 
  RWR    <- crossprod(R, w*R) 
  b1     <- coef(model)[terms]
  
  est    <- as.numeric(pp*t(b1)%*%RWR%*%b1)
  pvalue <- pchisq(est, df=ncol(X1), lower.tail=FALSE)
  return(list(est=est, pvalue=pvalue))
}

wald.test <- function(model = model, terms){
  if(class(model)[1]!="glm")
    stop(paste("\n 'model' is not an object from 'glm' class", "!!!\n"))
  
  fam <- summary(model)$family$family 
  if(!(fam %in% c("binomial","gaussian", "Gamma", "inverse.gaussian","poisson"))) 
    stop(paste("\n When the parameter phi is unknown, the model family must be 'gaussian', 'Gamma' or 'inverse.gaussian'","!!!\n"))
  
  rslt = wald.default(model = model, terms)
  
  structure(
    list(
      W = rslt$est,
      pvalue = rslt$pvalue
      ),
    class="wald.test"
    )   
}

lr.test <- function(fit1, fit2){
  if(class(fit1)[1]!="glm")
    stop(paste("\n 'fit1' is not an object from 'glm' class", "!!!\n"))
  
  if(class(fit2)[1]!="glm")
    stop(paste("\n 'fit2' is not an object from 'glm' class", "!!!\n"))
  
  if(!(all(fit1$y %in% fit2$y)))
    stop(paste("\n 'fit1' and 'fit2' must have the same response variable","!!!\n"))
  
  if(!(all(fit1$family %in% fit2$family)))
    stop(paste("\n 'fit1' and 'fit2' must have the same family and link function","!!!\n"))
  
  y         <- fit1$y
  statistic <- 2*(logLik(fit2)[1] - logLik(fit1)[1])
  g1 <- fit1$df.residual
  g2 <- fit2$df.residual
  df <- g1 - g2
  if(statistic < 0){
    statistic <- abs(statistic)
    df <- abs(df)
    warning(paste("\n 'fit1' must be under null hypothesis","!!!\n"))
  }
  pvalue    <- pchisq(statistic,df,lower.tail=FALSE)
  
  structure(
    list(
      "LR"         = statistic,
      "pvalue"     = pvalue
    ),
    class="lrt.test"
  )   
}