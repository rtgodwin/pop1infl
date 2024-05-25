hidden.values <- function(model, data, maxpred) {
  
  b <- model$beta
  g <- model$gamma
  
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  y <- cleandata$y
  
  l <- exp(X %*% b)
  if(class(model) == "oneinflmodel"){
    Z <- cleandata$Z
    t <- exp(-Z %*% g)
  }
  
  if (model$dist == "negbin") {
    a <- model$alpha
    th <- l / a
  }
  
  if(missing(maxpred)) {
    maxpred = max(y)
  }
  
  z <- matrix(, (maxpred + 1), 3)
  
  if(model$dist == "Poisson") {
    if(class(model) == "oneinflmodel"){
      w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
      z[2,1] <- sum(w + (1 - w) * l / (exp(l) - 1))
      for(pp in 2:(maxpred)){
        z[pp+1,1] <- sum((1 - w) * (l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
      }
    }

    for(pp in 0:(maxpred)){
      z[pp+1,2] <- sum((l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
    }
    
  } else if(model$dist == "negbin") {

    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + t)
    z[2,1] <- sum(w + (1 - w) * a * ((1 / (1 + th)) ^ a) * (th / (1 + th - (1 + th) ^ (1 - a))))
    for(pp in 2:(maxpred)){
      z[pp+1,1] <- sum((1 - w) * (gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    } 
    #z[2:maxpred+1,1] <- pred
    # counterfactual
    #pred.c <- rep(0, maxpred + 1)
    for(pp in 0:(maxpred)){
      z[pp+1,2] <- sum((gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    #z[1:maxpred+1,2] <- pred.c
  }
  z[2:(maxpred+1), 3] <- tabulate(y)[1:maxpred]
  z <- as.data.frame(z)
  rownames(z) <- c(0:maxpred)
  colnames(z) <- c("predicted.values", "counterfactual.values","actual.values")
  if(class(model) == "truncmodel"){
    z <- z[c("actual.values", "counterfactual.values")]
  }
  if(class(model) == "oneinflmodel"){
    z <- z[c("actual.values", "predicted.values", "counterfactual.values")]
  }
  z
}
