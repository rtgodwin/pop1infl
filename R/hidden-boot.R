hidden.boot <- function(model, data, maxpred, bootsize = 1000, bootseed = 1) {
  
  q <- list()
  
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}
  
  formula <- model$formula
  data <- data[, all.vars(formula), drop = FALSE]
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  
  set.seed(bootseed)
  boots <- matrix(, bootsize, maxpred + 1)
  
  if (model$dist == "Poisson") {
    for(i in 1:bootsize) {
      data[[all.vars(formula)[1]]] <- roipp(b, g, X, Z)
      modboot <- oneinfl(formula = formula, data = data, dist = "Poisson")
      boots[i, ] <- hiddencounter(model = modboot, data = data, maxpred = maxpred)$counterfactual.values
    }
  }
  
  q$se <- sqrt((colSums((boots - colMeans(boots)) ^ 2)) / (bootsize - 1))
  bootsrt <- apply(boots, 2, sort)
  q$lower95 <- bootsrt[bootsize*.025, ]
  q$upper95 <- bootsrt[bootsize*.975, ]
  q
}
