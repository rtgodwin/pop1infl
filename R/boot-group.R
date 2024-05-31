boot.group <- function(model, data, bootsize = 1000, bootseed = 1) {
  
  q <- list()
  
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}
  
  formula <- model$formula
  data <- data[, all.vars(formula), drop = FALSE]
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  
  is.dummy <- function(X) {length(unique(X)) == 2}
  dummies <- colnames(data[-1])[apply(data[-1], 2, is.dummy)]
  
  set.seed(bootseed)
  boots <- matrix(, bootsize, 2*length(dummies))
  
  if (model$dist == "Poisson") {
    for(i in 1:bootsize) {
      data[[all.vars(formula)[1]]] <- roipp(b, g, X, Z)
      modboot <- oneinfl(formula = formula, data = data, dist = "Poisson")
      boots[i, ] <- hidden.groups(model = modboot, data = data)
    }
  }
  
  q$se <- sqrt((colSums((boots - colMeans(boots)) ^ 2)) / (bootsize - 1))
  bootsrt <- apply(boots, 2, sort)
  q$lower95 <- bootsrt[bootsize*.025, ]
  q$upper95 <- bootsrt[bootsize*.975, ]
  q
}
