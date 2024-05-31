boot.continuous <- function(model, data, bootsize = 1000, bootseed = 1) {
  
  q <- list()
  
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}
  
  formula <- model$formula
  data <- data[, all.vars(formula), drop = FALSE]
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  
  # get only the X variables
  if(any(grepl("\\|", formula))) {
    lhs_formula <- strsplit(formula, "|", fixed = TRUE)[[1]][1]
    lhs_vars <- all.vars(as.formula(lhs_formula))
  } else
  {
    lhs_vars <- all.vars(formula[[3]])
  }  
  
  # Identify continuous or integer variables but not factors or binary
  continuous_or_integer_vars <- sapply(lhs_vars, function(var) {
    is.numeric(data[[var]]) && 
      !(is.factor(data[[var]]) || all(data[[var]] %in% c(0, 1)))
  })
  
  # Filter the variables
  numbers <- lhs_vars[continuous_or_integer_vars]
  
  set.seed(bootseed)
  boots <- matrix(, bootsize, length(numbers))
  
  if (model$dist == "Poisson") {
    for(i in 1:bootsize) {
      data[[all.vars(formula)[1]]] <- roipp(b, g, X, Z)
      modboot <- oneinfl(formula = formula, data = data, dist = "Poisson")
      for(j in 1:length(numbers)) {
        boots[i,j] <- hidden.continuous(model = modboot, data = data)[[j]]
      }
    }
  }
  
  q$se <- sqrt((colSums((boots - colMeans(boots)) ^ 2)) / (bootsize - 1))
  bootsrt <- apply(boots, 2, sort)
  q$lower95 <- bootsrt[bootsize*.025, ]
  q$upper95 <- bootsrt[bootsize*.975, ]
  q
}
