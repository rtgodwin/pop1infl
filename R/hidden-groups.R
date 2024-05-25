hidden.groups <- function(model, data) {
  
  #q <- list()
  
  b <- model$beta
  if (model$dist == "negbin") {a <- model$alpha}
  
  names(b) <- substring(names(b), 2)
  
  formula <- model$formula
  
  # Remove pipe and any variables after pipe
  remove_after_pipe <- function(formula) {
    # Convert the formula to a string
    formula_string <- deparse(formula)
    # Split the string at the pipe symbol
    parts <- strsplit(formula_string, "|", fixed = TRUE)[[1]]
    # Select the first part, which is before the pipe
    before_pipe <- parts[1]
    # Trim any excess whitespace and convert back to a formula
    as.formula(trimws(before_pipe))
  }
  formula <- remove_after_pipe(formula)
  
  data <- data[, all.vars(formula), drop = FALSE]
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  #Z <- cleandata$Z
  
  # determine which variables are dummies
  is.dummy <- function(X) {length(unique(X)) == 2}
  dummies <- colnames(data[-1])[apply(data[-1], 2, is.dummy)]
  
  ### Poisson
  
  if (model$dist == "Poisson") {
    
    l <- exp(X %*% b)
    group.est <- matrix(, 1, 2*length(dummies))
    group.names <- levels(as.factor(data[[dummies]]))
    colnames(group.est) <- group.names
    
    for(i in 1:length(dummies)) {
      j <- 2*i - 1
      group.est[j] <- sum(1 / (exp(l[data[[dummies[i]]] == group.names[j]]) - 1))
      group.est[j + 1] <- sum(1 / (exp(l[data[[dummies[i]]] == group.names[j + 1]]) - 1))
    }
  }
  
  if (model$dist == "negbin") {
    q$dfee <- colMeans(dfee_nb(b, g, a, X, Z, dummies))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(dfee_nb(b, g, a, X, Z, dummies)), c("b", "g", "a")), "gradient")))
    q$sefee <- sqrt(diag(J %*% model$vc %*% t(J)))
  }
  group.est
}
