hidden.continuous <- function(model, data) {
  
  q <- list()
  
  b <- model$beta
  if (model$dist == "negbin") {a <- model$alpha}
  
  names(b) <- substring(names(b), 2)
  
  formula <- model$formula
  
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  
  l <- exp(X %*% b)
  
  # replace the y values with their expected values from the counterfactual distribution
  if(model$dist == "Poisson") {
    evs <- l #* exp(l) / (exp(l) - 1)
  }
  data[[all.vars(formula)[1]]] <- evs
  
  # determine which variables are continuous or counts
  check.num.int <- function(X) {is.numeric(X) || is.integer(X)}
  numbers <- names(data)[sapply(data, check.num.int)]
  
  # Check if there is a pipe in the formula
  original_formula <- deparse(formula)
  if (grepl(" \\| ", original_formula)) {
    parts <- unlist(strsplit(original_formula, " \\| "))
    right_part <- unlist(strsplit(parts[2], " \\+ "))
  } else {
    parts <- original_formula
    right_part <- character(0)  # No variables on the right side
  }
  
  # Extract the dependent and independent variables from the left part
  left_part <- unlist(strsplit(parts[1], " ~ "))
  dependent_var <- left_part[1]
  independent_vars_left <- unlist(strsplit(left_part[2], " \\+ "))
  
  # Combine independent variables from both left and right parts (avoid duplications)
  all_independent_vars <- unique(c(independent_vars_left, right_part))
  
  # Generate new formulas for each variable in independent_vars_left,
  # but only if the variable is included in the 'numbers' vector
  new_formulas <- lapply(independent_vars_left, function(x) {
    if (x != dependent_var && x %in% numbers) {
      new_dependent_vars <- setdiff(all_independent_vars, x)
      formula(paste(x, "~", paste(c(dependent_var, new_dependent_vars), collapse = " + ")))
    }
  })
  
  # Filter out NULL entries if any
  new_formulas <- Filter(Negate(is.null), new_formulas)
  
  # weight the predicted values by the probability that each observation generates a 0
  weights <- (1 / (exp(l) - 1)) / sum(1 / (exp(l) - 1))
  
  # estimate the linear equations and get weighted predicted values
  for(i in 1:length(new_formulas)) {
    modlm <- lm(formula = new_formulas[[i]], data = data)
    q[[numbers[i]]] <-  sum(predict(modlm) * weights)
  }
  
  q
}