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
  
  # This code is specifically for the farm submissions data. Needs to be generalized to any model.
  size.mod <- lm(formula = new_formulas[[1]], data = data)
  dist.mod <- lm(formula = new_formulas[[2]], data = data)
  dairy.hat <- hidden.groups(model, data)[1,2] / (hidden.groups(model, data)[1,1] + hidden.groups(model, data)[1,2])
  q[[numbers[1]]] <- as.numeric(size.mod$coefficients[1] + size.mod$coefficients[3] * dist.mod$coefficients[1] + (size.mod$coefficients[3] * dist.mod$coefficients[4] + size.mod$coefficients[4]) * dairy.hat)
  q[[numbers[2]]] <- as.numeric(dist.mod$coefficients[1] + dist.mod$coefficients[3] * size.mod$coefficients[1] + (dist.mod$coefficients[3] * size.mod$coefficients[4] + dist.mod$coefficients[4]) * dairy.hat)
  
  # weight the predicted values by the probability that each observation generates a 0
  #weights <- (1 / (exp(l) - 1)) / sum(1 / (exp(l) - 1))
  
  # estimate the linear equations and get weighted predicted values
  #for(i in 1:length(new_formulas)) {
  #  modlm <- lm(formula = new_formulas[[i]], data = data)
  #  q[[numbers[i]]] <-  sum(predict(modlm) * weights)
  #}
  
  q
}
