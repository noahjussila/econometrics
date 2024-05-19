sim_linear_model <- function(beta, n, dist_vec, dist_params_list){
  K <- length(beta)
  
  # store model details
  args <- dist_params_list %>% 
    map(\(x) unlist(append(n, x))) %>% 
    paste() %>% 
    str_remove(., "c")
  
  funcs <- paste(substitute(dist_vec))[-1]
  
  model <- list(
    "beta" = beta,
    "distribution" = paste0(funcs, args)
  )
  
  # Draw (X, e) where e is the final vector given by provided distributions
  drawn <- map2(dist_vec, dist_params_list, \(x, y) do.call(x, append(n, y))) %>%
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>% 
    # name first k - 1 columns x2,...,xk
    rename_with(\(col) paste0("x", as.numeric(str_remove(col, "...")) + 1), ) %>%
    # add x1 = 1 for intercept
    add_column("x1" = 1, .before = 1) %>%
    # name last column e, for structural error
    rename("e" = K + 1)
  
  # Define table of the observed data (y,X)
  observed <- drawn %>%
    mutate(y = as.numeric(as.matrix(across(1:K)) %*% beta + e)) %>%
    select(-e)
  
  # Define design matrix
  X <- observed %>%
    select(-y) %>%
    as.matrix()
  
  output <- list(
    "sim_draws" = drawn,
    "observed_data" = observed,
    "e" = drawn$e,
    "X" = X,
    "y" = observed$y,
    "model" = model
  )
}