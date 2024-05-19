sim_endog_linear_model <- function(beta, n, dist_vec, dist_params_list){
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
  
  # Draw (X, Z, e) where e is the final vector given by provided distributions
  drawn <- map2(dist_vec, dist_params_list, \(x, y) do.call(x, append(n, y))) %>% 
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>% 
    rename_with(\(col) 
                ifelse(
                  as.numeric(str_remove(col, "...")) < K ,
                  paste0("x", as.numeric(str_remove(col, "...")) + 1),
                  paste0("z", as.numeric(str_remove(col, "...")) + 2 - K)
                )) %>% 
    add_column("x1" = 1, .before = "x2") %>%
    add_column("z1" = 1, .before = "z2") %>%
    # name last column e, for structural error
    rename("e" = ncol(.))
  
  # Define table of the observed data (y,X, Z)
  observed <- drawn %>%
    mutate(y = as.numeric(as.matrix(drawn[,1:K]) %*% beta + e)) %>%
    select(-e)
  
  # Define design matrix
  X <- observed %>%
    select(-y, -starts_with("z")) %>%
    as.matrix()
  
  #Matrix of Instruments 
  Z <- drawn %>% 
    select(starts_with("z")) %>% 
    as.matrix()
  
  output <- list(
    "sim_draws" = drawn,
    "observed_data" = observed,
    "e" = drawn$e,
    "X" = X,
    "Z" = Z,
    "y" = observed$y,
    "model" = model
  )
  return(output)
}

