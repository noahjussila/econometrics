library(mvnfast)
library(tidyverse)

sim_gen_linear_model <- function(beta, n, mu_X, mu_e = 0, cov_Xe, cov_e = NULL, sked_fun = NULL){
  # number of regressors (excluding intercept)
  K <- length(beta) - 1
  # number of instruments (excluding intercept)
  L <- nrow(cov_Xe) - K - 1
  
  # Are any regressors endogenous? 
  endog <- (cov_Xe[-(K+L+1), K+L+1]!= 0) %>% 
    sum() %>% 
    as.logical()
  
  # Do errors exhibit any autocorrelation
  if(is.null(cov_e)){
    auto_corr <- FALSE
  } else {
    auto_corr <- (cov_e / diag(cov_e) != diag(n)) %>% 
      sum() %>% 
      as.logical()
  }
  
  # Do errors exhibit any heteroskedasticity
  if(missing(cov_e) & missing(sked_fun)){
    het <- FALSE
  } else if(missing(sked_fun)){
    het <- (diag(cov_e) != 1) %>% 
      sum() %>% 
      as.logical()
  } else{
    het <- TRUE
  }
  
  # If errors exhibit any heteroskedasticity, overwrite cov(e,e) in cov_Xe to be 1
  if(het){
    cov_Xe[K+L+1, K+L+1] <- 1
  }
  
  # Endogeneity + Autocorrelation -- draw from the n*(K+L-1) dimensional joint distribution for the entire sample 
  if(endog & auto_corr){
    Sigma <- kronecker(cov_Xe, diag(n))
    Sigma[(n*(K+L)+1):(n*(K+L+1)), (n*(K+L)+1):(n*(K+L+1))] <- cov_e
    mu <- kronecker(c(mu_X, mu_e), rep(1, n))
    
    drawn <- rmvn(1, mu, Sigma)[1,] %>% 
      as_tibble() %>% 
      mutate(
        rn = row_number() - 1,
        rn2 = floor(rn /n) + 1,
        var = case_when(
          rn2 <= K  ~ paste0("x", rn2 + 1),
          rn2 <= K + L ~ paste0("z", rn2 - K + 1),
          rn2 == K + L + 1 ~ paste0("e")
        ),
        i = (rn %% n) + 1
      ) %>% 
      select(
        value, 
        i, 
        var
      ) %>% 
      pivot_wider(names_from = var, values_from = value) %>% 
      select(-i) %>% 
      add_column("x1" = 1, .before = "x2") %>%
      add_column("z1" = 1, .before = "z2")
    
  # Exogeneity + Autocorrelation -- draw n observations of X from K-1 dimensional joint distribution and 1 observation of e from n dimensional joint distribution 
  } else if(auto_corr) {
    X <- rmvn(n, mu_X, cov_Xe[-1,-1])
    e <- rmvn(1, rep(mu_e, n), cov_e)[1,]
    Xe <- cbind(X, e)
    colnames(Xe)[1:(K-1)] <- paste0("x", 2:(K+1))
    drawn <- Xe %>% 
      as_tibble() %>% 
      add_column("x1" = 1, .before = "x2")
    
  # Endogeneity + No Autocorrelation -- draw n Xe observations from K+L-1 dimensional joint distribution
  } else if(endog){
    Xe <- rmvn(n, c(mu_X, mu_e), cov_Xe)
    colnames(Xe) <- c(
      paste0("x",2:(K+1)),
      paste0("z",2:(L+1)),
      "e"
    )
    drawn <- Xe %>% 
      as_tibble() %>% 
      add_column("x1" = 1, .before = "x2") %>%
      add_column("z1" = 1, .before = "z2")
    
  # Exogeneity + Autocorrelation -- draw n Xe observations from K-1 dimensional joint distribution
  }  else {
    Xe <- rmvn(n, c(mu_X, mu_e), cov_Xe)
    colnames(Xe) <- c(
      paste0("x", 2:(K+1)),
      "e"
    )
    drawn <- Xe  %>%
      as_tibble() %>%
      add_column("x1" = 1, .before = "x2")
  }
  
  # Now introduce heteroskedasticity
  sked_fun <- str_replace(sked_fun, "u", "e")
  if(het){
    if(!auto_corr & !is.null(cov_e)){
      drawn <- drawn %>% 
        mutate(e = e*sqrt(diag(cov_e)))
    }
    if(!is.null(sked_fun)){
      drawn <- drawn %>% 
        mutate(e = eval(parse(text = sked_fun)))
    }
  }
  
  # Define table of the observed data (y, X, Z)
  observed <- drawn %>%
    mutate(y = as.numeric(as.matrix(drawn[,1:(K+1)]) %*% beta + e)) %>%
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
    "y" = observed$y
  )
  return(output)
}
