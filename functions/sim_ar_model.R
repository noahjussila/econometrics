ar1 <- function(n, phi, t){
  u <- rnorm(n)
  accumulate(2:n, ~ phi * .x + u[.y], .init = u[1])
  output <- tibble(
    time = 1:n,
    value = accumulate(2:n, ~ phi * .x + u[.y], .init = u[1]),
    iter_num = t
  )
  return(output)
}

ar <- function(p, n, phi){
}

phi <- 9:1/10
p <- length(phi)
n <- 1e3
u <- rnorm(n + p)


