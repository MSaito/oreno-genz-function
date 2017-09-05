#genz function

# oscillatory region [0, 1)^s
oscillatory.expected <- function(alpha, beta) {
  v1 <- cos(2 * pi * beta[1] + sum(alpha / 2))
  v2 <- prod(sin(alpha / 2) / (alpha / 2))
  v1 * v2
}

# productpeak region [0, 1)^s
productpeak.expected <- function(alpha, beta) {
  prod(alpha * atan((1.0 - beta) * alpha) + atan(beta * alpha))
}

# gaussian region [0, 1)^s
gaussian.expected <- function(alpha, beta) {
  ab <- sqrt(2.0)
  prod((sqrt(pi) / alpha) *
         (pnorm((1.0 - beta) * ab * alpha)
        - pnorm(-beta * ab * alpha)))
}

# c0function
c0function.expected <- function(alpha, beta) {
  ab <- alpha * beta
  prod((2.0 - exp(-ab) - exp(ab - alpha)) / alpha)
}

#discontinuous
discontinuous.expected <- function(alpha, beta) {
  prod((exp(alpha * beta) - 1.0) / alpha)
}

oscillatory.integrand <- function(alpha, beta) {
  
}