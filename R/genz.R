#genz function
setClass("genz", representation(expected = "function",
                                integrand = "function",
                                name = "character"))
oscillatory <- new(Class = "genz",
                   name = "Oscillatory",
                   expected = function(alpha, beta) {
                     v1 <- cos(2 * pi * beta[1] + sum(alpha / 2))
                     v2 <- prod(sin(alpha / 2) / (alpha / 2))
                     v1 * v2
                   },
                   integrand = function(alpha, beta, z) {
                     cos(2.0 * pi * beta[1] + sum(alpha * z))
                   })
productpeak <- new(Class = "genz",
                   name = "Product Peak",
                   expected = function(alpha, beta) {
                     prod(alpha * atan((1.0 - beta) * alpha) + atan(beta * alpha))
                   },
                   integrand = function(alpha, beta, z) {
                     1.0 / prod(1.0 / alpha^2 + (z - beta)^2)
                   }
                   )
gaussian <- new(Class = "genz",
                name = "Gaussian",
                expected = function(alpha, beta) {
                  ab <- sqrt(2.0)
                  prod((sqrt(pi) / alpha) *
                         (pnorm((1.0 - beta) * ab * alpha)
                          - pnorm(-beta * ab * alpha)))
                },
                integrand = function(alpha, beta, z) {
                  exp(-min(sum((alpha * (z - beta))^2), 100.0))
                })

c0function <- new(Class = "genz",
                  name = "C0 Function",
                  expected = function(alpha, beta) {
                    ab <- alpha * beta
                    prod((2.0 - exp(-ab) - exp(ab - alpha)) / alpha)
                  },
                  integrand = function(alpha, beta, z) {
                    exp(-sum(alpha * abs(z - beta)))
                  }
                  )
discontinuous <- new(Class = "genz",
                     name = "Discontinuous",
                     expected = function(alpha, beta) {
                       prod((exp(alpha * beta) - 1.0) / alpha)
                     },
                     integrand = function(alpha, beta, z) {
                       if (any(beta < z)) {
                         0.0
                       } else {
                         exp(sum(alpha * z))
                       }
                     })

#cornerpeak.expected <- function(alpha, beta) {
#  
#}
#cornerpeak.integrand <- function(alpha, beta, z) {
#  b <- beta < 0.5
#  1 / (sum(b * z + !b * (alpha - z))^(length(z) + 1))
#}

# 森
# cornerpeak
cornerpeak.integrand <- function(alpha, beta, z) {
  (1 + sum(alpla * z))^(-length(z)+ 1)
}
cornerpeak.expected <- function(alpha, beta) {
  
}

