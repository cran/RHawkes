pred.haz <-
function(x, tms, cens, par,
        h.fn = function(x, p) dexp(x, rate = 1 / p),
        mu.fn=function(x, p){
            exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
            pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
        },
        H.fn = function(x, p) pexp(x, rate = 1 / p),
        Mu.fn=function(x, p){
            - pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
        }){
  p.mu <- par[1:2]
  p.h <- par[3]
  eta <- par[4]
  n <- length(tms)
  ## auxillary functions
  mu <- function(t) mu.fn(t, p.mu)
  Mu <- function(t) Mu.fn(t, p.mu)
  h <- function(t) h.fn(t, p.h)
  H <- function(t) H.fn(t, p.h)
  phi <- function(s) eta * sum(h(s - tms[tms < s]))
  Phi <- function(t) eta * sum(H(t - tms[tms < t]))
  
  Stilde <- exp( - Mu(outer(cens + x, tms, "-")) +
                Mu(outer(rep(tms[n], length(x)), tms, "-")) -
                (sapply(cens + x, Phi) - Phi(tms[n]))
                )
  tmp <- mllRH1(tms, cens, par)
  pnp1 <- exp(tmp$log.p[n*(n-1)/2 + 1:n])
  
  dnp1 <- Stilde*(
            mu(outer(cens + x, tms, "-")) +
            sapply(cens + x, phi)
            )
  rowSums(dnp1 * outer(rep(1, length(x)), pnp1))/
    rowSums(Stilde*outer(rep(1, length(x)), pnp1))
}