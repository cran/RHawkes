EM1partial <- 
function(tms, cens, pars, maxiter = 1000, tol = 1e-8,
         h.fn = function(x, p) dexp(x, rate = 1 / p),
         mu.fn = function(x, p){
              exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
              pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
         },
         H.fn = function(x, p) pexp(x, rate = 1 / p),
         logg.fn = function(x, p){
              dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
              pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE) 
              - (x / p[2])^p[1]},
         Mu.fn = function(x, p){
              - pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
         }){
  p.mu <- pars[1:2]
  p.h <- pars[3]
  eta <- pars[4]
  n <- length(tms)
  iter<-1; diff<-10000
  
  while(diff > tol & iter <=maxiter) {
    ##auxiliary functions
    mu <- function(t) mu.fn(t, p.mu)
    Mu <- function(t) Mu.fn(t, p.mu)
    h <- function(t) h.fn(t, p.h)
    H <- function(t) H.fn(t, p.h)
    phi <- function(s) eta * sum(h(s - tms[tms < s]))
    mustar <- function(i,omega) sum(omega[i,1:(i - 1)] * 
                                      mu(tms[i] - tms[tms < tms[i]]))
    phi.fn <- function(s, p,eta) eta * sum(h.fn(s - tms[tms < s], p))
    
    #Now calculate the probabilities in the E-Step
    pik <- matrix(NA, nrow = n, ncol = n)
    for(i in 2:n) for(k in 1:(i - 1)) { 
      pik[i,k] <- mu(tms[i] - tms[k]) / (mu(tms[i] - tms[k]) + phi(tms[i]))
    }
    
    pijk <- array(NA, dim = c(n, n, n))
    for(i in 2:n) for(j in 1:(i - 1)) for(k in 1:(i - 1)){
      pijk <- eta * h(tms[i] - tms[j]) / (mu(tms[i] - tms[k]) + phi(tms[i]))
    }
    
    pi <- numeric(n)
    pi[1] <- 1
    pi[2] <- mu(tms[2] - tms[1]) / (mu(tms[2] - tms[1]) + phi(tms[2]))
    
    pij <- matrix(NA, nrow=n, ncol=n)
    pij[2,1] <- eta * h(tms[2]-tms[1]) / (mu(tms[2] - tms[1]) + phi(tms[2]))
    
    omega <- matrix(NA, nrow = n + 1, ncol = n + 1)
    omega[2,1] <- 1
    
    i <- 3
    while(i <= n){
      c <- 1 - pik[i-1,1:(i - 2)]
      omega[i,1:(i - 2)] <- omega[i-1,1:(i - 2)]*c
      omega[i,i - 1] <- pi[i - 1]
      pi[i] <- mustar(i, omega)/(mustar(i, omega) + phi(tms[i]))
      pij[i,1:(i - 1)] <- eta * h(tms[i] - tms[tms < tms[i]]) / 
        (mustar(i, omega) + phi(tms[i]))
      i <- i + 1
    }
    c <- 1 - pik[n, 1:(n - 1)]
    wnp1 <- c(omega[n, 1:(n - 1)] * c, pi[n])
    
    #Now the M-Step via maximizing the Q function into two independent problems. 
    Qmu <- function(p.mu){
      sum <- logg.fn(tms[1], p.mu)
      for(i in 2:n) for(j in 1:(i - 1)){
        sum <- sum + omega[i,j] * pik[i, j] * log(mu.fn(tms[i] - tms[j], p.mu)) - 
          omega[i, j] * pik[i, j] * Mu.fn(tms[i] - tms[j], p.mu)        
      }
      sum <- sum - sum(wnp1 * Mu.fn(cens - tms, p.mu))
      return(sum)
    }
    temp <- optim(par = p.mu, fn = Qmu, control = list(fnscale = -1))
    pars.mu <- temp$par
    
    Qh <- function(pars){
      tau0 <- pars[1]; eta <- pars[2]
      sum <- 0
      for(i in 2:n) for(j in 1:(i - 1)){
        sum <- sum + pij[i, j] * log(eta * h.fn(tms[i] - tms[j], tau0))
      }
      sum <- sum - eta * sum(H.fn(cens - tms,tau0))
      return(sum)
    }
    temp <- optim(par = c(p.h, eta), fn = Qh, control = list(fnscale = -1))
    pars.h <- temp$par[1]
    pars.eta <- temp$par[2]
    
    diff <- sum(abs(c(pars.mu, pars.h, pars.eta)-c(p.mu, p.h, eta))) 
    #Update Parameters
    p.mu <- pars.mu
    p.h <- pars.h
    eta <- pars.eta
    
    print(c(p.mu, p.h, eta))
    iter <- iter + 1
  }
  list(iterations = iter - 1, diff = diff, pars = c(p.mu, p.h, eta))
}