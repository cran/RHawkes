mllRH <-
function(tms, cens, par,
         h.fn = function(x, p) dexp(x, rate = 1 / p),
         mu.fn = function(x, p){
              exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
              pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
        },
        H.fn = function(x, p) pexp(x, rate = 1 / p),
        Mu.fn = function(x, p){
          - pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
        }){
  p.mu <- par[1:2]
  p.h <- par[3]
  eta <- par[4]
  n <- length(tms)
  if(n == 0) return(Mu.fn(cens, p.mu))
  ## auxiliary functions
  mu <- function(t) mu.fn(t, p.mu)
  Mu <- function(t) Mu.fn(t, p.mu)
  h <- function(t) h.fn(t, p.h)
  H <- function(t) H.fn(t, p.h)
  phi <- function(s) eta * sum(h(s - tms[tms < s]))
  Phi <- function(t) eta * sum(H(t - tms[tms < t]))
  ## vectors to hold log(p(tau_i|tau_{1:i-1})), log(d_{i,1:(i-1)}),
  ## log(p_{i,1:(i-1)}), log(p_{i-1,1:(i-2)})
  lden <- lcon.den <- numeric(n);
  lpi <- lpimo <- numeric(n);
  res <- 0;
  lden[1] <- - Mu(tms[1]) + log(mu(tms[1]));
  res <- res - lden[1]
  i <- 2;
  while(i <= n + 1){
    if(i == 2){
      lpi[1] <- 0;
      if(i <= n){
        lcon.den[1] <- - Mu(tms[2] - tms[1]) - Phi(tms[2]) +
          log(mu(tms[2] - tms[1]) + phi(tms[2]))
      }else{
        lcon.den[1] <- - Mu(cens - tms[1]) - Phi(cens)
      }
    }else{
      ph <- phi(tms[i - 1])
      m <- mu(tms[i - 1] - tms[1:(i - 2)])
      lpi[1:(i - 2)] <-
        log(ph) - log(ph + m) + lcon.den[1:(i - 2)] + lpimo[1:(i - 2)]-lden[i - 1]
      mx <- max(lpimo[1:(i - 2)] + lcon.den[1:(i - 2)])
      lpi[i - 1] <- log(weighted.mean(m / (ph + m),
                                    exp(lcon.den[1:(i - 2)] +
                                        lpimo[1:(i - 2)] - mx))
                      )
      if(i <= n){
        lcon.den[1:(i - 1)] <- - Mu(tms[i] - tms[1:(i - 1)]) +
                                Mu(tms[i - 1] - tms[1:(i - 1)]) -
                                Phi(tms[i]) + Phi(tms[i - 1]) +
                                log(mu(tms[i] - tms[1:(i - 1)]) + phi(tms[i]))
      }else{
        lcon.den[1:(i - 1)] <- - Mu(cens - tms[1:(i - 1)]) +
                                Mu(tms[i - 1] - tms[1:(i - 1)]) -
                                Phi(cens) + Phi(tms[i - 1])
      }
    }
    mlcon.den <- max(lcon.den[1:(i - 1)])
    mlpi <- max(lpi[1:(i - 1)])
    lden[i] <- log(weighted.mean(exp(lcon.den[1:(i - 1)] - mlcon.den),
                             exp(lpi[1:(i - 1)] - mlpi))) + mlcon.den
    res <- res-lden[i]
    lpimo[1:(i - 1)] <- lpi[1:(i - 1)]
    i <- i + 1
  }
  res
}
