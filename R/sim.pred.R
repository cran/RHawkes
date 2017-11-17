sim.pred <-
function(tms, re.dist = rweibull, par,
         par.redist = list(shape = par[1], scale = par[2]),
         h.fn = function(x, p) dexp(x, rate = 1 / p),
         p.ofd = par[3], branching.ratio = par[4],
         cens, cens.tilde = cens * 1.5,
         mu.fn=function(x, p){
              exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
              pweibull(x, shape = p[1], scale = p[2],
              lower.tail = FALSE, log.p = TRUE))
         }){
  n <- length(tms)
  mu <- function(t) mu.fn(t, unlist(par.redist))
  h <- function(t) h.fn(t, p.ofd)
  phi <- function(s) branching.ratio * sum(h(s - tms[tms < s]))
  
  # Step 1: Simulate the last immigrant from I(T)
  fit.mod <- mllRH1(tms, cens, par)
  last.im.probs <- exp(fit.mod$log.p[n*(n - 1)/2 + 1:n])
  last.im <- sample(x = 1:n, 1, replace = T, prob = last.im.probs) 
  last.im.tm <- tms[last.im]
  
  # Step2: Simulate offspring of individuals already in the 
  # population at time cens
  RH0fit <-  simHawkes1(nu = function(x){ 
                sapply(x + cens, function(s) phi(s))},                  
                g = function(x) branching.ratio * h(x), 
                cens = cens.tilde - cens)
  tms.pred <- sort(cens + unlist(RH0fit))
  
  # Step 3: Simulate the next immigrant given its waiting time
  # is greater than cens - last.im.tm
  wt.fni <- rweibull(1, shape = par[1], scale = par[2])
  while(wt.fni <= cens - last.im.tm){
    wt.fni <- rweibull(1, shape = par[1], scale = par[2])
  }

  re.tm <- last.im.tm + wt.fni #renewal time
  tms.pred <- sort(c(tms.pred, re.tm))
  
  # Step 4/5: Simulate an RHawkes process from the new immigrant as the process
  # renews. Need to include offspring at the begining as the RHawkes process
  # assumes the first arrival is an immigrant. 
  if(wt.fni > cens.tilde - last.im.tm){
    return(tms.pred = tms.pred)
  }else{
    #RHawkes process starting from the first renewal
    RH0fit <- simRHawkes(par.redist = list(shape = par[1], scale = par[2]),
                p.ofd = par[3], branching.ratio = par[4],
                cens = cens.tilde - re.tm)
    tms.pred <- sort(c(tms.pred, re.tm + unlist(RH0fit)))
    #Offspring between new renewal and immigrant of RHawkes
    RH0fit <- simHawkes1(nu = function(x) branching.ratio * h(x),
                g = function(x) branching.ratio * h(x), 
                cens = cens.tilde - re.tm)
    tms.pred <- sort(c(tms.pred, re.tm + unlist(RH0fit))) 
    return(tms.pred = tms.pred)
  }
}