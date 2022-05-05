simRHawkes <-
function(re.dist = rweibull, par.redist = list(shape = 1, scale = 1),
        ofspr.den = function(x,p.ofd) 1 / p.ofd * exp(-x/p.ofd),
        p.ofd = 1, branching.ratio = 0.5, cens = 1, B = 10, B0 = 50,
        max.ofspr.den = max(optimize(ofspr.den, c(0, cens), maximum = TRUE,
        p = p.ofd)$obj, ofspr.den(0, p.ofd), ofspr.den(cens, p.ofd)) * 1.1
        ){
  ## simulation of iid waiting times for immigrants
  wtms <- do.call(re.dist, args = c(n = B0, par.redist))
  n0 <- ceiling(((sqrt(4 * var(wtms) + 4 * mean(wtms) * cens) - 2 * sd(wtms)) /
                   (2 * mean(wtms))) ^ 2)
  wtms <- c(wtms, do.call(re.dist, args = c(n = n0, par.redist)))
  last.i.tm <- sum(wtms)
  while(last.i.tm <= cens){
    tmp <- do.call(re.dist, args = c(n = B, par.redist))
    wtms <- c(wtms,tmp)
    last.i.tm <- last.i.tm + sum(tmp)
  }
  tmp <- cumsum(wtms)
  i.tms <- tmp[tmp <= cens]
  ## for each immigrant simulate the corresponding offspring as an 
  ## inhomogenous poisson process using IHSEP
  of.tms <- sapply(i.tms,
              function(btm){
                tmp <- simHawkes1(nu = function(x) branching.ratio * ofspr.den(x, p.ofd),
                       nuM = max(max.ofspr.den * branching.ratio, 1e-5),
                       cens = cens - btm, 
                       g = function(x) branching.ratio * ofspr.den(x, p.ofd),
                       gM = max(max.ofspr.den * branching.ratio, 1e-5))
                if(length(tmp) == 0) return(numeric(0))
                sort(unlist(tmp)) + btm
                },
                simplify = FALSE)
  sort(c(i.tms, unlist(of.tms)))
}
