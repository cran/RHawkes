simRHawkes1 <-
    function(re.dist = rweibull, par.redist = list(shape = 1, scale = 1),
             of.dis = "exp",par.ofdis = list(rate=1), branching.ratio = 0.5, 
             cens = 1, B = 10, B0 = 50, flatten=TRUE){
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
        of.tms <- lapply(i.tms,
                         function(btm){
                             btm + simoffspring(br = branching.ratio, dis = of.dis, par.dis = par.ofdis,
                                                 cens = cens - btm, sorted=FALSE)
                         })
        if(flatten)
            return( sort(c(i.tms, unlist(of.tms))) )
        list(immitimes=i.tms, offspringtimes=of.tms) 
}
