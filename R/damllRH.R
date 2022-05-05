damllRH <- function(tms,cens,par,q=0.999,qe=0.999,
                    h.fn=function(x,p)dexp(x,rate=1/p),
                    mu.fn=function(x,p){
                        exp(dweibull(x,shape=p[1],scale=p[2],log=TRUE)-
                            pweibull(x,shape=p[1],scale=p[2],lower.tail=FALSE,
                                     log.p=TRUE))
                    },
                    H.fn=function(x,p)pexp(x,rate=1/p),
                    Mu.fn=function(x,p){
                        -pweibull(x,shape=p[1],scale=p[2],lower.tail=FALSE,log.p=TRUE)
                    },
                    keepB=FALSE,
                    H.inv=function(x,p)qexp(x,rate=1/p)){
    p.mu <- par[1:2]
    p.h <- par[3]
    eta <- par[4]
    n <- length(tms)
    if(n == 0)
        return(Mu.fn(cens, p.mu))
    mu <- function(t) mu.fn(t, p.mu)
    Mu <- function(t) Mu.fn(t, p.mu)

    lden <- lcon.den <- numeric(n)
    lpi <- lpimo <- numeric(n)
    res <- 0
    lden[1] <- -Mu(tms[1]) + log(mu(tms[1]))
    res <- res - lden[1]

    ph <- numeric(n);
    Q <- H.inv(qe,p.h);
    if(keepB){
        Bes <- numeric(n);
        i <- 2
        while(i<=n){
            Be <- n; if(tms[n]-tms[i-1] >Q)Be <- i-1+ min(which(tms[i:n]-tms[i-1] >Q));
            ph[i:Be] <- ph[i:Be]+eta*h.fn(tms[i:Be]-tms[i-1],p.h);
            Bes[i] <- Be-i+1
            i <- i+1
        }
    }else{
        i <- 2
        while(i<=n){
            Be <- n; if(tms[n]-tms[i-1] >Q)Be <- i-1+ min(which(tms[i:n]-tms[i-1] >Q));
            ph[i:Be] <- ph[i:Be]+eta*h.fn(tms[i:Be]-tms[i-1],p.h);
            i <- i+1
        }
    }
    Ph <- eta*sum(H.fn(cens-tms,p.h))

    lpi[1] <- 0
    if(2 <= n){
        lcon.den[1] <- -Mu(tms[2] - tms[1]) +
            log(mu(tms[2] - tms[1]) + ph[2])
    } else{
        lcon.den[1] <- -Mu(cens - tms[1])
    }
    lden[2] <- lcon.den[1];
    res <- res - lden[2]
    lpimo[1] <- lpi[1];

    Bs <- integer(n);
    Bs[2-1] <- 1;
    ## lq <- log(q);
    i <- 3;
    while(i <= n+1){
        B <- Bs[i-1-1];
        ## ph <- phi(tms[i - 1]);
        m <- mu(tms[i - 1] - tms[(i - B - 1):(i - 2)]);
        lpi[(i - B - 1):(i -2)] <- log(ph[i-1]) - log(ph[i-1] + m) + lcon.den[(i - B - 1):(i - 2)] +
            lpimo[(i - B - 1):(i - 2)] - lden[i - 1]
        mx <- max(lpimo[(i - B - 1):(i -2)] + lcon.den[(i - B - 1):(i - 2)])
        lpi[i - 1] <- log(weighted.mean(m/(ph[i-1] + m),
                                        exp(lcon.den[(i - B - 1):(i - 2)] + lpimo[(i - B - 1):(i - 2)] - mx)))
        ## j <- 1; lcp <- lpi[i-j];
        ## while(lcp < lq && j<=B+1){
        ##     j <- j+1;
        ##     lcp <- lcp+log(1+exp(lpi[i-j]-lcp))
        ## }
        j <- 1; cp <- exp(lpi[i-j]);
        while(cp < q && j<=B+1){
            j <- j+1;
            cp <- cp+exp(lpi[i-j])
        }
        if(j<=B+1){
            B <- j
        }else{
            B <- j-1;
            cp <- 1;
        }
        Bs[i-1] <- B;
        if(i <= n){
            lcon.den[(i - B ):(i - 1)] <- - Mu(tms[i] - tms[(i - B ):(i - 1)]) +
                Mu(tms[i - 1] - tms[(i - B ):(i - 1)]) +
                log(mu(tms[i] - tms[(i - B ):(i - 1)]) + ph[i])
        } else {
            lcon.den[(i - B ):(i - 1)] <- - Mu(cens - tms[(i - B ):(i - 1)]) +
                Mu(tms[i - 1] - tms[(i - B ):(i - 1)])
        }
        mlcon.den <- max(lcon.den[(i - B):(i - 1)])
        mlpi <- max(lpi[(i - B ):(i - 1)])
        lden[i] <- log(weighted.mean(exp(lcon.den[(i - B ):(i - 1)] - mlcon.den),
                                     exp(lpi[(i - B ):(i - 1)] - mlpi))) + mlcon.den
        lpimo[(i-B):(i - 1)] <- lpi[(i-B):(i - 1)]-log(cp);#lq;

        res <- res - lden[i]
        i <- i + 1
    }
    res <- res + Ph ## put Ph in
    if(!keepB) return(res);
    list(mll=res,Bs=Bs,Bes=Bes)
}
