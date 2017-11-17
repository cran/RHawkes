\name{mllRH}
\alias{mllRH}
\title{
Minus loglikelihood of a RHawkes model
}
\description{
Calculates the minus loglikelihood of a RHawkes model with given
immigration hazard function \eqn{\mu}, offspring density function 
\eqn{h} and bracnhing ratio \eqn{\eta} for event times \code{tms} 
on interval \eqn{[0,cens]}.
}
\usage{
mllRH(tms, cens, par, 
      h.fn = function(x, p) dexp(x, rate = 1 / p), 
      mu.fn = function(x, p) {
        exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
        pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
      }, 
      H.fn = function(x, p) pexp(x, rate = 1 / p), 
      Mu.fn = function(x, p) {
        -pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
      })
}

\arguments{
  \item{tms}{
    A numeric vector, with values sorted in ascending order. Event times
    to fit the RHawkes point process model.
}
  \item{cens}{
    A scalar. The censoring time.
}
  \item{par}{
    A numeric vector containing the parameters of the model, in order of the 
    immigration parameters \eqn{\mu(.)}, offspring parameters \eqn{h(.)} and 
    lastly the branching ratio \eqn{\eta(.)}.
}
  \item{h.fn}{
    A (vectorized) function. The offspring density function. 
}
  \item{mu.fn}{
    A (vectorized) function. The immigration hzard function. 
}
  \item{H.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the offspring density function from 0 to \code{t}. 
}
  \item{Mu.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the immigrant hazard function from 0 to \code{t}. 
}
}

\value{
  The value of the negative log-liklihood. 
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>
}

\examples{
    ## simulate data from a RHawkes process (Hawkes, by the defaulting
    ## parameters)
    set.seed(1)
    ce <- 50;
    tms <- simRHawkes(cens=ce,par.redist=list(shape=1,scale=1.5),branching.ratio=0.3)
    ## calculate the MLE when data is treated as Hawkes process event
    ## times on [0,ce]
    sim.est0 <- optim(c(1.5,1,0.3),
                      function(p){
                          mloglik1d(jtms=tms,TT=ce,
                                    nu=function(x){rep(1/p[1],length(x))},
                                    gcoef=c(p[3]*p[2],p[2]),
                                    Inu=function(x){x/p[1]})
                      },
                      control=list(maxit=5000,trace=TRUE,reltol=1e-9),
                      )
    ## calculate the MLE when data is treated as RHawkes process event
    ## times on [0,ce]
    sim.est1 <- optim(c(1,sim.est0$par),
                      mllRH,tms=tms,cens=ce,
                      control=list(maxit = 5000, trace = TRUE,reltol=1e-9),
                      hessian=TRUE)
    sim.est1$par
\donttest{
## earthquake times over 96 years
data(quake);
tms <- sort(quake$time);
# add some random noise to the simultaneous occuring event times
tms[213:214] <- tms[213:214] + 
                    sort(c(runif(1, -1, 1)/(24*60), runif(1, -1, 1)/(24*60)))


## calculate the minus loglikelihood of an RHawkes with some parameters 
## the default hazard function and density functions are Weibull and 
## exponential respectivley
mllRH(tms, cens = 96*365.25 , par = c(0.5, 20, 1000, 0.5))
## calculate the MLE for the parameter assuming known parametric forms
## of the immigrant hazard function and offspring density functions.  
system.time(est <- optim(c(0.5, 20, 1000, 0.5), 
                        mllRH, tms = tms, cens = 96*365.25,
                        control = list(maxit = 5000, trace = TRUE),
                        hessian = TRUE)
            )
## point estimate by MLE
est$par
## standard error estimates:
diag(solve(est$hessian))^0.5
}
\dontrun{
    set.seed(1)
    ce <- 2000;
    tms <- simRHawkes(cens=ce,par.redist=list(shape=1,scale=2),branching.ratio=0.4)
    ## calculate the MLE when data is treated as Hawkes process event
    ## times on [0,ce]
    sim.est0 <- optim(c(2,1,0.4),
                      function(p){
                          mloglik1d(jtms=tms,TT=ce,
                                    nu=function(x){rep(1/p[1],length(x))},
                                    gcoef=c(p[3]*p[2],p[2]),
                                    Inu=function(x){x/p[1]})
                      },
                      control=list(maxit=5000,trace=TRUE,reltol=1e-9),
                      )
    ## calculate the MLE when data is treated as Weibull RHawkes
    ## process event times on [0,ce]
    system.time(sim.est1 <- optim(c(1,sim.est0$par),
                      mllRH,tms=tms,cens=ce,
                      control=list(maxit = 5000, trace = TRUE,reltol=1e-9),
                      hessian=TRUE))
    sim.est1$par
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RHawkes }
\keyword{ likelihood }% __ONLY ONE__ keyword per line
