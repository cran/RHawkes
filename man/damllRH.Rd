\name{damllRH}
\alias{damllRH}
\title{
Dynamically approxomated minus loglikelihood of a RHawkes model
}

\description{ Calculates an apprximation to the minus loglikelihood of a
RHawkes model with given immigration hazard function \eqn{\mu},
offspring birth time density function \eqn{h} and branching ratio
\eqn{\eta} relative to event times \code{tms} on interval \eqn{[0,cens]}.  }

\usage{
damllRH(tms, cens, par, q=0.999, qe=0.999,
      h.fn = function(x, p) dexp(x, rate = 1 / p), 
      mu.fn = function(x, p) {
        exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
        pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
      }, 
      H.fn = function(x, p) pexp(x, rate = 1 / p), 
      Mu.fn = function(x, p) {
        -pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
      },
      keepB=FALSE,
      H.inv=function(x,p)qexp(x,rate=1/p) )
}

\arguments{
  
  \item{tms}{
    A numeric vector, with values sorted in ascending order. Event times
    to fit the RHawkes point process model.
  }
  
  \item{cens}{
    A numericl scalar. The censoring time.
  }
  
  \item{par}{
    A numeric vector containing the parameters of the model,
    in order of the immigration parameters, in \eqn{\mu(.)}, offspring
    distribution parameters, in \eqn{h(.)}, and lastly the branching ratio
    \eqn{\eta(.)}.  }

  \item{q}{
    A numeric scalar in (0,1] and close to 1, which controls how far
    we look back when truncating the distribution of the most recent immigrant. 
  }

  \item{qe}{
    A numeric scalar in (0,1] and close to 1, which controls how to
    truncation is used in the offspring birth time distribution.
  }
  
  \item{h.fn}{ A (vectorized) function. The offspring birth time density
    function.
  }
  
  \item{mu.fn}{
    A (vectorized) function. The immigrant waiting time hazard function. 
  }
  
  \item{H.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the offspring birth time density function from 0 to \code{t}. 
  }
  
  \item{Mu.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the immigrant waiting time hazard function from 0 to \code{t}. 
  }

  \item{keepB}{ A boolean scalar, indicating whether the looking back
    values \code{B_i} should be part of the output or not.
  }

  \item{H.inv}{
    A (vectorized) function, giving the inverse function of the integral
    of the excitation.
  }
}

\value{ A scalar giving the value of the (approximate) negative
  log-likelihood, when \code{keepB} is FALSE (the default); A list with
  components \code{mll}, whhich contains the value of the negative
  log-likelihood, \code{Bs}, which gives the look-back order of the
  truncation of the distribution of the last immigrant, and \code{Bes},
  which gives the look-forward order in determining how far into the
  future the excitation effect is allowed to last.
}


\author{
  Feng Chen <feng.chen@unsw.edu.au>
}

\examples{
\dontrun{
## earthquake times over 96 years
data(quake);
tms <- sort(quake$time);
# add some random noise to the simultaneous occurring event times
tms[213:214] <- tms[213:214] + 
                    sort(c(runif(1, -1, 1)/(24*60), runif(1, -1, 1)/(24*60)))


## calculate the minus loglikelihood of an RHawkes with some parameters 
## the default hazard function and density functions are Weibull and 
## exponential respectively
mllRH(tms, cens = 96*365.25 , par = c(0.5, 20, 1000, 0.5))
damllRH(tms, cens = 96*365.25 , par = c(0.5, 20, 1000, 0.5),q=1,qe=1)
## calculate the MLE for the parameter assuming known parametric forms
## of the immigrant hazard function and offspring density functions.  
system.time(est <- optim(c(0.5, 20, 1000, 0.5), 
                        mllRH, tms = tms, cens = 96*365.25,
                        mu.fn=function(x,p)p[1]/p[2]*(x/p[2])^(p[1]-1),
                        Mu.fn=function(x,p)(x/p[2])^p[1],
                        control = list(maxit = 5000, trace = TRUE),
                        hessian = TRUE)
            )
system.time(est1 <- optim(c(0.5, 20, 1000, 0.5), 
                        function(p){
                            if(any(p<0)||p[4]<0||p[4]>=1)
                                return(Inf);
                            damllRH(tms = tms, cens = 96*365.25,
                                    mu.fn=function(x,p)p[1]/p[2]*(x/p[2])^(p[1]-1),
                                    Mu.fn=function(x,p)(x/p[2])^p[1],
                                    par=p,q=0.999999,qe=0.999999)
                        },
                        control = list(maxit = 5000, trace = TRUE),
                        hessian = TRUE)
            )
## point estimate by MLE
est$par
est1$par
## standard error estimates:
diag(solve(est$hessian))^0.5
diag(solve(est1$hessian))^0.5
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ RHawkes }
\keyword{ likelihood }% __ONLY ONE__ keyword per line
