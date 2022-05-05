\name{EM2partial}
\alias{EM2partial}
\title{
Partial EM algorithm for the RHawkes process, version 2
}
\description{
Calculates the RHawkes model parameters via a partial Expectation-Maximization 
(EM2) algorithm of Wheatley, Filimonov and Sornette (2016).
}
\usage{    
EM2partial(tms, cens, pars, maxiter = 1000, tol = 1e-8,
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
  \item{pars}{
    A numeric vector containing the parameters of the model, in order of the 
    immigration parameters \eqn{\mu(.)}, offspring parameters \eqn{h(.)} and 
    lastly the branching ratio \eqn{\eta(.)}.
}
  \item{maxiter}{
  The maximum number of iterations to perform.
  }
  \item{tol}{
  The algorithm stops when the difference between the previous iteration and 
  current iteration parameters sum is less than \code{tol}.
  }

  \item{h.fn}{
    A (vectorized) function. The offspring density function. 
}
  \item{mu.fn}{
    A (vectorized) function. The immigration hazard function. 
}
  \item{H.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the offspring density function from 0 to \code{t}. 
}
  \item{Mu.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the immigrant hazard function from 0 to \code{t}. 
}
  \item{logg.fn}{
    A (vectorized) function. The log of the immigrant distribution function.
  }
}
\value{
  \item{iterations}{The number of iterations until convergence}
  \item{diff}{The absolute sum of the difference between the final two 
  parameter estimates}
  \item{pars}{The parameter estimates from the EM algorithm} 
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>
}

\examples{
\dontrun{
## simulated data
tms <- sort(runif(100,0,100))
## the quicker version on simulated data with default immigrant hazard function
## and offspring density
system.time( est2 <- EM2partial(tms, 101, c(2,1,0.5,1)) )
}
}

\keyword{ RHawkes }
\keyword{ EM algorithm }
