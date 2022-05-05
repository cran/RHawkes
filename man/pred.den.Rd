\name{pred.den}
\alias{pred.den}
\title{
RHawkes predictive density function
}
\description{
Calculates the predictive density of the next observed event time after the 
censoring time \code{cens} based on observations over the interval 
[0,\code{cens}].
}
\usage{
pred.den(x, tms, cens, par, 
        h.fn = function(x, p) dexp(x, rate = 1 / p), 
        mu.fn = function(x, p) {
          exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
          pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))}, 
          H.fn = function(x, p) pexp(x, rate = 1 / p), 
          Mu.fn = function(x, p) {
            -pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE)
          })
}
\arguments{
  \item{x}{
    A scalar. The amount of time after the censoring time \code{cens}.
}
  \item{tms}{
    A numeric vector, with values sorted in ascending order. The event times
    to fit the RHawkes point process model.
}
  \item{cens}{
    A scalar. The censoring time.
}
  \item{par}{
    A numeric vector. Contains the parameters of the model, in order of the 
    immigration parameters \eqn{\mu(.)}, offspring parameters \eqn{h(.)} and 
    lastly the branching ratio \eqn{\eta(.)}.
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
}

\value{
The predictive density of the next event evaluated at \code{x}.
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>
  }
  
\examples{
data(quake);
tms <- sort(quake$time);
# add some random noise to the identical event times
tms[213:214] <- tms[213:214] + 
                    sort(c(runif(1, -1, 1)/(24*60), runif(1, -1, 1)/(24*60)))
curve(pred.den(x, tms = tms, cens = 35064, par= c(0.314, 22.2, 1266, 0.512)) 
        ,0 ,2000, col = 2, lty = 2)
}

\keyword{ prediction }
\keyword{ point process }% __ONLY ONE__ keyword per line
