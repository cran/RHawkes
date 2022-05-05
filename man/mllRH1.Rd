\name{mllRH1}
\alias{mllRH1}
\title{
Minus loglikelihood of a RHawkes model with parent probabilities
}
\description{
Calculates the minus loglikelihood of a RHawkes model with given
immigration hazard function \eqn{\mu}, offspring density function 
\eqn{h} and branching ratio \eqn{\eta} for event times \code{tms} 
on interval \eqn{[0,cens]}. The same as \code{mllRH} although this
version also returns the parent probabilities.
}
\usage{
mllRH1(tms, cens, par, 
      h.fn = function(x, p) dexp(x, rate = 1/p), 
      mu.fn = function(x, p) {
        exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
        pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
        }, 
        H.fn = function(x, p) pexp(x, rate = 1/p), 
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
\item{mll}{minus log-likelihood}
\item{log.p}{parent probabilities}
\item{n}{number of events}
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>}

\seealso{
\code{mllRH}}
\examples{
tmp <- mllRH1(sort(runif(1000,0,1000)), 1001, c(2,1,0.5,1))
for(i in 1:tmp$n)
  cat(exp(tmp$log.p[i*(i - 1)/2 + 1:i]), "\n")    
}
\keyword{ parent probabilities }
\keyword{ likelihood }
