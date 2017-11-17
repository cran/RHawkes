\name{sim.pred}
\alias{sim.pred}
\title{
Simulate a fitted RHawkes process model
}
\description{
Simulate a fitted RHawkes process model after the censoring time \code{cens} 
to a future time point \code{cens.tilde}.
}
\usage{
sim.pred(tms, re.dist = rweibull, par, 
    par.redist = list(shape = par[1], scale = par[2]), 
    h.fn = function(x, p) dexp(x, rate = 1 / p), p.ofd = par[3], 
    branching.ratio = par[4], cens, cens.tilde = cens * 1.5, 
    mu.fn = function(x, p) {
        exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
        pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
    })
}
\arguments{
  \item{tms}{
    A numeric vector, with values sorted in ascending order. Event times
    to fit the RHawkes point process model.}
  \item{re.dist}{
    A (vectorized) function. The immigrant renewal distribution function.
}
  \item{par}{
    A numeric vector, giving the parameters of the model with the 
    immigration parameters \eqn{\mu(.)}, offspring parameters \eqn{h(.)} and lastly
    the branching ratio \eqn{\eta(.)}.
}
  \item{par.redist}{
    A numeric vector. The parameters of the immigrant renewal distribution.  
}
  \item{h.fn}{
    A (vectorized) function. The offspring density function. 
}
  \item{p.ofd}{
    A numeric vector. The paramters of the offspring density.
}
  \item{branching.ratio}{
    A scalar. The branching ratio parameter.    
}
  \item{cens}{
    A scalar. The censoring time. 
}
  \item{cens.tilde}{
    A scalar. The future time that the simulation run until. 
}
  \item{mu.fn}{
    A (vectorized) function. The immigration hazard function. 
 }
}

\value{
A numeric vector that contains the simulated event times from censoring time 
\code{cens} up until \code{cens.tilde}
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>
}


\examples{
N <- 5; i <- 0;
data(quake); tms <- sort(quake$time);
# add some random noise the simultaneous occuring event times
tms[213:214] <- tms[213:214] + 
                    sort(c(runif(1, -1, 1)/(24*60), runif(1, -1, 1)/(24*60)))
# simulate future event time based on MLE fitted Rhawkes model
times <- replicate(N,
                   {cat(i<<-i+1,'\n');
                   sim.pred(tms = tms, par = c(0.314, 22.2, 1266, 0.512),
                            cens=35063)
                   })
plot(NA,NA,xlim=c(0,35063*1.5),ylim=c(0,max(lengths(times))+nrow(quake)),
     xlab="time",ylab="Sample path")
lines(c(0,quake$time),0:nrow(quake),type="s")
for(i in 1:N)
    lines(c(tail(quake$time,1),times[[i]]),nrow(quake)+0:length(times[[i]]),
          type="s",lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ predict }
\keyword{ simulations }% __ONLY ONE__ keyword per line
