\name{sim.pred1}
\alias{sim.pred1}
\title{
  Simulate a fitted RHawkes process model for prediction purposes
}
\description{
Simulate a fitted RHawkes process model from the censoring time \code{cens} 
to a future time point \code{cens.tilde}, conditional on the observed
event times until the censoring time.
}
\usage{
sim.pred1(tms, par, re.dist = rweibull, 
    par.redist = list(shape = par[1], scale = par[2]), 
    of.dis="exp", par.ofdis = list(rate=par[3]), 
    branching.ratio = par[4], cens=tail(tms,1)+mean(diff(tms))/2, 
    cens.tilde = cens * 1.5, 
    mu.fn = function(x, p) {
        exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
        pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))
    })
}
\arguments{
  \item{tms}{
    A numeric vector, with values sorted in ascending order. Event times
    to fit the RHawkes point process model.}
  \item{par}{
    A numeric vector, giving the parameters of the model with the 
    immigration parameters \eqn{\mu(.)}, offspring parameters \eqn{h(.)} and lastly
    the branching ratio \eqn{\eta(.)}.
  }
  \item{re.dist}{
    A (vectorized) function. The function to simulate from the immigrant
    waiting times distribution.
  }
  \item{par.redist}{
    A (named) list, giving the parameters of the immigrant waiting time distribution.  
  }
  \item{of.dis}{
    A character string, for the name of the offspring birth time distribution. 
  }
  \item{par.ofdis}{
    A (named) list, giving the parameters of the offspring birth time distribution.
  }
  \item{branching.ratio}{
    A scalar in [0,1), the branching ratio parameter.    
  }
  \item{cens}{
    A scalar. The censoring time. 
  }
  \item{cens.tilde}{
    A scalar. The future time to run the simulation to. 
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

\seealso{
  \code{\link{sim.pred}}.
}

\examples{
N <- 5; i <- 0;
data(quake); tms <- sort(quake$time);
# add some random noise the simultaneous occurring event times
tms[213:214] <- tms[213:214] + 
                    sort(c(runif(1, -1, 1)/(24*60), runif(1, -1, 1)/(24*60)))
# simulate future event time based on MLE fitted RHawkes model
times <- replicate(N,
                   {cat(i<<-i+1,'\n');
                   sim.pred1(tms = tms, par = c(0.314, 22.2, 1266, 0.512),
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
