\name{simRHawkes}
\alias{simRHawkes}
\title{
Simulate a renewal Hawkes (RHawkes) process
}
\description{
Simulate a renewal Hawkes (RHawkes) process with given renewal immigration 
distribution function, offspring density function and branching ratio. 
}
\usage{
simRHawkes(re.dist = rweibull, par.redist = list(shape = 1, scale = 1), 
           ofspr.den = function(x, p.ofd) 1 / p.ofd * exp(-x / p.ofd), 
           p.ofd = 1, branching.ratio = 0.5, cens = 1, B = 10, B0 = 50, 
           max.ofspr.den = max(optimize(ofspr.den, c(0, cens), maximum = TRUE, 
           p = p.ofd)$obj, ofspr.den(0, p.ofd), ofspr.den(cens, p.ofd)) * 1.1)
}
\arguments{
  \item{re.dist}{
    A (vectorized) function. The immigrant renewal distribution function.
}
  \item{par.redist}{
    A numeric vector. The parameters of the immigrant renewal distribution.  
}
  \item{ofspr.den}{
    A (vectorized) function. The offspring density function.    
}
  \item{p.ofd}{
    A numeric vector. The parameters of the offspring density.
}
  \item{branching.ratio}{
    A scalar. The branching ratio parameter.    
}
  \item{cens}{
    A scalar. The censoring time. 
}
  \item{B}{
    A scalar. Tuning parameter for simulation of further immigrants.  
}
  \item{B0}{
    A scalar. Tuning parameter for simulation of initial immigrants.
}
  \item{max.ofspr.den}{
    A scalar. The maximum value of the offspring density function from 0 
    to \code{cens}.
}
}
\details{
  The function works by simulating the arrival times of immigrants
  according to the renewal immigration distribution. The birth times of 
  offspring from each immigrant are then simulated according to an 
  inhomogeneous Poisson processes with appropriate intensity functions.
}
\value{
  A numeric vector of all pooled events (immigration/birth) times of
  all generations 0, 1, ...
}

\author{
   Feng Chen <feng.chen@unsw.edu.au>
   Tom Stindl <t.stindl@unsw.edu.au>
}

\examples{
B <- 10; i <- 0;
tms <- replicate(B,
                  {cat(i<<-i+1,'\n');
                  simRHawkes(par.redist = list(shape = 3, scale = 1),
                      p.ofd = 0.5, branching.ratio = 0.5,
                      cens = 100)
                  })
}
\keyword{renewal process}
\keyword{self-exciting}
