\name{simRHawkes1}
\alias{simRHawkes1}
\title{
Simulate a renewal Hawkes (RHawkes) process
}
\description{
Simulate a renewal Hawkes (RHawkes) process with given renewal immigration 
distribution function, offspring density function and branching ratio. 
}
\usage{
simRHawkes1(re.dist = rweibull, par.redist = list(shape = 1, scale = 1), 
            of.dis = "exp", par.ofdis = list(rate=1),
            branching.ratio = 0.5, cens = 1, B = 10, B0 = 50, 
            flatten=TRUE)
}
\arguments{
  \item{re.dist}{
    A (vectorized) function. The immigrant renewal distribution function.
}
  \item{par.redist}{
    A numeric vector. The parameters of the immigrant renewal distribution.  
}
  \item{of.dis}{
    A character string indicating the distribution for the offspring
    birth times, which has to be the
    'distname' part of the \code{rdistname} functions for simulating
    postive random variables implemented in R, such as "exp", "weibull",
    "gamma", "lnorm", etc.
}
  \item{par.ofdis}{
    A list with named elements, giving the list of parameters of the
    offspring distribution, such as list(rate=1), list(shape=1,scale=1),
    etc. 
}
  \item{branching.ratio}{
    A scalar between 0 and 1, the branching ratio parameter. 
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
  \item{flatten}{
    A boolean scalar, which indicates whether the output events times
    should be flattened into an increasing sequence of times, or not (in
    which case the output is the immigrant arrival times, and the
    offspring birth times for different immigrants).
}
}
\details{
  The function works by simulating the arrival times of immigrants
  according to the renewal immigration distribution. The birth times of 
  offspring from each immigrant are then simulated according to an 
  inhomogeneous Poisson processes with appropriate intensity functions.
}
\value{
  A numeric vector of pooled event (immigration/offspring birth) times of
  all generations 0, 1, ..., if \code{flatten=TRUE}; A list with two
  components: \code{immitimes} for immigrant arrival times, and
  \code{offspringtimes} for birth times of offspring due to different
  immigrants. 
}

\author{
   Feng Chen <feng.chen@unsw.edu.au>
   Tom Stindl <t.stindl@unsw.edu.au>
}

\examples{

tms <- simRHawkes1(par.redist = list(shape = 3, scale = 1),
                   par.ofdis = list(rate=0.5), branching.ratio = 0.5,
                   cens = 50)
plot(stepfun(tms,0:length(tms)),do.points=FALSE,vertical=FALSE,xlim=c(0,50))

tms.clust <- simRHawkes1(par.redist = list(shape = 3, scale = 1),
                         par.ofdis = list(rate=0.5), branching.ratio = 0.5,
                         cens = 50,flatten=FALSE)
plot(c(0,50),c(0, 1+(nt <-length(it <- tms.clust$immitimes))),
     type="n",xlab="time",ylab="cluster")
segments(x0 = it,y0=-0.2,y1=0.2)
for(i in 1:nt)
    segments(x0 = c(it[i],it[i]+tms.clust$offspringtimes[[i]]),
             y0=i-0.2,y1=i+0.2)
abline(h=0:(nt+1),col="light gray",lty=2)
segments(x0=unlist(lapply(1:nt,function(i)c(it[i],it[i]+tms.clust$offspringtimes[[i]]))),
         y0=nt+1-0.2,y1=nt+1+0.2)
}
\keyword{renewal process}
\keyword{self-exciting}
