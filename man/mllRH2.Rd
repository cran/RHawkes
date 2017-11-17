\name{mllRH2}
\alias{mllRH2}
\title{
Minus loglikelihood of a RHawkes model with Rosenblatt residuals}
\description{
Calculates the minus loglikelihood of a RHawkes model with given
immigration hazard function \eqn{\mu}, offspring density function 
\eqn{h} and bracnhing ratio \eqn{\eta} for event times \code{tms} 
on interval \eqn{[0,cens]}. The same as \code{mllRH} although this
version also returns the Rosenblatt residuals.
}
\usage{
mllRH2(tms, cens, par, h.fn = function(x, p) dexp(x, rate = 1/p), 
       mu.fn = function(x, p) {
         exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) - 
         pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, log.p = TRUE))}, 
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
\details{
Calculate the RHawkes point process Rosenblatt residuals
}
\value{
\item{mll}{minus log-likelihood}
\item{U}{Rosenblatt residual of observed event time}
\item{n}{number of events}
}

\author{
  Feng Chen <feng.chen@unsw.edu.au>
  Tom Stindl <t.stindl@unsw.edu.au>}




\seealso{
\code{mllRH}}
\examples{
tmp <- mllRH2(sort(runif(50,0,100)),100,c(1,2,2,0.05))
qqunif<-function(dat,...){
  dat<-sort(as.numeric(dat));
  n<-length(dat);
  pvec<-ppoints(n);
  plot(pvec,dat,xlab="Theoretical Quantiles",
       ylab="Sample Quantiles",main="Uniform Q-Q Plot",...)
}
par(mfrow=c(1,2))
qqunif(tmp$U)
acf(tmp$U)
ks.test(tmp$U,"punif")
}
\keyword{ Rosenblatt residual }
\keyword{ goodness-of-fit }
\keyword{ point process }% __ONLY ONE__ keyword per line
