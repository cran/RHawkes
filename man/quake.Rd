\name{quake}
\alias{quake}
\docType{data}
\title{
An RHawkes earthquake data set
}
\description{
An earthquake data set containing the earthquake occurrence times near the Japan 
region previously examined by Ogata (1998).
}

\usage{data(quake)}
\format{
  The format is a vector of the arrival/birth times of earthquakes.
}
\details{
  Times of arrivals of earthquake occurrences in a vector in ascending order.
}
\source{
  Simulated by a call to the function \code{simHawkes1}.
}
\examples{
\dontrun{
data(quake)
## number of earthquake occurrences
nrow(quake)
}
}
\keyword{datasets}
