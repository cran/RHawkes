\name{quake}
\alias{quake}
\docType{data}
\title{
  Japanese Earthquake Data
}
\description{
  Times and magnitudes (Richter scale) of 483 shallow earthquakes in
  an offshore region east of Honshu and south of Hokkaido, for the
  period 1885-1980.
}
\usage{data("quake")}
\format{
  A data frame with 483 observations on the following 2 variables.
  \describe{
        \item{\code{time}}{a numeric vector giving earthquake times in
	  days since start of 1885}
	    \item{\code{mag}}{a numeric vector giving the magnitudes of
    the
    earthquakes}
}
}
\details{
  Times of arrivals of earthquake occurences in a vector in ascending order.
}
\source{
  Simulated by a call to the function \code{simHawkes1}.
}
\examples{
data(quake)
summary(quake)
}
\keyword{datasets}

