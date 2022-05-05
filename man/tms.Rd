\name{tms}
\alias{tms}
\docType{data}
\title{
mid-price change times of the AUD/USD exchange rate
}
\description{
A financial data set containing the mid-price changes of the AUD/USD foreign
exchange rate during the trading week from 20:00:00 GMT on Sunday 19/07/2015 to 
21:00:00 GMT Friday 24/07/2015.
}
\usage{data(tms)}
\format{
  The format is a list of the arrival times of mid price changes that occur 
  every hour in 121 non-overlapping windows.
}
\details{
  Times of arrivals of mid-price changes is listed
  together in ascending order.
}
\source{
  Simulated by a call to the function \code{simHawkes1}.
}
\examples{
data(tms)
## number of non over-lapping hourly windows
length(tms)
}
\keyword{datasets}
