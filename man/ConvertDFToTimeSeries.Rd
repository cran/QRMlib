\name{ConvertDFToTimeSeries}
\alias{ConvertDFToTimeSeries}
\title{ ConvertDFToTimeSeries() method }
\description{
Method to convert a data.frame object to a timeSeries object and
insure that the any eight- or nine-character date elements like
1/1/2001 or 1/10/2001 or 10/1/2001
are converted to the ten-character format required by timeSeries.
}
\details{
Insures that the month-day-year format has 2 monthly digits, two daily digits
and 4 annual digits plus the two separating backslashes.  otherwise the
RMetrics timeSeries class will balk at converting a data.frame to a timeSeries
}
\usage{
ConvertDFToTimeSeries(dataframe)
}
\arguments{
  \item{dataframe}{a data.frame object with DATE field in "m/d/Y" format }
}
\value{
  a timeSeries in the RMetrics fSeries package (270.60)
}
\author{ documentation by Scott Ulman for R-language distribution }
\examples{
data(danish.df);
danishTS <- ConvertDFToTimeSeries(danish.df);
save("danishTS",file="danishTS.R");
}
\keyword{ ts }
\keyword{models }
