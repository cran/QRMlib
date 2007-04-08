\name{kurtosisSPlus}
\alias{kurtosisSPlus}

\title{ S-Plus Version of Kurtosis which differs from the R-versions}
\description{
The values calculated by R and S-Plus differ when we use the 
call kurtosis(x, method="moment") which causes serious consequences
in the fit.NH() function call. Hence we introduce the S-Plus 
version here. S-Plus has only the "moment" and "fisher" methods. 
R has a 3rd type, the "excess" which should parallel the R "moment" 
type but fails.
}
\usage{
kurtosisSPlus(x, na.rm = FALSE, method = "fisher")
}
\arguments{
  \item{x}{ data vector }
  \item{na.rm}{TRUE or FALSE indicating whether to remove any NA values 
      from the data vector }
  \item{method}{ either the 'moment' or 'fisher' method }
}
\details{
  use R-code which reflects the way S-Plus calculates Kurtosis so we match
  up the answer regardless of whether using S-Plus or R
}
\value{
  a single number reflecting the kurtosis statistic as calculated via 
  the S-Plus method (either "moment" or "fisher")
}
\author{documentation by Scott Ulman for R-language distribution}
\keyword{methods}
