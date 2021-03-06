% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sens_spec_fmodel.R
\name{senspec}
\alias{senspec}
\alias{senspec.default}
\alias{senspec.afmodel}
\title{S3 methods to estimate diagnosis performance of an afmodel}
\usage{
senspec(object, ...)

\method{senspec}{default}(object, ...)

\method{senspec}{afmodel}(object, cutoff, ...)
}
\arguments{
\item{object}{with the data to calculate the sensitivity and specificity}

\item{...}{other parameters for the implementing functions}

\item{cutoff}{vector of cut-off points to make the estimations}
}
\value{
a matrix with the columns sensitivity and specificity,
ppv (positive predicted value) and npv (negative predicted value)

No return value. Raise an error.

a matrix with the columns sensitivity and specificity,
ppv (positive predicted value) and npv (negative predicted value)
}
\description{
Estimate sensitivity, specificity, positive predicted value and
negative predicted value negative predictive value from an afmodel.
The estimated "true" negative and "true" positive are estimated using
the estimated overall attributable fraction and the predictive positive value
associated with each cut-off point as described by
Smith, T., Schellenberg, J.A., Hayes, R., 1994.
Attributable fraction estimates and case definitions for malaria
in endemic areas. Stat Med 13, 2345–2358.
}
\examples{
{
# Get the sample data
head(malaria_df1)
fit <- logitexp(malaria_df1$fever, malaria_df1$density)
fit
senspec(fit,  c(1,100,500,1000,2000,4000,8000,16000, 32000,54000,100000))
}
}
\seealso{
\code{\link{logitexp}}
}
