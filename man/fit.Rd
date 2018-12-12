\name{Fit}
\alias{PLSDA}
\alias{fit}
\alias{print.PLSDA}
\alias{predict.linmod}
\alias{summary.PLSDA}
\title{PLS fit}
\description{The function fit performs partial least squares regression for the multivariate case (i.e. more than one response variable)}
\usage{ fit(formula,data,ncomp=2,scale=T)
\method{fit}{default}(formula, data, ...)
\method{print}{PLSDA}(x, ...)
\method{summary}{PLSDA}(object, ...)
\method{predict}{PLSDA}(object, newdata=NULL, ...) }
\arguments{
\item{formula}{ a symbolic description of the model to be fit. }
\item{data}{ an  data frame containing explivatives and explained variables in the model. }
\item{scale}{ boleean. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE) }
\item{object}{ an object of class
\code{"PLSDA"}, i.e., a fitted model. }
\item{\dots}{ not used. }
}
\value{ An object of class
\code{PLSDA}, basically a list including elements
\item{x.scores}{ components of the predictor variables (also known as T-components)}
\item{x.loads}{ loadings of the predictor variables }
\item{y.scores}{ components of the response variables (also known as U-components) }
\item{y.loads}{ loadings of the response variables }
\item{expvar}{ table with R-squared coefficients }}

\author{Houda Quentin Faly}
\examples{  data(iris)
  mod1 <- fit(Hwt~Bwt*Sex, data=iris)
  mod1
  summary(mod1) }
\keyword{regression}
