\name{setPar}
\alias{setPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sets a display parameter}
\description{
  Sets a display parameter
}
\usage{
setPar(obj, name, val, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object, usually a gdObject.}
  \item{name}{Name of display parameter to set.}
  \item{val}{Value of display parameter.}
  \item{...}{Ignored}
}
\examples{
a <- new("GenomeAxis")
setPar(a, "size", 100)
gdPlot(a, minBase = 10, maxBase = 10000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

