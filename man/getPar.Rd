\name{getPar}
\alias{getPar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Retrieves a display parameter from an object.}
\description{
  Retrieves a display parameter from an object. 
}
\usage{
getPar(obj, name, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A gdObject or DisplayPars object.}
  \item{name}{Name of parameter to return.}
  \item{...}{Ignored}
}
\examples{
a <- new("GenomeAxis")
getPar(a, "size")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

