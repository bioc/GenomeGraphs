\name{drawGD}
\alias{drawGD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generic called on each gdObject to do the plotting.}
\description{
  This generic gets called on each of the gdObjects (obviously the
  generic is implemented by a method for each object) and thus if a user
  wishes to implement new gdObjects they need to have access to this
  generic. 
}
\usage{
drawGD(gdObject, minBase, maxBase, vpPosition, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gdObject}{}
  \item{minBase}{}
  \item{maxBase}{}
  \item{vpPosition}{}
  \item{...}{Ignored}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
