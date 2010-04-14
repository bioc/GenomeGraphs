\name{geneRegionBiomart}
\alias{geneRegionBiomart}
\title{Construct an AnnotationTrack object from biomaRt.}
\description{
  This function constructs an AnnotationTrack object from Biomart. It is
  a convenience function.
}
\usage{
geneRegionBiomart(chr, start, end, strand, biomart, dp = NULL, chrFunction = function(x) x, strandFunction = function(x) x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{\code{chr} An integer}
  \item{start}{\code{start} The start location }
  \item{end}{\code{end} The end location }
  \item{strand}{\code{strand} An integer -1, 0, 1}
  \item{biomart}{\code{biomart} A mart}
  \item{dp}{\code{dp} DisplayPars object}
  \item{chrFunction}{\code{chrFunction} A function which takes
    \code{chr} and converts it into the correct representation for
    biomaRt. For instance yeast likes to have chromosomes as roman
    numerals so you can use as.roman here.}
  \item{strandFunction}{\code{strandFunction} Analagous to chrFunction,
    but for strand. The internal representation of strand is -1,0,1
    whereas biomaRt has different species dependent representations.}
}
\value{
  An AnnotationTrack object. 
}
\author{James Bullard}
\keyword{hplot}
