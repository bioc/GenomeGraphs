\name{TranscriptRegion-class}
\docType{class}
\alias{TranscriptRegion-class}
\alias{show,TranscriptRegion-method}

\title{Class "TranscriptRegion", representing  a genomic region with transcripts }
\description{Upon creation of this object, transcripts present in a specified region will be retrieved from Ensembl}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TranscriptRegion", ...)} 
}
\section{Slots}{
	 \describe{
    \item{\code{start}:}{Object of class \code{"numeric"}, the start base of the genomic region}
    \item{\code{end}:}{Object of class \code{"numeric"}, the end base of the genomic region }
    \item{\code{chromosome}:}{Object of class \code{"character"}, the chromosome}
%     \item{\code{size}:}{Object of class \code{"numeric"}, the size of the plot for this object}
    \item{\code{biomart}:}{Object of class \code{"Mart"}, contains link to Ensembl and should be created using the useMart function of the biomaRt package }
    \item{\code{ens}:}{Object of class \code{"data.frame"}, users should not specify this, it contains the output of the query to Ensembl }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "TranscriptRegion")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}

\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
showClass("TranscriptRegion")
}
\keyword{classes}
