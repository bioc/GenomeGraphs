\name{GeneRegion-class}
\docType{class}
\alias{GeneRegion-class}
\alias{initialize,GeneRegion-method}
\alias{show,GeneRegion-method}

\title{Class "GeneRegion", representing gene structures in a defined genomic region }
\description{Given a start and end position and a chromosome name, all gene structures in this region will be retrieved from Ensembl upon creation of the object.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneRegion", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{start}:}{Object of class \code{"numeric"}, start position}
    \item{\code{end}:}{Object of class \code{"numeric"}, end position }
    \item{\code{chromosome}:}{Object of class \code{"character"}, chromosome name }
    \item{\code{strand}:}{Object of class \code{"character"}, represents the strand from which the gene structures should be retrieved.  Value is either + or -}
    \item{\code{size}:}{Object of class \code{"numeric"}, represents the size of the GeneRegion in the final plot }
    \item{\code{color}:}{Object of class \code{"character"}, represents the color to be used to plot the exons }
    \item{\code{biomart}:}{Object of class \code{"Mart"}, containing the link to the Ensembl database.  This should be created by the useMart function from the biomaRt package}
    \item{\code{ens}:}{Object of class \code{"data.frame"}, output of the biomaRt query, should not be used by users}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GeneRegion")}: ... }
    \item{show}{\code{signature(object = "GeneRegion")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}
\seealso{
       objects to See Also as \code{\link{gdPlot}}
}
\examples{
if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
plusStrand = new("GeneRegion", chromosome = "17", start = 30450000, end = 30550000, strand = "+", biomart = mart)
genomeAxis = new("GenomeAxis", add53=TRUE)
gdPlot(list(genomeAxis, plusStrand), minBase = 30450000, maxBase =  30550000)
}
}
\keyword{classes}
