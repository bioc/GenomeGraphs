\name{Transcript-class}
\docType{class}
\alias{Transcript-class}
\alias{initialize,Transcript-method}
\alias{show,Transcript-method}
\alias{drawGD,Transcript-method}

\title{Represent known transcript isoforms as annoted by Ensembl}
\description{Represent known transcript isoforms as annoted by Ensembl}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Transcript", ...)}. 
}
\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"character"}, represents the gene identifier that should be used to retrieve the tanscript level annotation }
    \item{\code{type}:}{Object of class \code{"character"}, represents the type of identifiers used to specify the gene e.g. hgnc\_symbol, entrezgene and ensembl\_gene\_id}
%     \item{\code{size}:}{Object of class \code{"numeric"}, represents the size of the plot of this object}
    \item{\code{transcriptSize}:}{Object of class \code{"numeric"}, represents the size of the transcripts in the plot}
    \item{\code{numOfTranscripts}:}{Object of class \code{"numeric"}, should not be used bu users }
%     \item{\code{color}:}{Object of class \code{"character"}, color of the exons}
    \item{\code{biomart}:}{Object of class \code{"Mart"}, containing the links to the Ensembl database.  This object should be created with the useMart function of the biomaRt package }
    \item{\code{ens}:}{Object of class \code{"data.frame"}, should not be used by the users.  Contains the output from the biomaRt query }
  }
}
\section{Methods}{
  \describe{
  \item{drawGD}{\code{signature(.Object = "Transcript")}: ... }    
    \item{initialize}{\code{signature(.Object = "Transcript")}: ... }
    \item{show}{\code{signature(object = "Transcript")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}

\seealso{
	objects to See Also as \code{\link{gdPlot}}

}
\examples{
if(interactive()){
data("unrData", package="GenomeGraphs")
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
transcript = new("Transcript", id ="ENSG00000009307" , biomart = mart)
gdPlot(list(transcript), minBase = min(exon@probeStart), maxBase=max(exon@probeEnd))
}
}
\keyword{classes}
