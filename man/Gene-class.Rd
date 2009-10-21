\name{Gene-class}
\docType{class}
\alias{Gene-class}
\alias{initialize,Gene-method}
\alias{drawGD,Gene-method}
\alias{show,Gene-method}

\title{Class "Gene" represents the Ensembl Gene level annotation }
\description{Class "Gene" represents the Ensembl Gene level annotation.  Upon creation of an object of this class, intron and exon boundaries are retrieved from Ensembl}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Gene", ...)}. 
}
\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"character"}, representing a unique identifier for the gene or a vector of identifiers for genes that are located near each other (or at least on the same chromosome)}
    \item{\code{type}:}{Object of class \code{"character"}, representing the type of identifier used, e.g. hgnc\_symbol, entrezgene and ensembl\_gene\_id, check the listFilters function of the biomaRt package for more identifier options}
%     \item{\code{size}:}{Object of class \code{"numeric"}, specifies the size of the plot}
%     \item{\code{color}:}{Object of class \code{"character"}, specifies the color of the exons }
    \item{\code{biomart}:}{Object of class \code{"Mart"}, contains the link to the Ensembl database and should be created using the useMart function from the biomaRt package }
    \item{\code{ens}:}{Object of class \code{"data.frame"}, contains the output from the Ensembl query, users don't need to give a value to this }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Gene")}: ... }
    \item{drawGD}{\code{signature(.Object = "Gene")}: ... }    
    \item{show}{\code{signature(object = "Gene")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Jim Bullard and Steffen Durinck}

\seealso{
objects to See Also as \code{\link{gdPlot}}
}

\examples{
if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene = new("Gene", id = "ENSG00000095203", type="ensembl_gene_id", biomart = mart)
gdPlot(list(gene), minBase= 110974000, maxBase = 111122900)
}
}
\keyword{classes}
