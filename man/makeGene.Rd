\name{makeGene}
\alias{makeGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class Gene}
\description{
Creates an object of class Gene.  This represents a gene structure as annotated in Ensembl.
}
\usage{
makeGene(id, type, biomart, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{An identifier used to specify of which gene the intron-exon structure should be retrieved}
  \item{type}{The type of identifiers used, examples are ensembl\_gene\_id, hgnc\_symbol,entrezgene. See listAttributes function of the biomaRt package for more info}
  \item{biomart}{Mart object, created by the useMart function of biomaRt}
  \item{dp}{object of class DisplayPars, determines the display of features on the plot}
}
\value{
An object of class Gene
}
\references{ ~put references to the literature/web site here ~ }
\author{Steffen Durinck and Jim Bullard}
\seealso{\code{\link{gdPlot}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (id, type, biomart, dp = NULL) 
{
    if (missing(id)) 
        stop("Need to specify a gene identifier for creating a Gene")
    pt <- getClass("Gene")@prototype
    if (is.null(dp)) 
        dp <- pt@dp
    if (missing(type)) 
        type = pt@type
    new("Gene", id = id, type = type, biomart = biomart, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
