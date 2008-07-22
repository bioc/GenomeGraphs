\name{makeGeneRegion}
\alias{makeGeneRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class Gene containing the intron-exon structures of genes}
\description{
Creates an object of class Gene containing the intron-exon structures of genes.
Given a start and end position, strand and chromosome, all the intron-exon strcutures of all genes laying in this region will be retrieved.
}
\usage{
makeGeneRegion(start, end, chromosome, strand, biomart, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{Start position on chromosome}
  \item{end}{End position on chromosome}
  \item{chromosome}{Chromosome name}
  \item{strand}{Strand either + or -}
  \item{biomart}{Mart object, created by the useMart function of biomaRt}
  \item{dp}{Object of class DisplayPars, determines the display of features on the plot}
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
function (start, end, chromosome, strand, biomart, dp = NULL) 
{
    if (missing(start)) 
        stop("Need to specify a start for creating a GeneRegion")
    pt <- getClass("GeneRegion")@prototype
    if (is.null(dp)) 
        dp <- pt@dp
    if (is.numeric(chromosome)) 
        chromosome = as.character(chromosome)
    new("GeneRegion", start = start, end = end, chromosome = chromosome, 
        strand = strand, biomart = biomart, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
