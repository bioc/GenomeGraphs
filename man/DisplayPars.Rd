\name{DisplayPars}
\alias{DisplayPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DisplayPars constructs objects of type DisplayPars which are used
  to effect the display of gdObjects}
\description{
  DisplayPars takes any number of named arguments which will be used by
  the drawGD method of the gdObject. These arguments are analagous to
  both par and gp of the traditional and grid graphics systems
  respectively. Different functions support different graphical
  parameters - thus it is necessary to consult the documentation of the
  particular gdObject to determine which DisplayPars will be used.
}
\usage{
DisplayPars(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ name value pairs }
}
\value{
  Returns an object of type DisplayPars, generally this will be called
  during a call to the \code{new} function for a particular gdObject. 
}
\details{
  It is not recommended to call \code{new("DisplayPars", ...)} directly;
  rather this function \code{DisplayPars()} should be called instead. If
  a gdObject has already been instantiated then the appropriate method
  for changing graphical parameters is: \code{setPar}. 
}
\author{ James Bullard }
\examples{
minbase = 10000
maxbase = 15000
mart <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
genesplus <- new("GeneRegion", start = minbase, end = maxbase, biomart = mart,
                strand = "+", chromosome = "I", dp = DisplayPars(color =
                "red"))
gaxis <- new("GenomeAxis", add53 = TRUE, add35 = TRUE)
genesminus <- new("GeneRegion", start = minbase, end = maxbase, biomart = mart,
                strand = "-", chromosome = "I", dp = DisplayPars(color =
                "purple", size = 2))
title <- new("Title", title = "genes in a region")
gdPlot(list(genesplus, gaxis, genesminus, title), minbase, maxbase)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}


