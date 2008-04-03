\name{DisplayPars-class}
\docType{class}
\alias{DisplayPars-class}
\alias{initialize,DisplayPars-method}
\alias{show,DisplayPars-method}
\alias{getPar,DisplayPars-method}
\alias{setPar,DisplayPars-method}

\title{Class "DisplayPars" is used to specify graphical parameters to
  gdObjects.
}
\description{The DisplayPars functions analagously to par and
  gp. Generally the class is instantiated using the \code{DisplayPars}
  function rather than directly.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{DisplayPars(...)}
  rather than calls to: \code{new("DisplayPars", ...)} by calling the
  \code{DisplayPars} function directly in the constructor the gdObjects
  are guaranteed to have the appropriate defaults.
}
\section{Slots}{
  \describe{
    \item{\code{pars}:}{Object of class \code{"environment"} Generally
      this slot is not accessed directly.
    }
  }
}
\section{Methods}{
  \describe{
    \item{getPar}{\code{signature(obj = "DisplayPars")}: gets a
      graphical parameter by name}
    \item{initialize}{\code{signature(.Object = "DisplayPars")}: This
      constructor should not be called directly.}
    \item{setPar}{\code{signature(obj = "DisplayPars")}: sets a
      graphical parameter - see the example below. Often it is easier to set
      the graphical parameter from within the gdObject.}
    \item{show}{\code{signature(object = "DisplayPars")}: prints current
      graphical parameters}
  }
}
\section{Warning}{
  The \code{DisplayPars} class should not be manipulated directly. The
  preferred method for interacting with the class can be seen in the
  example below. 
}
\author{James Bullard}
\examples{
showClass("DisplayPars")

if (interactive()) {
minbase = 10000
maxbase = 15000
mart <- useMart("ensembl", dataset = "scerevisiae_gene_ensembl")
genesplus <- new("GeneRegion", start = minbase, end = maxbase, biomart = mart,
                strand = "+", chromosome = "I", dp = DisplayPars(size = 2))
## plot it. 
gdPlot(list(genesplus, new("Title", title = "genes")), minbase, maxbase)

## to obtain a list of the current graphical parameters:
print(genesplus@dp)

## to set a parameter:
setPar(genesplus, "protein_coding", "pink")
gdPlot(list(genesplus, new("Title", title = "genes")), minbase, maxbase)
}
}
\keyword{classes}
