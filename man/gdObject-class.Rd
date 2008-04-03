\name{gdObject-class}
\docType{class}
\alias{gdObject-class}
\alias{getCex,gdObject-method}
\alias{getColor,gdObject-method}
\alias{getLty,gdObject-method}
\alias{getLwd,gdObject-method}
\alias{getPch,gdObject-method}
\alias{getPointSize,gdObject-method}
\alias{getSize,gdObject-method}
\alias{initialize,gdObject-method}
\alias{getPar,gdObject-method}
\alias{setPar,gdObject-method}

\title{Class "gdObject" is the parent class of all of the objects in the
system.}
\description{The gdObject is the superclass of all the classes in the
  system and provides some basic functionality for displaying and
  managing graphical parameters.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("gdObject",
    ...)}. Generally, this class is meant to be subclassed and not
  created directly. 
}
\section{Slots}{
  \describe{
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getCex}{\code{signature(obj = "gdObject")}: ... }
    \item{getColor}{\code{signature(obj = "gdObject")}: ... }
    \item{getLty}{\code{signature(obj = "gdObject")}: ... }
    \item{getLwd}{\code{signature(obj = "gdObject")}: ... }
    \item{getPar}{\code{signature(obj = "gdObject")}: ... }
    \item{getPch}{\code{signature(obj = "gdObject")}: ... }
    \item{getPointSize}{\code{signature(obj = "gdObject")}: ... }
    \item{getSize}{\code{signature(obj = "gdObject")}: ... }
    \item{initialize}{\code{signature(.Object = "gdObject")}: ... }
    \item{setPar}{\code{signature(obj = "gdObject")}: ... }
  }
}
\author{James Bullard}
\examples{
showClass("gdObject")
}
\keyword{classes}
