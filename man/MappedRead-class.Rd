\name{MappedRead-class}
\docType{class}
\alias{MappedRead-class}
\alias{show,MappedRead-method}

\title{Represents mapped reads}
\description{ Represents mapped reads}
\section{Slots}{
         \describe{
    \item{\code{start}:}{Object of class \code{"numeric"}, containing start position of the reads}
    \item{\code{end}:}{Object of class \code{"numeric"}, containing end position of the reads}
    \item{\code{strand}:}{Object of class \code{"numeric"}, containing strand  to which the reads map}
    \item{\code{chromosome}:}{Object of class \code{"numeric"}, containing chromosome to which the reads map}

   }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MappedRead")}: ... }
         }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}



\examples{
## maybe str(MappedRead) ; plot(MappedRead) ...
}
\keyword{datasets}
