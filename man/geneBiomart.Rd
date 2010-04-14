\name{geneBiomart}
\alias{geneBiomart}
\title{AnnotationTrack objects from biomaRt}
\description{
  Convenience function to construct an AnnotationTrack object from biomaRt.
}
\usage{
geneBiomart(id, biomart, type = "ensembl_gene_id", dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{Gene identifier}
  \item{biomart}{Mart object connected to BioMart database, use useMart function to generate}
  \item{type}{Type of identifier used, this should be a filter of the BioMart database e.g. ensembl_gene_id, hgnc_symbol}
  \item{dp}{Display parameters}
}
\value{
  An AnnotationTrack object
}
\author{James Bullard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

