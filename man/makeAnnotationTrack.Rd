\name{makeAnnotationTrack}
\alias{makeAnnotationTrack}
\title{Create objects of class AnnotationTrack}
\description{
  Convenience function for constructing objects of class
  AnnotationTrack. 
}
\usage{
makeAnnotationTrack(regions = NULL, chr = NULL, strand = NULL, start = NULL, end = NULL, feature = NULL, group = NULL, ID = NULL, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regions}{A dataframe with columns start, end, feature,
    group, ID. start and end delineate the boundaries of the boxes
    feature can be used to color the boxes. Group denotes linking so
    generally exons from a gene form a group. Finally, ID can be used to
    plot names on boxes.}
  \item{chr}{The chromosome of the regions (can be ignored)}
  \item{strand}{The strand of the regions (can be ingored) }
  \item{start}{If regions is missing then we construct a dataframe from
   the remaining parameters.}
  \item{end}{Construct regions with this vector}
  \item{feature}{Construct regions with this feature vector or scalar}
  \item{group}{Defines a grouping }
  \item{ID}{Defines an ID for each annotation bit }
  \item{dp}{DisplayPars, in this case we can create a mapping between
    feature and color. So lets say in the feature column you have:
    gene, transcript, gene, then in the dp you can say gene = 'blue' and
    transcript = 'green'}
}
\details{
  
}
\value{
  Returns an object of class AnnotationTrack
}
\examples{
a <-  makeAnnotationTrack(start = c(10, 15, 25), end = c(12, 19, 31),
                      group = c(1,1,2), feature = c("gene", "gene", "tf"),
                      ID = paste("id", 1:3, sep = ""), dp = DisplayPars(gene = 'blue'))
gdPlot(a, minBase = 0, maxBase = 40)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

