\name{seqDataEx}
\alias{seqDataEx}
\docType{data}
\title{This is an example data set from chromosome 4 of yeast from
  various publicly available datasets.}
\description{
  This was a small dataset constructed from publicly available
  datasets. Please see references for details.
}
\usage{data(seqDataEx)}
\format{
  The format is:
List of 5
 $ nislow      : num [1:1751, 1:4] 4 4 4 4 4 4 4 4 4 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1751] "1" "2" "3" "4" ...
  .. ..$ : chr [1:4] "chr" "location" "strand" "evalue"
 $ conservation: num [1:5829, 1:4] 4 4 4 4 4 4 4 4 4 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:5829] "1" "2" "3" "4" ...
  .. ..$ : chr [1:4] "chr" "location" "strand" "score"
 $ david       : num [1:1633, 1:4] 4 4 4 4 4 4 4 4 4 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1633] "1" "2" "3" "4" ...
  .. ..$ : chr [1:4] "chr" "location" "strand" "expr"
 $ snyder      : int [1:1392, 1:4] 4 4 4 4 4 4 4 4 4 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1392] "1" "2" "3" "4" ...
  .. ..$ : chr [1:4] "chr" "location" "strand" "counts"
 $ miura       : num [1:10282, 1:4] 4 4 4 4 4 4 4 4 4 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10282] "1" "2" "3" "4" ...
  .. ..$ : chr [1:4] "chr" "location" "strand" "transcrid"
}
\source{
  @article{Nagalakshmi:2008ve,
	Author = {Nagalakshmi, Ugrappa and Wang, Zhong and Waern, Karl and Shou, Chong and Raha, Debasish and Gerstein, Mark and Snyder, Michael},
	Journal = {Science},
	Number = {5881},
	Pages = {1344--1349},
	Title = {The transcriptional landscape of the yeast genome defined by RNA sequencing.},
	Volume = {320},
	Year = {2008}}

@article{David:2006kl,
	Author = {David, Lior and Huber, Wolfgang and Granovskaia, Marina and Toedling, Joern and Palm, Curtis J and Bofkin, Lee and Jones, Ted and Davis, Ronald W and Steinmetz, Lars M},
	Journal = {Proc Natl Acad Sci U S A},
	Number = {14},
	Pages = {5320--5325},
	Title = {A high-resolution map of transcription in the yeast genome.},
	Volume = {103},
	Year = {2006}}


@article{Lee:2007cl,
	Author = {Lee, William and Tillo, Desiree and Bray, Nicolas and Morse, Randall H and Davis, Ronald W and Hughes, Timothy R and Nislow, Corey},
	Journal = {Nat Genet},
	Number = {10},
	Pages = {1235--1244},
	Title = {A high-resolution atlas of nucleosome occupancy in yeast.},
	Volume = {39},
	Year = {2007}}


@article{Siepel:2005ud,
	Author = {Siepel, Adam and Bejerano, Gill and Pedersen, Jakob S and Hinrichs, Angie S and Hou, Minmei and Rosenbloom, Kate and Clawson, Hiram and Spieth, John and Hillier, Ladeana W and Richards, Stephen and Weinstock, George M and Wilson, Richard K and Gibbs, Richard A and Kent, W James and Miller, Webb and Haussler, David},
	Journal = {Genome Res},
	Number = {8},
	Pages = {1034--1050},
	Title = {Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes.},
	Volume = {15},
	Year = {2005}}
}
\references{
  
}
\examples{
data(seqDataEx)
}
\keyword{datasets}
