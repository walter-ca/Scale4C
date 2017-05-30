\name{liverDataVP}

\alias{liverDataVP}

\docType{data}

\title{
Example 4C-seq data set of fetal liver data, with added VP 
}

\description{
This data set contains an instance of a \code{Scale4C} object.

The 4C-seq data was taken from Stadhouders et al's fetal liver data set.
It contains a manually added viewpoint peak and a peak with a manual
coordinate correction.
}

\usage{
data("liverDataVP")
}

\format{
Formal class 'Scale4C'
}

\value{A pre-computed instance of a \code{Scale4C} object with fingerprint 
map and singularities}


\source{
Shortened version of Stadhouders et al's fetal liver data:

Stadhouders, R., Thongjuea, S., et al. (2012): Dynamic long-range chromatin
interactions control Myb proto-oncogene transcription during erythroid 
development. EMBO, 31, 986-999.
}

\examples{
    data("liverDataVP")
    liverDataVP
}

\keyword{datasets}
