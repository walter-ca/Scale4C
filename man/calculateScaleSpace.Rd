\name{calculateScaleSpace}

\alias{calculateScaleSpace}
\alias{calculateScaleSpace,Scale4C-method}

\title{Calculate scale-space representation for given 4C-seq data}

\description{
\pkg{Scale4C} uses Gauss kernels of increasing (square) sigma to smooth the 
original 4C-seq signal. The resulting data is stored in a 2D matrix 
(position x sigma). 
}

\usage{calculateScaleSpace(data, maxSQSigma = 5000)}

\arguments{
\item{data}{Scale-space object for the 4C-seq data}
\item{maxSQSigma}{Maximum square sigma used to calculate the scale space}
}

\details{
The central idea of the scale-space transformation is to smooth the original 
signal with increasing strength, identify inflection points, track those 
inflection points throughout the different smoothing layers, and find 
singularities in those inflection point 'lines'. In case of 4C-seq data, this 
corresponds to smoothing the signal gradually, while making notes when 
features such as 'peaks' or 'valleys' disappear by merging with other 
features. \code{calculateScaleSpace} smoothes the original signal up to 
a provided smoothing factor square sigma (Gauss kernel).
}

\value{A \code{SummarizedExperiment} that contains in its first assay the scale 
space representation for a \code{Scale4C} object}

\examples{
    # read prepared example data
    data(liverData)
    # use small maxSQSigma for a fast example
    scaleSpace(liverData)<-calculateScaleSpace(liverData, maxSQSigma = 10)
    head(t(assay(scaleSpace(liverData), 1))[,1:20])
}

\author{Carolin Walter}

\keyword{calculateScaleSpace}
