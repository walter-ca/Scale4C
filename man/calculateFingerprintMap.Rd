\name{calculateFingerprintMap}

\alias{calculateFingerprintMap}
\alias{calculateFingerprintMap,Scale4C-method}

\title{Calculate scale-space fingerprint map for given 4C-seq data}

\description{
\code{calculateFingerprintMap} uses the scale space map to calculate the 
inflection points of the smoothed signals.
}

\usage{calculateFingerprintMap(data, maxSQSigma = 5000, epsilon = 0.0000001)}

\arguments{
\item{data}{Scale-space object for the 4C-seq data}
\item{maxSQSigma}{Maximum square sigma used to calculate the fingerprint map}
\item{epsilon}{Small numeric value (can also be zero); used to test for 
inflection points}
}

\details{
\pkg{Scale4C} uses Gauss kernels of increasing (square) sigma to smooth the 
original 4C-seq signal. The resulting inflection points for a chosen sigma are 
stored in the corresponding line of the fingerprint map, i.e. a 2D matrix 
(position x sigma). 
}

\value{A Scale4C object containing the fingerprint map for a Scale4C object, 
i.e. a matrix with quite a lot of zeros and the occasional 2 or -1 as symbols 
for inflection points. The fingerprint map is included as second assay of 
the Scale4C object's scaleSpace SummarizedExperiment slot.}

\examples{
    # read prepared example data
    data(liverData)
    # use small maxSQSigma for a fast example
    liverData<-calculateFingerprintMap(liverData, maxSQSigma = 50)
    head(t(assay(scaleSpace(liverData), 2))[1:10,1:20])
}

\author{Carolin Walter}

\keyword{calculateFingerprintMap}
