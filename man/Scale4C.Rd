\name{Scale4C}

\alias{Scale4C}
\alias{Scale4C,numeric,character,GRanges-method}


\title{Creating a Scale4C object}

\description{
This function creates a Scale4C object. Data on the 4C-seq experiment, i.e. 
read counts per fragment and viewpoint coordinates, are stored and checked 
for plausibility.
}

\usage{Scale4C(viewpoint, viewpointChromosome, rawData)}

\arguments{
\item{viewpoint}{The experiment's viewpoint (start, single coordinate)}
\item{viewpointChromosome}{The experiment's viewpoint Chromosome}
\item{rawData}{Reads of the 4C-seq experiment per fragment on an interval 
of interest (GRanges object with position and read data)}
}

\details{
A \code{Scale4C} object contains the basic information on a 4C-seq experiment 
for a certain interval of interest, i.e. read counts at given positions. 
See \code{\link{Scale4C-class}} for details. Scale-space features such as 
fingerprint maps or tesselation are calculated during further steps of the
analysis by the appropriate functions.   

Scale4C expects the raw data to be in a simple data frame consisting of 
'position' and 'reads'. \code{importBasic4CseqData} allows to import fragment 
data from Basic4Cseq for convenience, however, preparing and importing a 
simple table with two columns into R is sufficient.
}

\value{
An instance of the Scale4C class.
}

\seealso{
    \code{\link{Scale4C-class}}
}

\examples{
    # create a Scale4C object from a Basic4Cseq export table with added 
    # viewpoint data 
    csvFile <- system.file("extdata", "liverData.csv", package="Scale4C")
    liverReads <- importBasic4CseqData(csvFile, viewpoint = 21160072, 
        viewpointChromosome = "chr10", distance = 1000000)
    liverData = Scale4C(rawData = liverReads, viewpoint = 21160072,
        viewpointChromosome = "chr10")
    liverData
}

\author{Carolin Walter}

\keyword{Scale4C}
