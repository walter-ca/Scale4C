\name{importBasic4CseqData}

\alias{importBasic4CseqData}
\alias{importBasic4CseqData,character,numeric,character,numeric-method}

\title{Import raw data from a provided Basic4Cseq output file}

\description{
A convenience function to easily include \pkg{Basic4Cseq} output data into 
\pkg{Scale4C}, \code{importBasic4CseqData} extracts valid fragments or valid
fragment ends from a typical \pkg{Basic4Cseq} output table.
}

\usage{importBasic4CseqData(rawFile, viewpoint, viewpointChromosome, 
    distance, useFragEnds = TRUE)}

\arguments{
\item{rawFile}{Name for the raw file}
\item{viewpoint}{Viewpoint position: only fragments around a certain point 
of interest are imported (doesn't have to be the actual viewpoint of the 
experiment, though)}
\item{viewpointChromosome}{Viewpoint chromosome of the experiment}
\item{distance}{Distance from the viewpoint: ony fragments within a certain 
distance of the viewpoint are imported}
\item{useFragEnds}{If TRUE, use full fragment end data, if FALSE, merge 
fragmentStart and fragmentEnd to a single item per fragment}
}

\details{
\code{importBasic4CseqData} is a convenience function to import data from 
\pkg{Basic4Cseq}. It can be ignored altogether if raw experimental data is 
imported from another source or with another function into R.
}

\value{A GRanges object that includes the experiment's raw data for 
further processing}

\examples{
    csvFile <- system.file("extdata", "liverData.csv", package="Scale4C")
    liverReads <- importBasic4CseqData(csvFile, viewpoint = 21160072, 
                viewpointChromosome = "chr10", distance = 1000000)
    head(liverReads)
}

\author{Carolin Walter}

\keyword{importBasic4CseqData}
