\name{importBedData}

\alias{importBedData}
\alias{importBedData,character,character,numeric,numeric-method}

\title{Import raw data from a provided bed file}

\description{
A convenience function to easily include bed data into \pkg{Scale4C}, 
\code{importBedData} imports a bed-like file into R. 
}

\usage{importBedData(rawFile, chromosome, viewpoint, distance)}

\arguments{
\item{rawFile}{Name for the raw file}
\item{chromosome}{Name of the viewpoint's chromosome}
\item{viewpoint}{Viewpoint position: only fragments around a certain point 
of interest are imported (doesn't have to be the actual viewpoint of the 
experiment, though)}
\item{distance}{Distance from the viewpoint: ony fragments within a certain 
distance of the viewpoint are imported}
}

\details{
\code{importBedData} is a convenience function to import bed-like data. 
It is assumed that column 1 contains information on the chromosome, columns 
2 and 3 refer to the start and end of the intervals in question, and column 
4 has data on the interval (or fragment) raw read number.
It can be ignored altogether if raw experimental data is imported from 
another source or with another function into R.
}

\value{A data frame with two columns, consisting of 'position' and 'reads', 
that includes the experiment's raw data for further processing}

\examples{
    csvFile <- system.file("extdata", "demoData.bed", package="Scale4C")
    liverReads <- importBedData(csvFile, chromosome = "chr10", 
                    viewpoint = 21160072, distance = 1000000)
    head(liverReads)
}

\author{Carolin Walter}

\keyword{importBedData}
