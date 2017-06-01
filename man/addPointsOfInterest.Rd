\name{addPointsOfInterest}

\alias{addPointsOfInterest}
\alias{addPointsOfInterest,Scale4C,data.frame-method}

\title{Add points of interest to a \code{Scale4C} object}

\description{
The function \code{addPointsOfInterest} adds marker points to a \code{Scale4C}
object, which are subsequently used to mark points of interest in created plots.
}

\usage{addPointsOfInterest(data, poi)}

\arguments{
\item{data}{Scale4C object the points are to be added to}
\item{poi}{Points of interest data, in a GRanges object. 
Important: column names must be specified and include 'colour' and
'name' for each point of interest with appropriate values}
}

\details{
The function \code{addPointsOfInterest} adds predefined points of interest to 
a \code{Scale4C} object. Each point of interest is defined by 'chr', 'start', 'end', 
'colour', and 'name'. A bed file or text file can be used to store the 
information, however, column names have to be added before import. Other 
additional columns are ignored by the function. The function then converts the
information to a GRanges object.
}

\value{A data frame that contains the data for all points of interest}

\examples{
    # import provided point of interest example, and check if import was 
    # successful
    data(liverData)
    poiFile <- system.file("extdata", "vp.txt", package="Scale4C")
    pointsOfInterest(liverData) <- addPointsOfInterest(liverData, 
        read.csv(poiFile, sep = "\t", stringsAsFactor = FALSE))
    head(pointsOfInterest(liverData))
}

\author{Carolin Walter}

\keyword{addPointsOfInterest}
