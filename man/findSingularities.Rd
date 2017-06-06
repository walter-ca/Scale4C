\name{findSingularities}

\alias{findSingularities}
\alias{findSingularities,Scale4C-method}

\title{Find singularities in a provided scale-space fingerprint map}

\description{
This function allows to identify singular points in a scale-space fingerprint 
map. 
}

\usage{findSingularities(data, minSQSigma = 5, outputTrackingInfo = FALSE, 
guessViewpoint = FALSE, useIndex = TRUE)}

\arguments{
\item{data}{Scale-space object for the 4C-seq data}
\item{minSQSigma}{Minimum square sigma used to calculate singularities; for a 
square sigma of 1, the data can be quite chaotic and identified singularities 
are less prone to error when a minSQSigma of 2 or higher is used}
\item{outputTrackingInfo}{If TRUE, notify the user that a certain position / 
singularity causes problems during tracking in the fingerprint map}
\item{guessViewpoint}{If TRUE, add another 'peak' at the coordinates of the 
viewpoint, if provided. Extra singularities can also be added manually. The 
idea is to decrease running speed significantly by not actually calculating the 
largest singularity for a typical 4C-seq experiment, i.e. the main viewpoint peak. 
Its inflection point contours should easily be visible in the fingerprint map, 
provided that the viewpoint position is actually included in the raw data, but 
calculating the full contours requires a very high sigma that should usually not 
be needed to identify other singularities in the area. 
Cave: Viewpoint contours don't have to start directly next to the viewpoint 
coordinates.}
\item{useIndex}{If TRUE, use fragment index instead of genomic position data}
}

\details{
\code{findSingularities} identifies possible singular points in the fingerprint
map's contours, i.e. points where a line of '2' and '-1' in the matrix meet. 
Starting from those points in scale-space, the contours are traced back down. 
This 'localization step' ensures that the coordinates for a feature ('peak' or 
'valley') corresponding to a given singular point are as accurate as possible: 
Smoothing with a high-sigma Gauss kernel distorts the original signal somewhat, 
so that the inflection points identifying the start and the end of a certain 
feature 'move outwards'. 
}

\value{A data frame that lists the position where a singular point occurs 
(genomic position and scale-space sigma), plus the size of the feature as given 
by its minimal / left and maximal / right position.}

\examples{
    # read prepared example data
    data(liverData)
    singularities(liverData) = findSingularities(liverData, 5, useIndex = TRUE)
    singularities(liverData)
}

\author{Carolin Walter}

\keyword{findSingularities}
