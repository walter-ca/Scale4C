\name{plotTraceback}

\alias{plotTraceback}
\alias{plotTraceback,Scale4C-method}

\title{Draw the traceback results for a list of singular points on a fingerprint
map}

\description{
This method plots the traceback results together with fingerprint data, allowing
to check for possible errors during tracking. Problems during tracking can occur
if contours are very close, have holes, or if the singularity in question is not
recognized at all due to holes at the meeting point of both contours that form 
a singular point. Each singular point is marked with a grey triangle, and the 
traced left and right end of the corresponding feature are connected with grey 
lines. If a contour's end doesn't match the traceback line, manual correction is
possible in the singularity list.
}

\usage{plotTraceback(data, maxSQSigma = -1, fileName = "tracebackPlot.pdf", 
width = 15, height = 15, useIndex = TRUE)}

\arguments{
\item{data}{\code{Scale4C} object with singularity data}
\item{maxSQSigma}{Maximum square sigma (i.e. maximum y value) to be drawn; 
if -1 then all available rows in the fingerprint map are used}
\item{fileName}{Optional name for export file (pdf)}
\item{width}{Width of the plot}
\item{height}{Height of the plot}
\item{useIndex}{If TRUE, use fragment index for x-axis}
}

\value{A traceback plot, showing the traced singular points with their points of 
origin throughout different smoothing layers)}

\note{
PDF export is supported. If no plot file name is provided, the result is 
plotted on screen.
}

\examples{
    if(interactive()) {
        data(liverData)
        plotTraceback(liverData)
    }
}

\author{Carolin Walter}

\keyword{plotTraceback}
