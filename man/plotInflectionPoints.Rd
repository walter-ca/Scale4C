\name{plotInflectionPoints}

\alias{plotInflectionPoints}
\alias{plotInflectionPoints,Scale4C,numeric-method}

\title{Draw a smoothed near-cis profile with marked inflection points for 
a 4C-seq signal}

\description{
\code{plotInflectionPoints} plots the inflection points for a given square 
sigma (i.e. a row of the fingerprint map) onto a corresponding smoothed 
near-cis plot for the 4C-seq signal. This allows to check problematic parts 
of the fingerprint map in more detail (e.g. unclear tracking areas with 
close contours), and to improve possible corrections in the singularity list. 
Plotting the smoothed signal for a given square sigma before calculation of 
the fingerprint map is also possible.
}

\usage{plotInflectionPoints(data, sqsigma, fileName = "inflectionPlot.pdf", 
width = 9, height = 5, maxVis = 5000, useIndex = TRUE, plotIP = TRUE)}

\arguments{
\item{data}{\code{Scale4C} object with experimental 4C-seq data to be 
smoothed and processed}
\item{sqsigma}{Chosen square sigma, i.e. row of the fingerprint map to pick 
the inflection points from}
\item{fileName}{Optional name for export file (pdf)}
\item{width}{Width of the plot}
\item{height}{Height of the plot}
\item{maxVis}{Maximum y-axis value (read number, not sigma!) for 
visualization}
\item{useIndex}{If TRUE, use fragment index for x-axis}
\item{plotIP}{If TRUE, then mark chosen inflection points, if FALSE, simply 
plot smoothed data}
}

\value{A near-cis plot of the smoothed data with (optional) marked inflection 
points in darker or lighter grey, depending on their direction}

\note{
PDF export is supported. If no plot file name is provided, the result is 
plotted on screen.
}

\examples{
    data(liverData)
    plotInflectionPoints(liverData, 50)
}

\author{Carolin Walter}

\keyword{plotInflectionPoints}
