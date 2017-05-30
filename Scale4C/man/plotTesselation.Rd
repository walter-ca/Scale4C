\name{plotTesselation}

\alias{plotTesselation}
\alias{plotTesselation,Scale4C-method}

\title{Draw the final scale space tesselation}

\description{
This method draws the final scale space tesselation, as specified by the list 
of singularities identified for a \code{Scale4C} object. Features are marked 
with different colours; for the default colour scheme, brown corresponds to 
'peaks' and blue to 'valleys', while slightly darker colours mark features 
originating from singularities ('central' features in a set of three features, 
e.g. 'valley-peak-valley' or 'peak-valley-peak') and lighter colours the two 
adjacent features. Different colours for 'central' and 'adjacent' features allow
for optical quality control of the tesselation: a 'central' / dark feature's 
direct predecessor or successor (y-axis) can't be of the same colour (i.e. a 
'peak' that passes through a singularity is smoothed out into a 'valley'), and 
neighbouring intervals have to be of the opposing (but lighter) colour (i.e. 
each 'peak' is surrounded by two 'valleys' for Gauss kernel smoothing). The 
same is not necessarily true for an 'adjacent' / light feature, however.
}

\usage{plotTesselation(data, minSQSigma = 5, maxSQSigma = -1, maxVis = -1, 
fileName = "tesselationPlot.pdf", width = 5, height = 5, xInterval = 100, 
yInterval = 50, chosenColour = c("grey50", "moccasin", "lightskyblue1", 
"beige", "azure"), useIndex = TRUE)}

\arguments{
\item{data}{\code{Scale4C} object with singularity data}
\item{minSQSigma}{Minimum square sigma to consider}
\item{maxSQSigma}{Maximum square sigma to consider; if -1 then the number of 
rows in the fingerprint map is used}
\item{maxVis}{Maximum y value for visualization (doesn't have to be maxSQSigma);
if -1 also defaults to number of rows in the fingerprint map}
\item{fileName}{Optional name for export file (pdf)}
\item{width}{Width of the plot}
\item{height}{Height of the plot}
\item{xInterval}{Interval length for x-axis}
\item{yInterval}{Interval length for y-axis}
\item{chosenColour}{Chosen colours for the tesselation plot, five in total. 
Colour 1 is used for the actual lines of the plot, colour 2 for 'central 
peaks', colour 3 for 'central valleys', colour 4 for 'adjacent peaks', and 
colour 5 for 'adjacent valleys'}
\item{useIndex}{If TRUE, use fragment index for x-axis}
}

\value{A tesselation plot, showing different features of the scale space with 
their range of existance (square sigma) and position)}

\note{
PDF export is supported. If no plot file name is provided, the result is 
plotted on screen.
}

\examples{
    if(interactive()) {
        data(liverData)
        plotTesselation(liverData)
    }
}

\author{Carolin Walter}

\keyword{plotTesselation}
