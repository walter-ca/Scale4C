\name{outputScaleSpaceTree}

\alias{outputScaleSpaceTree}
\alias{outputScaleSpaceTree,Scale4C-method}

\title{Output list of all features for a given scale-space map}

\description{
This function provides a list of features for a given fingerprint map in 
scale-space, with position and range of sigma for which the feature in 
question exists 
}

\usage{outputScaleSpaceTree(data, outputPeaks = TRUE, useLog = TRUE, 
    useIndex = TRUE)}

\arguments{
\item{data}{a \code{Scale4C} object with singularity data}
\item{outputPeaks}{If TRUE, output GRanges peak list only, if FALSE, 
    also output valley data in a larger table}
\item{useLog}{If TRUE, use a log2 transformation on the square sigma values 
(fewer changes and fewer singularities for high sigma, in contrast to low 
sigma)}
\item{useIndex}{If TRUE, use fragment position}
}

\details{
Similar to \code{plotTesselation}, \code{outputScaleSpaceTree} analyzes a 
list of singular points and calculates corresponding features, i.e. 'peaks' 
and 'valleys'. Each singular point marks the disappearance (or occurrence, 
depending on the view) of a feature in scale space: With increasing square 
sigma as smoothing parameter for the Gauss kernel, smaller features are 
merged into larger features. In case of Gauss smoothing, one feature is 
always surrounded by two features of the opposite type, e.g. a 'peak' is 
surrounded by two 'valleys'. If a 'peak' is smoothed out, it is replaced by 
a new valley formed of the former peak's adjacent valleys. The singularity 
list contains only direct information on those 'central' features; 
\code{outputScaleSpaceTree} adds data on the direct neighbours / adjacent 
features and also provides the sigma ranges for the features as a measure 
of their stability throughout the smoothing process. Mean read counts for
the identitied features are also provided ("signal").
If \code{outputPeaks} is true, a reduced list of peaks is printed, while
omitting valleys or the central-left-right structural information.
}

\value{A GRanges object that includes all features as identified through singular 
points, plus 'neighbour features' at each side (each 'peak' is surrounded by 
two 'valleys' and vice versa for Gauss kernel smoothing), with positions and 
range of sigma for which the feature in question remains stable}

\examples{
    # read prepared example data
    data(liverDataVP)
    output = outputScaleSpaceTree(liverDataVP, useLog = FALSE)
    head(output) 
}

\author{Carolin Walter}

\keyword{outputScaleSpaceTree}
