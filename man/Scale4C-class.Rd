\name{Scale4C-class}

\docType{class}
\alias{Scale4C-class}
\alias{viewpoint}
\alias{viewpoint,Scale4C-method}
\alias{viewpoint<-}
\alias{viewpoint<-,Scale4C,numeric-method}
\alias{viewpointChromosome}
\alias{viewpointChromosome,Scale4C-method}
\alias{viewpointChromosome<-}
\alias{viewpointChromosome<-,Scale4C,character-method}
\alias{pointsOfInterest}
\alias{pointsOfInterest,Scale4C-method}
\alias{pointsOfInterest<-}
\alias{pointsOfInterest<-,Scale4C,GRanges-method}
\alias{rawData}
\alias{rawData,Scale4C-method}
\alias{rawData<-}
\alias{rawData<-,Scale4C,GRanges-method}
\alias{scaleSpace}
\alias{scaleSpace,Scale4C-method}
\alias{scaleSpace<-}
\alias{scaleSpace<-,Scale4C,SummarizedExperiment-method}
\alias{singularities}
\alias{singularities,Scale4C-method}
\alias{singularities<-}
\alias{singularities<-,Scale4C,GRanges-method}

\title{Class \code{"Scale4C"}}
\description{
This class is a container for information on a specific 4C-seq scale-space 
transformation. Stored information includes raw read data, the experiment's
viewpoint location (optional), possible points of interest, the scale-space
fingerprint map, and a list of identified singularities in scale-space.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Scale4C", ...)}.
}

\section{Slots}{
\describe{
\item{\code{viewpoint}:}{Object of class \code{"numeric"} representing the 
viewpoint's location}
\item{\code{viewpointChromosome}:}{Object of class \code{"character"} 
representing the viewpoint's chromosome}
\item{\code{pointsOfInterest}:}{Object of class \code{"GRanges"} 
representing any points of interest to be marked in the visualizations 
(usually near-cis based, i.e. close to the viewpoint)}
\item{\code{rawData}:}{Object of class \code{"GRanges"} representing 
the 4C-seq reads (or signal strength) of the experiment at given genomic 
positions}
\item{\code{scaleSpace}:}{Object of class \code{"SummarizedExperiment"} 
representing the gradually smoothed 4C-seq signal ('scale space') in its
first assay and the corresponding fingerprint map in its second assay.}
\item{\code{singularities}:}{Object of class \code{"GRanges"} 
representing singularities in the fingerprint map for the given 4C-seq signal}
}
}

\section{Methods}{
\describe{
\item{viewpoint<-}{\code{signature(object = "Scale4C", value = "numeric")}: 
Setter-method for the viewpoint slot.}
\item{viewpoint}{\code{signature(object = "Scale4C")}: 
Getter-method for the viewpoint slot.}
\item{viewpointChromosome<-}{\code{signature(object = "Scale4C", value = "character")}: 
Setter-method for the viewpointChromosome slot.}
\item{viewpointChromosome}{\code{signature(object = "Scale4C")}: 
Getter-method for the viewpointChromosome slot.}
\item{pointsOfInterest<-}{\code{signature(object = "Scale4C", 
    value = "GRanges")}: 
Setter-method for the pointsOfInterest slot.}
\item{pointsOfInterest}{\code{signature(object = "Scale4C")}: 
Getter-method for the pointsOfInterest slot.}
\item{rawData<-}{\code{signature(object = "Scale4C", value = "GRanges")}: 
Setter-method for the rawData slot.}
\item{rawData}{\code{signature(object = "Scale4C")}: 
Getter-method for the rawData slot.}
\item{scaleSpace<-}{\code{signature(object = "Scale4C", value = "matrix")}: 
Setter-method for the scaleSpace slot.}
\item{scaleSpace}{\code{signature(object = "Scale4C")}: 
Getter-method for the scaleSpace slot.}
\item{singularities<-}{\code{signature(object = "Scale4C", 
    value = "GRanges")}: 
Setter-method for the singularities slot.}
\item{singularities}{\code{signature(object = "Scale4C")}: 
Getter-method for the singularities slot.}
}
}

\author{
Carolin Walter
}

\examples{
    showClass("Scale4C")
}

\keyword{classes}
