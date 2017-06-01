setClass("Scale4C",
    representation=list(           
        viewpoint="numeric",
        viewpointChromosome="character",
        pointsOfInterest="GRanges",
        rawData="GRanges",
        scaleSpace="matrix",
        fingerprint="matrix",
        singularities="GRanges"
    ),

    prototype=prototype(
        viewpoint=numeric(),
        viewpointChromosome=character(),
        pointsOfInterest=GRanges(),
        rawData=GRanges(),
        scaleSpace=matrix(),
        fingerprint=matrix(),
        singularities=GRanges()
    ),

    validity=function(object) {
        if (ncol(mcols(rawData(object))) != 2)
            return("Please check the format of your fragment data")
        return(TRUE)
    }
)





