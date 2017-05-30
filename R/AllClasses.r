setClass("Scale4C",
    representation=list(           
        viewpoint="numeric",
        viewpointChromosome="character",
        pointsOfInterest="data.frame",
        rawData="data.frame",
        scaleSpace="matrix",
        fingerprint="matrix",
        singularities="data.frame"
    ),

    prototype=prototype(
        viewpoint=numeric(),
        viewpointChromosome=character(),
        pointsOfInterest=data.frame(),
        rawData=data.frame("position"=c(10,20,30),"reads"=c(1,2,3)),
        scaleSpace=matrix(),
        fingerprint=matrix(),
        singularities=data.frame()
    ),

    validity=function(object) {
        if ((ncol(rawData(object)) != 2) | !(is.numeric(rawData(object)$position))
            | !(is.numeric(rawData(object)$reads)))
            return("Please check the format of your fragment data: two numeric 
                columns (position, reads) are required")
        return(TRUE)
    }
)





