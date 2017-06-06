# constructor for class "Scale4C"

setMethod("Scale4C", signature(viewpoint="numeric", viewpointChromosome="character", rawData="GRanges"),
    function(viewpoint, viewpointChromosome, rawData) {

        if (length(rawData) < 10) {
            message("The input raw data seems to be very short, please check if the data is correct")
        }

        if (viewpoint < 1000000) {
            message("The viewpoint seems to be close to the start of the chromosome, please check if the data is correct")
        }

        newScale4CObject = new("Scale4C", viewpoint = viewpoint, viewpointChromosome = viewpointChromosome, rawData = rawData)

        return(newScale4CObject)
    }
)


# "set" and "replace" methods for class "Scale4C"

setMethod("viewpoint", signature(object="Scale4C"),
    function(object) {
        return(object@viewpoint)
    }
)

setReplaceMethod("viewpoint",
    signature=signature(object="Scale4C", value="numeric"),
    function(object, value) {
        object@viewpoint = value
        return(object)
    }
)

setMethod("viewpointChromosome", signature(object="Scale4C"),
    function(object) {
        return(object@viewpointChromosome)
    }
)

setReplaceMethod("viewpointChromosome",
    signature=signature(object="Scale4C", value="character"),
    function(object, value) {
        object@viewpointChromosome = value
        return(object)
    }
)

setMethod("pointsOfInterest", signature(object="Scale4C"),
    function(object) {
        return(object@pointsOfInterest)
    }
)

setReplaceMethod("pointsOfInterest",
    signature=signature(object="Scale4C", value="GRanges"),
    function(object, value) {
        object@pointsOfInterest = value
        return(object)
    }
)

setMethod("rawData", signature(object="Scale4C"),
    function(object) {
        return(object@rawData)
    }
)

setReplaceMethod("rawData",
    signature=signature(object="Scale4C", value="GRanges"),
    function(object, value) {
        object@rawData = value
        return(object)
    }
)

setMethod("scaleSpace", signature(object="Scale4C"),
    function(object) {
        return(object@scaleSpace)
    }
)

setReplaceMethod("scaleSpace",
    signature=signature(object="Scale4C", value="SummarizedExperiment"),
    function(object, value) {
        object@scaleSpace = value
        return(object)
    }
)

setMethod("singularities", signature(object="Scale4C"),
    function(object) {
        return(object@singularities)
    }
)

setReplaceMethod("singularities",
    signature=signature(object="Scale4C", value="GRanges"),
    function(object, value) {
        fullSings = value
        fullSings = sort(fullSings, by=~sqsigma)
        for (i in (length(fullSings):1)) {
            tl = min(fullSings$left[i], fullSings$right[i])
            tr = max(fullSings$left[i], fullSings$right[i])
            fullSings$left[i] = tl
            fullSings$right[i] = tr
        }
        for (i in (length(fullSings):2)) {
            if (abs(fullSings[i]$sqsigma - fullSings[i-1]$sqsigma) < fullSings[i]$sqsigma/10 
                        && abs(start(ranges(fullSings))[i] - start(ranges(fullSings))[i-1]) < fullSings[i]$sqsigma/100) {
                    fullSings[i-1]$left = -42
            }
        }
        fullSings = subset(fullSings, fullSings$left != -42)
        maxDetected = max(fullSings$sqsigma)
        singList = subset(fullSings, fullSings$sqsigma < maxDetected)
        maxSings = subset(fullSings, fullSings$sqsigma == maxDetected)
        for (i in (length(singList):1)) {
            tempSing = singList[i,]
            tempLarger = subset(fullSings, fullSings$sqsigma > tempSing$sqsigma)
            for (j in 1:length(tempLarger)) {
                tempLarge = tempLarger[j,]
                if ((tempSing$left <= tempLarge$left && tempSing$right >= tempLarge$left) 
                || (tempSing$left <= tempLarge$right && tempSing$right >= tempLarge$right)) {
                    start(ranges(singList))[i] = -42
                }
            }    
        }
        outputSings = c(singList, maxSings)
        outputSings = sort(outputSings, by=~sqsigma)

        object@singularities = subset(outputSings, start(ranges(outputSings)) != -42)
        return(object)
    }
)



setMethod("show", "Scale4C",
    function(object){
        cat("4C-seq scale space data\n")
        cat("Type:", class(object), "\n")
        cat("Viewpoint: ", viewpointChromosome(object), ":", viewpoint(object), "\n")
        cat("Number of total fragments: ", length(rawData(object)), "\n")
        cat("Points of interest: ", length(pointsOfInterest(object)), "\n")
        cat("Maximum sigma of fingerprint map: ", ncol(scaleSpace(object)) - 1, "\n")
        cat("Number of singularities: ", length(singularities(object)), "\n")
    }
)
