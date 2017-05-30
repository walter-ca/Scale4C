# constructor for class "Scale4C"

setMethod("Scale4C", signature(viewpoint="numeric", viewpointChromosome="character", rawData="data.frame"),
    function(viewpoint, viewpointChromosome, rawData) {

        if (nrow(rawData) < 10) {
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
    signature=signature(object="Scale4C", value="data.frame"),
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
    signature=signature(object="Scale4C", value="data.frame"),
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
    signature=signature(object="Scale4C", value="matrix"),
    function(object, value) {
        object@scaleSpace = value
        return(object)
    }
)

setMethod("fingerprint", signature(object="Scale4C"),
    function(object) {
        return(object@fingerprint)
    }
)

setReplaceMethod("fingerprint",
    signature=signature(object="Scale4C", value="matrix"),
    function(object, value) {
        object@fingerprint = value
        return(object)
    }
)

setMethod("singularities", signature(object="Scale4C"),
    function(object) {
        return(object@singularities)
    }
)

setReplaceMethod("singularities",
    signature=signature(object="Scale4C", value="data.frame"),
    function(object, value) {
        fullSings = value
        fullSings$position1 = as.numeric(fullSings$position1)
        fullSings$position2 = as.numeric(fullSings$position2)
        fullSings$sqsigma = as.numeric(fullSings$sqsigma)
        fullSings$left = as.numeric(fullSings$left)
        fullSings$right = as.numeric(fullSings$right)
        fullSings = fullSings[order(fullSings$sqsigma),]
        for (i in (nrow(fullSings):1)) {
            tl = min(fullSings$left[i], fullSings$right[i])
            tr = max(fullSings$left[i], fullSings$right[i])
            fullSings$left[i] = tl
            fullSings$right[i] = tr
        }
        for (i in (nrow(fullSings):2)) {
            if (abs(fullSings[i,]$sqsigma - fullSings[i-1,]$sqsigma) < fullSings[i,]$sqsigma/10 
                        && abs(fullSings[i,]$position1 - fullSings[i-1,]$position1) < fullSings[i,]$sqsigma/100) {
                    fullSings[i-1,4] = -42
            }
        }
        fullSings = subset(fullSings, fullSings[,4] != -42)
        maxDetected = max(fullSings$sqsigma)
        singList = subset(fullSings, fullSings$sqsigma < maxDetected)
        maxSings = subset(fullSings, fullSings$sqsigma == maxDetected)
        for (i in (nrow(singList):1)) {
            tempSing = singList[i,]
            tempLarger = subset(fullSings, as.numeric(fullSings$sqsigma) > as.numeric(tempSing$sqsigma))
            for (j in 1:nrow(tempLarger)) {
                tempLarge = tempLarger[j,]
                if ((as.numeric(tempSing$left) <= as.numeric(tempLarge$left) && as.numeric(tempSing$right) >= as.numeric(tempLarge$left)) 
                || (as.numeric(tempSing$left) <= as.numeric(tempLarge$right) && as.numeric(tempSing$right) >= as.numeric(tempLarge$right))) {
                    singList[i,3] = -42
                }
            }    
        }
        outputSings = rbind(singList, maxSings)
        outputSings = outputSings[order(outputSings$sqsigma),]

        object@singularities = subset(outputSings, outputSings[,3] != -42)
        return(object)
    }
)



setMethod("show", "Scale4C",
    function(object){
        cat("4C-seq scale space data\n")
        cat("Type:", class(object), "\n")
        cat("Viewpoint: ", viewpointChromosome(object), ":", viewpoint(object), "\n")
        cat("Number of total fragments: ", nrow(rawData(object)), "\n")
        cat("Points of interest: ", nrow(pointsOfInterest(object)), "\n")
        cat("Maximum sigma of fingerprint map: ", nrow(fingerprint(object)) - 1, "\n")
        cat("Number of singularities: ", nrow(singularities(object)), "\n")
    }
)
