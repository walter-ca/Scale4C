.findSingularities <- function(data, minSQSigma = 5, outputTrackingInfo = FALSE, guessViewpoint = FALSE, useIndex = TRUE) {

    rawFP = (t(assay(scaleSpace(data), 2)))[(minSQSigma+1):ncol(scaleSpace(data)),]

    # track singularities
    singularities = list(0)

    for (i in 1:nrow(rawFP)) {

        tempRow = rawFP[i,]
        tempSing = 0
        for (j in 4:(ncol(rawFP)-4)) {
            if (tempRow[j] * tempRow[j-1] * tempRow[j-2] != 0  | tempRow[j] * tempRow[j-1] * tempRow[j-3] != 0
            | tempRow[j] * tempRow[j-2] * tempRow[j-3] != 0 | tempRow[j-1] * tempRow[j-2] * tempRow[j-3] != 0) {
                
                if (!is.null(tempSing)) {
                    last = tempSing[length(tempSing)]
                    if (j - last == 1) {
                        tempSing[length(tempSing)] = j
                    } else {
                        tempSing = c(tempSing, j)
                    }
                } else {
                    tempSing = c(tempSing, j)
                }
            }
        }
        singularities[[i]] = tempSing

        if (i > 1) {
            for (k in 1:length(singularities[[i]])) {
                if ( (singularities[[i]][k] %in% singularities[[i-1]]) | ((singularities[[i]][k] - 1) %in% singularities[[i-1]]) | ((singularities[[i]][k] + 1) %in% singularities[[i-1]]) ) {
                    singularities[[i-1]] = subset(singularities[[i-1]], !(singularities[[i-1]] %in% c(singularities[[i]][k], (singularities[[i]][k] - 1), (singularities[[i]][k] + 1) )))
                    if (length(singularities[[i-1]]) == 0) {
                        singularities[[i-1]] = 0
                    }
                }
            }
        }
    }

    # make list of singularities
    singList = GRanges(character(0), IRanges(integer(0), integer(0)), "sqsigma" = integer(0))

    for (i in 1:length(singularities)) {
        if (length(singularities[[i]]) > 1) {
            for (j in 1:length(singularities[[i]])) {
                singList = c(singList, GRanges(viewpointChromosome(data), IRanges(max(1, singularities[[i]][j] - 3), max(singularities[[i]][j]-2, 1)), "sqsigma" = i))
            }
        } else if (singularities[[i]] != 0) {
            singList = c(singList, GRanges(viewpointChromosome(data), IRanges(max(1, singularities[[i]][1] - 3), max(singularities[[i]]-2, 1)), "sqsigma" = i))
        }
    }

    if (guessViewpoint) {

        maxFP = rawFP[nrow(rawFP),]
        oneFound = FALSE

        for (i in 1:length(maxFP)) {
            
            if (maxFP[i] == "-1" & !(oneFound)) {
                oneFound = TRUE
                tempLeft = i
            }
            if (maxFP[i] == "2" & oneFound) {
                oneFound = FALSE
                tempRight = i
                singList = c(singList, GRanges(viewpointChromosome(data), IRanges(tempLeft, tempRight), "sqsigma" = nrow(rawFP))) 
            }
        }

    }

    singList = traceContour(data, singList, outputTrackingInfo)

    # clean up singularity list: a singularity with a smaller sqsigma that partly overlaps the interval of a singularity with a higher sqsigma is probably not a true singularity
    maxDetected = max(singList$sqsigma)
    fullSings = singList
    singList = subset(fullSings, as.numeric(fullSings$sqsigma) < maxDetected)
    maxSings = subset(fullSings, as.numeric(fullSings$sqsigma) == maxDetected)

    for (i in (length(singList):1)) {
        tempSing = singList[i,]
        tempLarger = subset(fullSings, as.numeric(fullSings$sqsigma) > as.numeric(tempSing$sqsigma))
        for (j in 1:length(tempLarger)) {
            tempLarge = tempLarger[j,]
            if ((as.numeric(tempSing$left) <= as.numeric(tempLarge$left) & as.numeric(tempSing$right) >= as.numeric(tempLarge$left)) 
            | (as.numeric(tempSing$left) <= as.numeric(tempLarge$right) & as.numeric(tempSing$right) >= as.numeric(tempLarge$right))) {
                singList$sqsigma[i] = -42
            }
        }    
    }
    singList = subset(singList, singList$sqsigma != -42)

    for (i in 1:(length(singList))) {
        tempSing = singList[i,]
        tempLarger = subset(singList, as.numeric(singList$sqsigma) > as.numeric(tempSing$sqsigma))
        if (length(tempLarger) > 0) {
            if (as.numeric(tempSing$left) %in% c(as.numeric(tempLarger$left), as.numeric(tempLarger$right), as.numeric(tempLarger$left) - 1, 
            as.numeric(tempLarger$left) + 1, as.numeric(tempLarger$right) - 1, as.numeric(tempLarger$right) + 1)
            | as.numeric(tempSing$right) %in% c(as.numeric(tempLarger$left), as.numeric(tempLarger$right), as.numeric(tempLarger$left) - 1, 
            as.numeric(tempLarger$left) + 1, as.numeric(tempLarger$right) - 1, as.numeric(tempLarger$right) + 1)) {
                singList$sqsigma[i] = -42
            }
        }
    }
    singList = subset(singList, singList$sqsigma != -42)
    singList = c(singList, maxSings)
    singList = sort(singList, by=~sqsigma)

    # convert position if necessary
    if (!(useIndex)) {
        ranges(singList) = IRanges(rawData(data)$meanPosition[start(ranges(singList))], rawData(data)$meanPosition[end(ranges(singList))])
        singList$left = rawData(data)$meanPosition[singList$left]
        singList$right = rawData(data)$meanPosition[singList$right]
    }

    singList$sqsigma = singList$sqsigma + minSQSigma - 1
    row.names(singList) = NULL

    return(singList)

}



setMethod("findSingularities",
    signature=signature(data="Scale4C"),
    .findSingularities)

