.traceContour <- function(data, singList, outputTrackingInfo = FALSE) {

    rawFP = t(assay(scaleSpace(data), 2))
    singListType = NULL

    for (i in 1:length(singList)) {

        sigma = singList$sqsigma[i]

        posLeft = max(1, start(ranges(singList))[i])
        valueLeft = as.numeric(rawFP[sigma, posLeft])

        if (valueLeft == 0) {
            posLeft = max(1, start(ranges(singList))[i] - 1)
            valueLeft = as.numeric(rawFP[sigma, posLeft])
        }
        posRight = max(1, end(ranges(singList))[i])
        valueRight = as.numeric(rawFP[sigma, posRight])
        if (valueRight == 0) {
            posRight = end(ranges(singList))[i] + 1
            valueRight = as.numeric(rawFP[sigma, posRight])
        }

        if (valueLeft*valueRight == -2) {

            if (valueLeft == -1 & valueRight == 2) {
                singListType[i] = "peak"
            } else if (valueLeft == 2 & valueRight == -1) {
                singListType[i] = "valley"
            }

            if (sigma != 1) {
                for (j in (sigma-1):1) {
                    newValueLeft = as.numeric(rawFP[j, posLeft])
                    newPosLeft = posLeft

                    if (newValueLeft != valueLeft) {
                        newValueLeft = as.numeric(rawFP[j, posLeft+1])
                        newPosLeft = posLeft+1
                        if (newValueLeft != valueLeft) {
                            newValueLeft = as.numeric(rawFP[j, posLeft-1])
                            newPosLeft = max(1, posLeft-1)
                            if (newValueLeft != valueLeft) {
                                if (outputTrackingInfo) {
                                    print(paste("problems while tracking ", posLeft, " ", sigma, sep = ""))
                                    newValueLeft = valueLeft
                                }
                            }
                        }
                    }
                    valueLeft = newValueLeft
                    posLeft = newPosLeft
                    singList$left[i] = posLeft

                    newValueRight = as.numeric(rawFP[j, posRight])
                    newPosRight = posRight

                    if (newValueRight != valueRight) {
                        newValueRight = as.numeric(rawFP[j, posRight-1])
                        newPosRight = max(1, posRight-1)
                        if (newValueRight != valueRight) {
                            newValueRight = as.numeric(rawFP[j, posRight+1])
                            newPosRight = posRight+1
                            if (newValueRight != valueRight) {
                                if (outputTrackingInfo) {
                                    print(paste("problems while tracking ", posRight, " ", sigma, sep = ""))
                                    newValueRight = valueRight
                                }
                            }
                        }
                    }
                    valueRight = newValueRight
                    posRight = newPosRight
                    singList$right[i] = posRight
                }
            } else {
                singList$left[i] = posLeft
                singList$right[i] = posRight
            }

        } else {
            if (outputTrackingInfo) {
                print(paste("problems while tracking ", posLeft, " ", sigma, sep = ""))
            }
            singList$left[i] = posLeft
            singList$right[i] = posRight
            singListType[i] = "tracking problem"
        }   
    }

    singList$type = singListType

    return(singList)
}


setMethod("traceContour",
    signature=signature(data="Scale4C", singList = "GRanges"),
    .traceContour)


