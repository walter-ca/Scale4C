.outputScaleSpaceTree <- function(data, outputPeaks = TRUE, useLog = TRUE, useIndex = TRUE) {

    singData = singularities(data)
    position = rawData(data)
    position$index = 1:nrow(position)
    colnames(position) = c("pos","reads","index")
    scaleSpace = scaleSpace(data)
    if (useIndex) {
        position$pos = position$index
    }

    singData$sqsigma = as.numeric(singData$sqsigma)
    singData$left = as.numeric(singData$left)
    singData$right = as.numeric(singData$right)
    singData = singData[order(singData$sqsigma),]

    if (useLog) {
        singData$sqsigma = log2(as.numeric(singData$sqsigma))
    }

    output = data.frame("centre_maxSSQ" = numeric(), "centre_minSSQ" = numeric(), "centre_leftPos" = numeric(), 
        "centre_rightPos" = numeric(), "centre_type" = character(), "centre_length" = numeric(), "centre_signal" = numeric(),
        "left_maxSSQ" = numeric(), "left_minSSQ" = numeric(), "left_leftPos" = numeric(), "left_rightPos" = numeric(), "left_type" = character(), 
        "left_length" = numeric(), "left_signal" = numeric(),
        "right_maxSSQ" = numeric(), "right_minSSQ" = numeric(), "right_leftPos" = numeric(), "right_rightPos" = numeric(), "right_type" = character(),
        "right_length" = numeric(), "right_signal" = numeric(),
        stringsAsFactors = FALSE)

    for (k in nrow(singData):1) {

        largerSings = subset(singData, singData$sqsigma > singData$sqsigma[k])
        leftSings = subset(largerSings, (largerSings$right <= singData$left[k]))  
        rightSings = subset(largerSings, (largerSings$left >= singData$right[k]))
        topSings = subset(largerSings, (largerSings$left <= singData$left[k]) & (largerSings$right >= singData$right[k]))
        bottomSings = subset(singData, singData$sqsigma < singData$sqsigma[k] & singData$left >= singData$left[k] & singData$right <= singData$right[k])

        if (nrow(leftSings) == 0 & nrow(topSings) == 0) {
            tempLeft = min(position$pos)
        } else {
            tempLeft = max(leftSings$right, topSings$left)
        }
        if (nrow(rightSings) == 0 & nrow(topSings) == 0) {
            tempRight = max(position$pos)
        } else {
            tempRight = min(rightSings$left, topSings$right)
        }
        if (nrow(bottomSings) == 0) {
            bottom = 0
        } else {
            bottom = max(bottomSings$sqsigma)
        }
        if (singData$type[k] == "peak" | singData$type[k] == "peak (VP)") {
            invType = "valley"
        } else if (singData$type[k] == "valley") {
            invType = "peak"
        } else {
            invType = "tracking problem"
        }
        bottomLeftSings = subset(singData, singData$sqsigma < singData$sqsigma[k] & singData$left >= tempLeft & singData$right <= singData$left[k])
        bottomRightSings = subset(singData, singData$sqsigma < singData$sqsigma[k] & singData$left >= singData$right[k] & singData$right <= tempRight)
        if (nrow(bottomLeftSings) == 0) {
            bottomLeft = 0
        } else {
            bottomLeft = max(bottomLeftSings$sqsigma)
        }
        if (nrow(bottomRightSings) == 0) {
            bottomRight = 0
        } else {
            bottomRight = max(bottomRightSings$sqsigma)
        }

        cl = min(singData$left[k], singData$right[k])
        cr = max(singData$left[k], singData$right[k])
        centreData = subset(position, position$pos >= cl & position$pos < cr)
        if (nrow(centreData) > 0) {
            centreLength = nrow(centreData)
            centreSignal = mean(position$reads[centreData[1,3]:centreData[nrow(centreData),3]])
        } else {
            centreLength = -1
            centreSignal = -1
        }
        ll = min(singData$left[k], tempLeft)
        lr = max(singData$left[k], tempLeft)
        leftData = subset(position, position$pos >= ll & position$pos < lr)
        if (nrow(leftData) > 0) {
            leftLength = nrow(leftData)
            leftSignal = mean(position$reads[leftData[1,3]:leftData[nrow(leftData),3]])
        } else {
            leftLength = -1
            leftSignal = -1
        }
        rl = min(singData$right[k], tempRight)
        rr = max(singData$right[k], tempRight)
        rightData = subset(position, position$pos >= rl & position$pos < rr)
        if (nrow(rightData) > 0) {
            rightLength = nrow(rightData)
            rightSignal = mean(position$reads[rightData[1,3]:rightData[nrow(rightData),3]])
        } else {
            rightLength = -1
            rightSignal = -1
        }
        output[k,] = c(singData$sqsigma[k], bottom, singData$left[k], singData$right[k], singData$type[k], centreLength, round(centreSignal,2),
            singData$sqsigma[k], bottomLeft, ll, lr, invType, leftLength, round(leftSignal,2),
            singData$sqsigma[k], bottomRight, rl, rr, invType, rightLength, round(rightSignal,2))

    }

    if (outputPeaks) {
        fullTable = output
        output = fullTable[,1:7]
        temp1 = fullTable[,8:14]
        temp2 = fullTable[,15:21]
        colnames(output) = gsub("centre_", "", colnames(output))
        colnames(temp1) = gsub("left_", "", colnames(output))
        colnames(temp2) = gsub("right_", "", colnames(output))
        output = rbind(output, temp1)
        output = rbind(output, temp2)
        output = subset(output, output[,5] != "valley")
        output[,1] = viewpointChromosome(data)
        output[,6] = round(as.numeric(output[,6]) * as.numeric(output[,7]))
        output = output[,c(1,3,4,6)]
        colnames(output) = c("chr", "start", "end", "total")
    }

    return(output)
}



setMethod("outputScaleSpaceTree",
    signature=signature(data="Scale4C"),
    .outputScaleSpaceTree)

