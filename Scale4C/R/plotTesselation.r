.plotTesselation <- function(data, minSQSigma = 5, maxSQSigma = -1, maxVis = -1, fileName = "tesselationPlot.pdf", width = 5, height = 5, xInterval = 100, yInterval = 100, chosenColour = c("grey50", "moccasin", "lightskyblue1", "beige", "azure"), useIndex = TRUE) {

    sample = rawData(data)$reads
    positions = rawData(data)$position
    pos = (positions - min(positions)) / (max(positions) - min(positions)) * length(positions)

    if (useIndex) {
        pos = 1:length(positions)
        positions = 1:length(positions)
    }

    singData = singularities(data)
    singData$sqsigma = as.numeric(singData$sqsigma)
    singData$left = as.numeric(singData$left)
    singData$right = as.numeric(singData$right)

    if (maxSQSigma == -1) {
        maxSQSigma = max(singData$sqsigma)
    }
    if (maxVis == -1) {
        maxVis = maxSQSigma
    }

    singData = subset(singData, singData$sqsigma <= maxSQSigma & singData$sqsigma >= minSQSigma)

    singData[singData$sqsigma == maxSQSigma,]$sqsigma = 0.0001 * (1:nrow(singData[singData$sqsigma == maxSQSigma,])) + maxSQSigma

    singData$sqsigma = log2(singData$sqsigma)

    if (fileName != "") {  
        pdf(file = fileName, title = "test", width = width, height = height, useDingbats = FALSE)
    }
    plot(c(positions[1],positions[max(pos)]), c(-0.5,log2(maxVis)), type = "n", axes = FALSE, xlab = "fragments", ylab = "sigma square")

    singData = singData[order(singData$sqsigma),]

    for (k in nrow(singData):1) {

        largerSings = subset(singData, singData$sqsigma > singData$sqsigma[k])
        leftSings = subset(largerSings, (largerSings$right <= singData$left[k]))  
        rightSings = subset(largerSings, (largerSings$left >= singData$right[k]))
        topSings = subset(largerSings, (largerSings$left <= singData$left[k]) & (largerSings$right >= singData$right[k]))

        bottomSings = subset(singData, singData$sqsigma < singData$sqsigma[k] & singData$left >= singData$left[k] & singData$right <= singData$right[k])

        if (nrow(leftSings) == 0 & nrow(topSings) == 0) {
            tempLeft = positions[1]
        } else {
            tempLeft = max(leftSings$right, topSings$left)
        }
        if (nrow(rightSings) == 0 & nrow(topSings) == 0) {
            tempRight = positions[max(pos)]
        } else {
            tempRight = min(rightSings$left, topSings$right)
        }                 

        if (singData$type[k] == "peak" | singData$type[k] == "peak (VP)") {
            if (nrow(bottomSings) == 0) {
                bottom = -0.5
            } else {
                bottom = max(bottomSings$sqsigma)
            }
            rect(tempLeft, -0.5, tempRight, singData$sqsigma[k], col = chosenColour[5], border=chosenColour[1])
            rect(singData$left[k], bottom, singData$right[k], singData$sqsigma[k], col = chosenColour[2], border=chosenColour[1])
        }
        if (singData$type[k] == "valley") {
            if (nrow(bottomSings) == 0) {
                bottom = -0.5
            } else {
                bottom = max(bottomSings$sqsigma)
            }
            rect(tempLeft, -0.5, tempRight, singData$sqsigma[k], col = chosenColour[4], border=chosenColour[1])
            rect(singData$left[k], bottom, singData$right[k], singData$sqsigma[k], col = chosenColour[3], border=chosenColour[1])
        }
        if (singData$type[k] == "tracking problem") {
            if (nrow(topSings) %% 2 == 0) {
                tempCol = chosenColour[2]
            } else {
                tempCol = chosenColour[3]
            }
            if (nrow(bottomSings) == 0) {
                bottom = -0.5
            } else {
                bottom = max(bottomSings$sqsigma)
            }
            rect(positions[singData$left[k]], (bottom), positions[singData$right[k]], (singData$sqsigma[k]), col = tempCol, border=chosenColour[1])
        }
        lines(c(tempLeft, tempRight), rep(singData$sqsigma[k], 2), type = "l", col = chosenColour[1])
        lines(rep(singData$left[k],2), c(singData$sqsigma[k], -0.5), type = "l", col = chosenColour[1])
        lines(rep(singData$right[k],2), c(singData$sqsigma[k], -0.5), type = "l", col = chosenColour[1])
    }

    poi = pointsOfInterest(data)
    if (nrow(poi) != 0) {
        if (useIndex) {
            for (i in 1:nrow(poi)) {
                points(poi$index[i], (log2(maxVis * 0.75)), cex = 1.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(poi$index[i], (log2(maxVis * 0.85)), poi$name[i], cex = 1.0)
                lines(rep(poi$index[i],2), c(-0.5,(log2(maxVis * 0.75))), type = "l", col = as.vector(poi$colour[i]), lty = "dashed")
            }
        } else {
            for (i in 1:nrow(poi)) {
                points(poi$position[i], (log2(maxVis * 0.75)), cex = 1.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(poi$position[i], (log2(maxVis * 0.85)), poi$name[i], cex = 1.0)
                lines(rep(poi$position[i],2), c(-0.5,(log2(maxVis * 0.75))), type = "l", col = as.vector(poi$colour[i]), lty = "dashed")
            }
        }
    }

    axis(side = 1, at = c(1,positions[(0:(round(max(pos)/xInterval))*xInterval)], positions[max(pos)]))
    axis(side = 2, at = 0:log2((round(maxVis/yInterval)+1)*yInterval), labels = 2^(0:log2((round(maxVis/yInterval)+1)*yInterval)))

    if (fileName != "") {  
        dev.off()
    }

}



setMethod("plotTesselation",
    signature=signature(data="Scale4C"),
    .plotTesselation)

