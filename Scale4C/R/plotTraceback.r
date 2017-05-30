.plotTraceback <- function(data, maxSQSigma = -1, fileName = "tracebackPlot.pdf", width = 15, height = 15, useIndex = TRUE) {

    if (maxSQSigma == -1) {
        maxSQSigma = nrow(fingerprint(data))
    }
    coordPos = rawData(data)
    coordPos$index = 1:nrow(coordPos)
    pos = rawData(data)$position
    sl = singularities(data)

    if (useIndex) {
        pos = 1:length(pos)
    }

    if (fileName != "") {  
        pdf(file = fileName, title = "test", width = width, height = height, useDingbats = FALSE)
    }
    plot(c(min(pos),max(pos)), c(1,maxSQSigma), type = "n", xlab = "position", ylab = "sigma square")#, xaxt = "n")

    for (k in 1:maxSQSigma) {
        zcptemp = (fingerprint(data))[k,]
        indexPlus = (pos)[grepl("2", zcptemp)]
        indexMinus = (pos)[grepl("-1", zcptemp)]
        lines(indexPlus, rep(k, length(indexPlus)), type = "p", col = "grey80", cex = 0.3)
        lines(indexMinus, rep(k, length(indexMinus)), type = "p", col = "grey60", cex = 0.3)
    }

    for (i in 1:nrow(sl)) {
        lines(c(sl$left[i], sl$position1[i], sl$position2[i], sl$right[i]), c(1, as.numeric(sl$sqsigma[i]) + 1, as.numeric(sl$sqsigma[i]) + 1, 1), type = "l", col = "grey30", lwd = 1.5)
        if (abs(sl$position1[i] - sl$position2[i]) < 5) {
            lines(round((sl$position1[i] + sl$position2[i]) / 2), as.numeric(sl$sqsigma[i]) + 2, type = "p", col = "black", bg = "grey70", pch = 25, cex = 1.5)
        }
    }

    poi = pointsOfInterest(data)
    if (nrow(poi) != 0) {
        if (useIndex) {
            for (i in 1:nrow(poi)) {
                points(poi$index[i], round(maxSQSigma * 0.95), cex = 2.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(poi$index[i], round(maxSQSigma * 0.99), poi$name[i], cex = 2.0)
                lines(rep(poi$index[i],2), c(0,round(maxSQSigma * 0.95)), type = "l", col = as.vector(poi$colour[i]), lwd = 2, lty = "dashed")
            }
        } else {
            for (i in 1:nrow(poi)) {
                points(poi$position[i], round(maxSQSigma * 0.95), cex = 2.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(poi$position[i], round(maxSQSigma * 0.99), poi$name[i], cex = 2.0)
                lines(rep(poi$position[i],2), c(0,round(maxSQSigma * 0.95)), type = "l", col = as.vector(poi$colour[i]), lwd = 2, lty = "dashed")
            }
        }
    }

    if (fileName != "") {  
        dev.off()
    }

}



setMethod("plotTraceback",
    signature=signature(data="Scale4C"),
    .plotTraceback)

