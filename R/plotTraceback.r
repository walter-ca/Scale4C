.plotTraceback <- function(data, maxSQSigma = -1, fileName = "tracebackPlot.pdf", width = 15, height = 15, useIndex = TRUE) {

    if (maxSQSigma == -1) {
        maxSQSigma = nrow(t(assay(scaleSpace(data), 2)))
    }
    coordPos = rawData(data)
    coordPos$index = 1:length(coordPos)
    pos = rawData(data)$meanPosition
    sl = singularities(data)

    if (useIndex) {
        pos = 1:length(pos)
    }

    if (fileName != "") {  
        pdf(file = fileName, title = "test", width = width, height = height, useDingbats = FALSE)
    }
    plot(c(min(pos),max(pos)), c(1,maxSQSigma), type = "n", xlab = "position", ylab = "sigma square")

    for (k in 1:maxSQSigma) {
        zcptemp = (t(assay(scaleSpace(data), 2)))[k,]
        indexPlus = (pos)[grepl("2", zcptemp)]
        indexMinus = (pos)[grepl("-1", zcptemp)]
        lines(indexPlus, rep(k, length(indexPlus)), type = "p", col = "grey80", cex = 0.3)
        lines(indexMinus, rep(k, length(indexMinus)), type = "p", col = "grey60", cex = 0.3)
    }

    for (i in 1:length(sl)) {
        lines(c(sl$left[i], start(ranges(sl))[i], end(ranges(sl))[i], sl$right[i]), c(1, sl$sqsigma[i] + 1, sl$sqsigma[i] + 1, 1), type = "l", col = "grey30", lwd = 1.5)
        if (abs(start(ranges(sl))[i] - end(ranges(sl))[i]) < 5) {
            lines(round((start(ranges(sl))[i] + end(ranges(sl))[i]) / 2), as.numeric(sl$sqsigma[i]) + 2, type = "p", col = "black", bg = "grey70", pch = 25, cex = 1.5)
        }
    }

    poi = pointsOfInterest(data)
    if (length(poi) != 0) {
        if (useIndex) {
            for (i in 1:length(poi)) {
                points(poi$index[i], round(maxSQSigma * 0.95), cex = 2.0, pch = 25, col = poi$colour[i], bg = poi$colour[i])
                text(poi$index[i], round(maxSQSigma * 0.99), poi$name[i], cex = 2.0)
                lines(rep(poi$index[i],2), c(0,round(maxSQSigma * 0.95)), type = "l", col = poi$colour[i], lwd = 2, lty = "dashed")
            }
        } else {
            for (i in 1:length(poi)) {
                points(start(ranges(poi))[i], round(maxSQSigma * 0.95), cex = 2.0, pch = 25, col = poi$colour[i], bg = poi$colour[i])
                text(start(ranges(poi))[i], round(maxSQSigma * 0.99), poi$name[i], cex = 2.0)
                lines(rep(start(ranges(poi))[i],2), c(0,round(maxSQSigma * 0.95)), type = "l", col = poi$colour[i], lwd = 2, lty = "dashed")
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

