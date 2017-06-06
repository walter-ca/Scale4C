.plotInflectionPoints <- function(data, sqsigma, fileName = "inflectionPlot.pdf", width = 9, height = 5, maxVis = 5000, useIndex = TRUE, plotIP = TRUE) {

    if (plotIP) {
        ip = (t(assay(scaleSpace(data), 2)))[sqsigma-1,] 
    }
    pos = rawData(data)$meanPosition
    index = 1:length(pos)

    if (useIndex) {
        pos = 1:length(pos)
    }

    gaussKernel = matrix(0, 1, 1000)
    for (k in 1:length(gaussKernel)) {
        gaussKernel[1,k] = 1/sqrt(2*pi*(sqsigma)) * exp(1) ^ (-((k-ceiling(length(gaussKernel)/2))^2/(2*(sqsigma))))
    }
    smoothedData = kernel2dsmooth(as.matrix(t(rawData(data)$reads)), K = gaussKernel)

    if (plotIP) {
        indexPlusTemp = index[grepl("2", ip)]
        indexMinusTemp = index[grepl("-1", ip)]

        # remove doubled markers
        indexPlus = indexPlusTemp[1]
        jp = 1
        indexMinus = indexMinusTemp[1]
        jm = 1
        if (length(indexPlusTemp) > 1) {
            for (i in 1:(length(indexPlusTemp)-1)) {
                if (indexPlusTemp[i+1] != indexPlus[jp]+1) {
                    indexPlus[jp+1] = indexPlusTemp[i+1]
                    jp = jp+1
                }
            }
        }
        if (length(indexMinusTemp) > 1) {
            for (i in 1:(length(indexMinusTemp)-1)) {
                if (indexMinusTemp[i+1] != indexMinus[jm]+1) {
                    indexMinus[jm+1] = indexMinusTemp[i+1]
                    jm = jm+1
                }
            }
        }
        indexPlus = (pos)[indexPlus]
        indexMinus = (pos)[indexMinus]
    }

    if (fileName != "") {  
        pdf(file = fileName, title = "test", width = width, height = height, useDingbats = FALSE)
    }
    plot(c(min(pos), max(pos)), c(0, maxVis), type = "n", xlab = "position", ylab = "signal strength")
    lines(pos, smoothedData, type = "l")  
    if (plotIP) {
        lines(indexPlus, rep(maxVis, length(indexPlus)), type = "h", col = "grey80")
        lines(indexMinus, rep(maxVis, length(indexMinus)), type = "h", col = "grey60")
    }

    poi = pointsOfInterest(data)
    if (length(poi) != 0) {
        if (useIndex) {
            for (i in 1:length(poi)) {
                points(poi$index[i], round(maxVis * 0.95), cex = 1.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(poi$index[i], round(maxVis * 0.99), poi$name[i], cex = 1.0)
                lines(rep(poi$index[i],2), c(0,round(maxVis * 0.95)), type = "l", col = as.vector(poi$colour[i]), lty = "dashed")
            }
        } else {
            for (i in 1:length(poi)) {
                points(start(ranges(poi))[i], round(maxVis * 0.95), cex = 1.0, pch = 25, col = as.vector(poi$colour[i]), bg = as.vector(poi$colour[i]))
                text(start(ranges(poi))[i], round(maxVis * 0.99), poi$name[i], cex = 1.0)
                lines(rep(start(ranges(poi))[i],2), c(0,round(maxVis * 0.95)), type = "l", col = as.vector(poi$colour[i]), lty = "dashed")
            }
        }
    }

    if (fileName != "") {  
        dev.off()
    }

}



setMethod("plotInflectionPoints",
    signature=signature(data="Scale4C", sqsigma="numeric"),
    .plotInflectionPoints)

