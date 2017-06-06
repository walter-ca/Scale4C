.calculateScaleSpace <- function(data, maxSQSigma = 5000) {

    scaleSpace = matrix(0, maxSQSigma+1, length(rawData(data)$reads))

    # gauss with different sigma^2
    for (j in 1:(maxSQSigma+1)) {

        gaussKernel = matrix(0, 1, 1000)
        if (j == 1) {
            for (k in 1:length(gaussKernel)) {
                gaussKernel[1,k] = 1/sqrt(2*pi*(0.5)) * exp(1) ^ (-((k-ceiling(length(gaussKernel)/2))^2/(2*(0.5))))
            }
        }
        else {
            for (k in 1:length(gaussKernel)) {
                gaussKernel[1,k] = 1/sqrt(2*pi*(j-1)) * exp(1) ^ (-((k-ceiling(length(gaussKernel)/2))^2/(2*(j-1))))
            }
        }
        scaleSpace[j,] = kernel2dsmooth(as.matrix(t(rawData(data)$reads)), K = gaussKernel)
    }

    fullScaleSpace = SummarizedExperiment(assays = list(scaleSpace = t(scaleSpace)), rowRanges = rawData(data), colData = data.frame("sigma" = c(0.5, 1:maxSQSigma)))

    return(fullScaleSpace)
}



setMethod("calculateScaleSpace",
    signature=signature(data="Scale4C"),
    .calculateScaleSpace)

