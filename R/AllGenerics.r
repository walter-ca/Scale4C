
# generics for main calculations
setGeneric("calculateFingerprintMap", signature=c("data"),
    function(data, maxSQSigma = 5000, epsilon = 0.0000001)
        standardGeneric("calculateFingerprintMap"))

setGeneric("calculateScaleSpace", signature=c("data"),
    function(data, maxSQSigma = 5000)
        standardGeneric("calculateScaleSpace"))

setGeneric("findSingularities", signature=c("data"),
    function(data, minSQSigma = 5, outputTrackingInfo = FALSE, guessViewpoint = FALSE, useIndex = TRUE)
        standardGeneric("findSingularities"))

setGeneric("outputScaleSpaceTree", signature=c("data"),
    function(data, outputPeaks = TRUE, useLog = TRUE, useIndex = TRUE)
        standardGeneric("outputScaleSpaceTree"))

setGeneric("traceContour", signature=c("data", "singList"),
    function(data, singList, outputTrackingInfo = FALSE)
        standardGeneric("traceContour"))


# generics for export / import
setGeneric("importBasic4CseqData", signature=c("rawFile", "viewpoint", "viewpointChromosome", "distance"),
    function(rawFile, viewpoint, viewpointChromosome, distance, useFragEnds = TRUE)
        standardGeneric("importBasic4CseqData"))

setGeneric("addPointsOfInterest", signature=c("data", "poi"),
    function(data, poi)
        standardGeneric("addPointsOfInterest"))


# generics for visualization
setGeneric("plotInflectionPoints", signature=c("data", "sqsigma"),
    function(data, sqsigma, fileName = "inflectionPlot.pdf", width = 9, height = 5, maxVis = 5000, useIndex = TRUE, plotIP = TRUE)
        standardGeneric("plotInflectionPoints"))

setGeneric("plotTesselation", signature=c("data"),
    function(data, minSQSigma = 5, maxSQSigma = -1, maxVis = -1, fileName = "tesselationPlot.pdf", width = 5, height = 5, xInterval = 100, yInterval = 50, chosenColour = c("grey50", "moccasin", "lightskyblue1", "beige", "azure"), useIndex = TRUE)
        standardGeneric("plotTesselation"))

setGeneric("plotTraceback", signature=c("data"),
    function(data, maxSQSigma = -1, fileName = "tracebackPlot.pdf", width = 15, height = 15, useIndex = TRUE)
        standardGeneric("plotTraceback"))


# class generics
setGeneric("Scale4C", function(viewpoint, viewpointChromosome, rawData) {
    standardGeneric("Scale4C")})

setGeneric("viewpoint", function(object) {
    standardGeneric("viewpoint")})

setGeneric("viewpoint<-", function(object, value) {
    standardGeneric("viewpoint<-")})

setGeneric("viewpointChromosome", function(object) {
    standardGeneric("viewpointChromosome")})

setGeneric("viewpointChromosome<-", function(object, value) {
    standardGeneric("viewpointChromosome<-")})

setGeneric("pointsOfInterest", function(object) {
    standardGeneric("pointsOfInterest")})

setGeneric("pointsOfInterest<-", function(object, value) {
    standardGeneric("pointsOfInterest<-")})

setGeneric("rawData", function(object) {
    standardGeneric("rawData")})

setGeneric("rawData<-", function(object, value) {
    standardGeneric("rawData<-")})

setGeneric("scaleSpace", function(object) {
    standardGeneric("scaleSpace")})

setGeneric("scaleSpace<-", function(object, value) {
    standardGeneric("scaleSpace<-")})

setGeneric("singularities", function(object) {
    standardGeneric("singularities")})

setGeneric("singularities<-", function(object, value) {
    standardGeneric("singularities<-")})


