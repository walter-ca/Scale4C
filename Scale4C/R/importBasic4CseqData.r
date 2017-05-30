.importBasic4CseqData <- function(rawFile, viewpoint, distance, useFragEnds = TRUE) {

    if (rawFile != "") {

        fullData = read.csv(rawFile, sep = "\t", stringsAsFactors = FALSE)
        fullData = subset(fullData,  fullData$fragmentStart > (viewpoint - distance) & fullData$fragmentEnd < (viewpoint + distance))
        if (useFragEnds) {
            fullDataLeft = subset(fullData, fullData$leftFragEndValid == TRUE)
            fullDataRight = subset(fullData, fullData$rightFragEndValid == TRUE)
            fullDataLeft$pos = fullDataLeft$fragmentStart
            fullDataLeft$reads = fullDataLeft$leftFragEndReads
            rawTableLeft = data.frame("position" = fullDataLeft$pos, "reads" = fullDataLeft$reads)
            fullDataRight$pos = fullDataRight$fragmentEnd
            fullDataRight$reads = fullDataRight$rightFragEndReads
            rawTableRight = data.frame("position" = fullDataRight$pos, "reads" = fullDataRight$reads)
            rawTable = rbind(rawTableLeft, rawTableRight)
            rawTable = rawTable[order(rawTable$pos),]
        } else {
            fullData = subset(fullData, fullData$leftFragEndValid == TRUE & fullData$rightFragEndValid == TRUE)
            rawTable = data.frame("position" = (fullData$fragmentStart + fullData$fragmentEnd)/2, "reads" = (fullData$leftFragEndReads + fullData$rightFragEndReads)/2)
        }
    } else {
        message("Please provide a valid file name for the raw data")
        rawTable = data.frame("position" = numeric(), "reads" = numeric())
    }

    return(rawTable)
}


.importBedData <- function(rawFile, chromosome, viewpoint, distance) {

    if (rawFile != "") {

        fullData = read.csv(rawFile, sep = "\t", stringsAsFactors = FALSE)
        fullData = subset(fullData, fullData[,1] == chromosome)
        fullData = subset(fullData,  fullData[,2] > (viewpoint - distance) & fullData[,3] < (viewpoint + distance))
        rawTable = data.frame("position" = as.numeric(fullData[,2]), "reads" = as.numeric(fullData[,4]))
    } else {
        message("Please provide a valid file name for the raw data")
        rawTable = data.frame("position" = numeric(), "reads" = numeric())
    }

    return(rawTable)
}



.addPointsOfInterest <- function(data, poi) {

    fullPositions = data.frame("position" = poi$position, "index" = numeric(length(poi$position)), "colour" = poi$colour, "name" = poi$name)
    tempRaw = rawData(data)
    tempRaw$index = 1:nrow(tempRaw)

    for (i in 1:nrow(fullPositions)) {
        tempPos = subset(tempRaw, tempRaw$position >= fullPositions$position[i])
        if (nrow(tempPos) > 0) {
            tempIndex = tempPos[tempPos$position == min(tempPos$position),]$index
            fullPositions$index[i] = tempIndex
        } else {
            fullPositions$index[i] = nrow(tempRaw)
        }
    }

    return(fullPositions)
}



setMethod("importBasic4CseqData",
    signature=signature(rawFile="character", viewpoint="numeric", distance="numeric"),
    .importBasic4CseqData)

setMethod("importBedData",
    signature=signature(rawFile="character", chromosome="character", viewpoint="numeric", distance="numeric"),
    .importBedData)

setMethod("addPointsOfInterest",
    signature=signature(data="Scale4C", poi="data.frame"),
    .addPointsOfInterest)
