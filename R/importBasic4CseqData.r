.importBasic4CseqData <- function(rawFile, viewpoint, viewpointChromosome, distance, useFragEnds = TRUE) {

    if (rawFile != "") {

        fullData = read.csv(rawFile, sep = "\t", stringsAsFactors = FALSE)
        fullData = subset(fullData,  fullData$fragmentStart > (viewpoint - distance) & fullData$fragmentEnd < (viewpoint + distance))
        if (useFragEnds) {
            fullDataLeft = subset(fullData, fullData$leftFragEndValid == TRUE)
            fullDataRight = subset(fullData, fullData$rightFragEndValid == TRUE)
            fullDataLeft$pos = fullDataLeft$fragmentStart
            fullDataLeft$reads = fullDataLeft$leftFragEndReads
            rawTableLeft = GRanges(viewpointChromosome, IRanges(as.numeric(fullDataLeft$pos), as.numeric(fullDataLeft$pos+1)), "reads" = fullDataLeft$reads)
            fullDataRight$pos = fullDataRight$fragmentEnd
            fullDataRight$reads = fullDataRight$rightFragEndReads
            rawTableRight = GRanges(viewpointChromosome, IRanges(as.numeric(fullDataRight$pos), as.numeric(fullDataRight$pos+1)), "reads" = fullDataRight$reads)
            rawTable = c(rawTableLeft, rawTableRight)
            rawTable$meanPosition = round((start(ranges(rawTable)) + end(ranges(rawTable)))/2)
        } else {
            fullData = subset(fullData, fullData$leftFragEndValid == TRUE & fullData$rightFragEndValid == TRUE)
            rawTable = GRanges(viewpointChromosome, IRanges(fullData$fragmentStart, fullData$fragmentEnd),
                "reads" = (fullData$leftFragEndReads + fullData$rightFragEndReads)/2, "meanPosition" = round((fullData$fragmentStart+fullData$fragmentEnd)/2))
        }
    } else {
        message("Please provide a valid file name for the raw data")
        rawTable = GRanges(character(0), IRanges(numeric(0), numeric(0)), "reads" = numeric(0), "meanPosition" = numeric(0))
    }

    rawTable = sort(rawTable)
    return(rawTable)
}


.addPointsOfInterest <- function(data, poi) {

    fullPositions = makeGRangesFromDataFrame(poi, keep.extra.columns = TRUE, ignore.strand = TRUE)
    tempRaw = rawData(data)
    tempRaw$index = 1:length(tempRaw)

    for (i in 1:length(fullPositions)) {
        tempPos = subset(tempRaw, tempRaw$meanPosition >= start(ranges(fullPositions))[i])
        if (length(tempPos) > 0) {
            fullPositions$index[i] = tempPos[1,]$index
        } else {
            fullPositions$index[i] = length(tempRaw)
        }
    }

    return(fullPositions)
}



setMethod("importBasic4CseqData",
    signature=signature(rawFile="character", viewpoint="numeric", viewpointChromosome="character", distance="numeric"),
    .importBasic4CseqData)


setMethod("addPointsOfInterest",
    signature=signature(data="Scale4C", poi="data.frame"),
    .addPointsOfInterest)
