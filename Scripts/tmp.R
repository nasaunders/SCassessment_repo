#-*- coding: latin-1 -*-

### File: tmp.R
### Time-stamp: <2018-11-01 11:08:31 yreecht>
###
### Created: 20/06/2016	12:47:55
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

head(catchData)
summary(catchData$SpeciesID)
head(biologicalData)

## Identification of fields for which interpolations are calculated:
interpMetrics <- metrics[is.element(metrics[ , "metrics"],
                                    c("AbundanceDensity", "WeightDensity",
                                      "AbundanceDensityCorr", "WeightDensitCorr")), ]


## Interpolations:
## concaveHull(bufferedTracks=, SelectionMask=NULL, unsuitableArea=NULL, suitableArea=NULL, suitableAreaForce=FALSE,
## alpha=500, bufferWidth=20)

system.time(concaveHullList <- concaveHull.f(bufferedTracks = gBuffer(stationData,
                                                                      byid=TRUE,
                                                                      width=getOption("surveyIDW.bufferWidth"),
                                                                      capStyle="ROUND"),
                                             SelectionMask = mask,
                                             unsuitableArea = unsuitableGround,
                                             suitableArea = suitableGround,
                                             suitableAreaForce = getOption("surveyBedForce"),
                                             alpha = getOption("surveyIDW.alpha"),
                                             bufferWidth =  getOption("surveyIDW.bufferWidth")))

length(concaveHullList)

concaveHullList <- sapply(seq_len(length(concaveHullList)),
                          function(i, concaveHullList)
                   {
                       return(spChFIDs(concaveHullList[[i]],
                                       paste(letters[i],
                                             seq_len(length(concaveHullList[[i]])),
                                             sep = "")))
                   }, concaveHullList = concaveHullList, simplify = FALSE)

if (length(concaveHullList) > 1)
{
    concaveHull <- do.call(rbind, concaveHullList)
}else{
    concaveHull <- concaveHullList[[1]]
}

length(concaveHull)

plot(concaveHull[2, ])

test <- idw.byZone.alt(id=interpMetrics[6 , "field"],
                       data=stationData,
                       concaveHull=concaveHull,
                       alpha = getOption("surveyIDW.alpha"),
                       nmin = getOption("surveyIDW.nmin"),
                       nmax = getOption("surveyIDW.nmax"),
                       idp = getOption("surveyIDW.idp"),
                       cellsize = getOption("surveyIDW.cellsize"))



densInterp <- sapply(interpMetrics[ , "field"],
                     idw.byZone.alt,
                     data = stationData,
                     SelectionMask = mask,
                     unsuitableArea = unsuitableGround,
                     suitableArea = suitableGround,
                     suitableAreaForce = getOption("surveyBedForce"),
                     alpha = getOption("surveyIDW.alpha"),
                     nmin = getOption("surveyIDW.nmin"),
                     nmax = getOption("surveyIDW.nmax"),
                     idp = getOption("surveyIDW.idp"),
                     bufferWidth =  getOption("surveyIDW.bufferWidth"),
                     cellsize = getOption("surveyIDW.cellsize"),
                     simplify = FALSE)

class(test)

system.time()
width <- 6.5 * 5 / 6
X11(width = width,
    height = width * diff(bboxMapZone[2, ]) / diff(bboxMapZone[1, ]),
    pointsize = 11)

## Base map:
if (is.element("SpatialGridDataFrame",
               class(coastline)))
{
    image(coastline, red = 1, green = 2, blue = 3,
          xlim=bboxMapZone[1, ],
          ylim=bboxMapZone[2, ]## ,
          ## setParUsrBB = TRUE
          )
}else{
    plot(coastline,
         col = getOption("surveyPlot.colGround"),
         bg = getOption("surveyPlot.colSea"),
         xlim=bboxMapZone[1, ],
         ylim=bboxMapZone[2, ])
}

sapply(test, image, add = TRUE)

test2 <- Reduce(f = function(x, y)gUnion(x, y, byid = TRUE), x = test, )

length(test2)

plot(test2, add = TRUE, col = "green")



ls()[grepl("metric", ls())]
row.names(metrics[metrics$metrics %in% densUnits, ])

row.names(stationDate2)[! row.names(stationDate2) %in% row.names(stationDate)]

unique(row.names(stationDate2))
unique(row.names(stationDate))

head(stationData@data)

sapply(stationData@data, class)

levels(stationData@data$Station_ID)
nlevels(stationData@data$Station_ID)

all( ! grepl("/", stationData@data$TIMESTART))

traceback()

head(surveyMapZones@data)

bbox(stationData)

x <- sizeSelections[[1]]
head(biologicalData)

sapply(biologicalData, class)
unique(biologicalData$LD_HaulNo)

sum( ! is.na(biologicalData$Size) & ! is.na(biologicalData$Weight))

prettyQuantileBreaks(rnorm(100), n.categ=getOption("surveyIDW.n.categ"))

i <- names(densInterp)[2]

summary(values)
length(values)


 breaks <- quantile(values,
                    probs=seq(from=0, to=1, length.out=6 + 1))
duplicated(breaks)
n.digits <- ceiling(log10(1 / min(diff(breaks[ ! duplicated(breaks)])))) + 1

round(x=breaks, digits=n.digits)

prettyQuantileBreaks(values, n.categ=getOption("surveyIDW.n.categ"))

plot(mask, add = TRUE, col="red")

proj4string(mask)
par("usr")

dim(stationData)

test <- as.vector(gIntersects(stationData,
                              gBuffer(mask, width = 0),
                              byid = TRUE))

summary(test)

plot(stationData, add = TRUE, col="black")
n.digits <- 1
breaks <- c(0, 0, 4, 8, 12)

stationData@data[stationData$AbundanceDensity_Cerastoderma_edule > 0 &
                 stationData$AbundanceDensity_Cerastoderma_edule <= 1.1,]

dim(stationsTmp)
stationsTmp[343, ]@data

head(biologicalData)

head(catchData)

catchData[catchData$LD_HaulNo == "149b", ]
biologicalData[biologicalData$LD_HaulNo == "149b", ]

head(stationData@data)

head(stationData$Den_R)

testEff <- as.numeric(as.character(stationData$Den_R)) / as.numeric(as.character(stationData$Den_Q))
testEff[is.infinite(testEff)] <- NA
summary(testEff)

mean(testEff, na.rm = TRUE)

mean(as.numeric(as.character(stationData$Den_R)), nam.rm = TRUE) /
    mean(as.numeric(as.character(stationData$Den_Q)), nam.rm = TRUE)

testEff18 <- as.numeric(as.character(stationData$Den_R18)) / as.numeric(as.character(stationData$Den_Q18))
testEff18[is.infinite(testEff18)] <- NA
summary(testEff18)

mean(as.numeric(as.character(stationData$Den_R18)), nam.rm = TRUE) /
    mean(as.numeric(as.character(stationData$Den_Q18)), nam.rm = TRUE)

mean(testEff18, na.rm = TRUE)

testEff22 <- as.numeric(as.character(stationData$Den_R22)) / as.numeric(as.character(stationData$Den_Q22))
testEff22[is.infinite(testEff22)] <- NA
summary(testEff22)

mean(testEff22, na.rm = TRUE)

mean(as.numeric(as.character(stationData$Den_R22)), nam.rm = TRUE) /
    mean(as.numeric(as.character(stationData$Den_Q22)), nam.rm = TRUE)


head(subset(biologicalData, Size == 22))

v1 <- 105.41
v2 <- 45.41

n1 <- 18
n2 <- 18

## vTot <-

2 * sqrt(12.83)

1.67 * sqrt(12.83)

class(densInterpPolyFusion[[1]])

head(densInterpPolyFusion[[1]]@data)

par(mfrow = c(3, 4))
tmp <- biolCat[["AbundanceDensity_Ensis_arcuatus"]]
head(tmp)
by(tmp, tmp$stationID, function(x)hist(x$Size, main = unique(as.character(x$stationID))))

by(tmp, tmp$stationID,
   function(x)
{
    hist(x$Size, main = unique(as.character(x$stationID)))
    abline(v = mean(x$Size, na.rm = TRUE), col = "red", lwd = 2)
    abline(v = mean(tmp$Size, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    ## abline(v = mean(subset(tmp$Size,
    ##                        tmp$categDens == unique(as.character(x$categDens))),
    ##                 na.rm = TRUE), col = "blue", lwd = 2, lty = 2)
    mtext(unique(as.character(x$categDens)),
          adj = 0, col = "blue")
})

X11()
par(mfrow = c(3, 4))
by(tmp, tmp$stationID,
   function(x)
{
    hist(x$Weight, main = unique(as.character(x$stationID)))
    abline(v = mean(x$Weight, na.rm = TRUE), col = "red", lwd = 2)
    abline(v = mean(tmp$Weight, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
    ## abline(v = mean(subset(tmp$Size,
    ##                        tmp$categDens == unique(as.character(x$categDens))),
    ##                 na.rm = TRUE), col = "blue", lwd = 2, lty = 2)
    mtext(unique(as.character(x$categDens)),
          adj = 0, col = "blue")
})


################################################
length(densInterp[[3]][[1]])
head(densInterp[[3]][[1]]@data)

names(stationsCat)

stationsCat[[3]]@data[ , c("AbundanceDensity_Ensis_arcuatus", "categDens")]

nlevels(stationsCat[[3]]@data$categDens)

save(stationsCat, breaksList, file = "../../Other/Surveys/Data_DL16_BP_Razor.RData")

png(filename="../../Other/Surveys/Distri_density_Ensis_arcuatus.png",
    width=600*3, height=400*3, pointsize=12*3)
par(mfrow = c(2, 2),
    mar = c(3, 3, 3, 1) + 0.1,
    mgp = c(2, 0.7, 0), lwd = 2)

hist(stationsCat[[3]]@data[ , "AbundanceDensity_Ensis_arcuatus"],
     breaks = 10,
     main = "All stations",
     xlab = expression(Density~italic(Ensis~arcuatus)))
abline(v = tail(head(breaksList[[3]], -1), -1), col = "red", lwd = 2, lty = 2)
legend(x = "topright", inset = 0.05, legend = "Breaks", lwd = 2, lty = 2, col = "red")

by(data = stationsCat[[3]]@data,
   INDICES = stationsCat[[3]]@data$categDens,
   FUN = function(x)
{
    hist(x[ , "AbundanceDensity_Ensis_arcuatus"],
         main = paste("Stratum",
                      unique(as.character(x$categDens))),
         xlab = expression(Density~italic(Ensis~arcuatus)))
})
dev.off()

save(stationsCat, breaksList, file = "../../Other/Surveys/Data_DL16_BP_Razor.RData")

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
