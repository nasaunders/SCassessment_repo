#-*- coding: latin-1 -*-

### File: 1_Summary_spatial_data.R
### Time-stamp: <2017-06-01 11:04:39 yreecht>
###
### Created: 15/03/2016	12:02:05
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Create a summary of station data useful for filling the inshore upload template.
####################################################################################################


## Finding information on the projection system (may take a while):
infoProj <- projEPSGinfo(crs = proj4string(stationData), SP = stationData[1:3, ])

## Effort calculation (class specific):
if (grepl("^SpatialLines", class(stationData)))
{
    stationData$length <- gLength(stationData, byid = TRUE)
    stationData$effort <- stationData$length * gearWidth
}else{
    if (grepl("^SpatialPoints", class(stationData)))
    {
        if (length(gearWidth) > 1 &&
            ! is.null(names(gearWidth)))
        {
            message("\n## Multi-gears point data: effort calculated later!")
        }else{
            stationData$effort <- gearWidth[1] ^ 2
        }
    }else{
        warning("Effort calculation not implemented for this class of object.")
    }
}

## stationData$density <- stationData$NBOYSTERS / stationData$effort

## Formatting time:
if (all( ! grepl("/", stationData@data$TIMESTART)) &&
    any(grepl("-", stationData@data$TIMESTART)))
{
    stationData$TIMESTART <- as.POSIXct(as.character(stationData$TIMESTART),
                                        format = "%Y-%m-%d %H:%M:%S")
}else{
    stationData$TIMESTART <- as.POSIXct(as.character(stationData$TIMESTART),
                                        format = "%d/%m/%Y %H:%M:%S")
}

stationData$startDate <- format(stationData$TIMESTART, format = "%d/%m/%Y")
stationData$startTime <- format(stationData$TIMESTART, format = "%H:%M:%S")

stationData$TRACKID <- as.character(stationData$TRACKID)

## Barycentres of tracks:
coordsTmp <- gCentroid(stationData, byid = TRUE)

## Binding and ordering information:
stationsSummary <- stationData
stationsSummary@data <- cbind(stationsSummary@data,
                              coordsTmp@coords,
                              EPSG = rep(infoProj[ , c("code")],
                                         times = nrow(stationsSummary)),
                              projName = rep(infoProj[ , c("note")],
                                             times = nrow(stationsSummary)))


stationsSummary <- stationsSummary[order(stationsSummary$TRACKID), ]

write.csv(stationsSummary@data,
          file = file.path(ResultsPath, "stationInfo.csv"))


message("You can use the output file \n(",
        file.path(getwd(), ResultsPath, "stationInfo.csv"),
        ") \nas a support for filling the inshore database upload template.\n",
        "Then export the \"CatchData\" and \"BiologicalData\" tabs as .csv file\n",
        "(those specified as *catchDataFile* and *biologicalDataFile* respectively).")



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
