#-*- coding: latin-1 -*-

### File: 2_Load_table_data_multiGear2.R
### Time-stamp: <2018-09-27 17:41:40 yreecht>
###
### Created: 29/06/2016	14:27:57
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

if ( ! all(is.element(catchData$LD_GearID, names(gearWidth))))
{
    Stop("Effort undefined for some gear(s): \"",
         paste(unique(catchData$LD_GearID[ ! is.element(catchData$LD_GearID, names(gearWidth))]),
               collapse = "\", \""), "\"")
}

catchDataTmp <- by(data = catchData,
                   INDICES = as.list(catchData[ , c("LD_HaulNo", "CatchUnitID")]),
                   FUN = function(x, gearWidth)
                {
                    ## In case effort is undefined for some gear:
                    x <- x[ , is.element(x[ , "LD_GearID"], names(gearWidth))]

                    ## Aggregation of catches and efforts in one entry:
                    res <- x[1, ]

                    dens <- x[ , "Caught"] / (gearWidth[as.character(x[ , "LD_GearID"])] ^ 2)
                    res[ , "Caught"] <- mean(dens, na.rm = TRUE)
                    res[ , "effort"] <- 1

                    res[ , "LD_GearID"] <- "MixtUnweight"

                    return(res)
                },
                gearWidth = gearWidth,
                simplify=FALSE)

catchData <- do.call(rbind, catchDataTmp)

## Adding the effort to station data:
catchDataSub <- subset(catchData,
                       is.element(as.character(LD_HaulNo),
                                  as.character(stationData@data$TRACKID)))

stationData@data$effort[match(as.character(catchDataSub$LD_HaulNo),
                              as.character(stationData@data$TRACKID))] <- catchDataSub$effort







### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
