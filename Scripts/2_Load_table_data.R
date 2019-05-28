#-*- coding: latin-1 -*-

### File: 2_Load_table_data.R
### Time-stamp: <2018-11-01 11:08:26 yreecht>
###
### Created: 16/03/2016	12:12:30
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Script for loading catch and biological data.
####################################################################################################


catchData <- read.csv(file = file.path(otherDataDir, catchDataFile))

biologicalData <- read.csv(file = file.path(otherDataDir, biologicalDataFile))
biologicalData$LD_EventStartDate <- as.Date(as.character(biologicalData$LD_EventStartDate),
                                            format = "%d/%m/%Y")

stationData@data[ , "TRACKID"] <- as.character(stationData@data[ , "TRACKID"])

catchData[ , "LD_HaulNo"] <- as.character(catchData[ , "LD_HaulNo"])
catchData <- subset(catchData,
                    is.element(catchData[ , "LD_HaulNo"],
                               stationData@data[ , "TRACKID"]))

biologicalData[ , "LD_HaulNo"] <- as.character(biologicalData[ , "LD_HaulNo"])

biologicalDataTot <- biologicalData     # May be needed for W-L relationship.

biologicalData <- subset(biologicalData,
                         is.element(biologicalData[ , "LD_HaulNo"],
                                    stationData@data[ , "TRACKID"]))
biologicalData <- subset(biologicalData,
                         is.element(biologicalData[ , "LD_HaulNo"],
                                    catchData[ , "LD_HaulNo"]))


if (exists("LWparam"))
{
    if (! is.matrix(LWparam) && is.character(LWparam))
    {
        LWparam <- as.matrix(read.csv(file = file.path(otherDataDir, LWparam[1]), row.names = 1))
    }
}else{
    LWparam <- NULL
}

if (length(gearWidth) > 1 &&
    ! is.null(names(gearWidth)))
{
    if (is.element("SpatialLinesDataFrame", class(stationData)))
    {
        ## Match gear width and station data to calculate effort:
        stationData@data$gear <- catchData[match(as.character(stationData@data$TRACKID),
                                                 as.character(catchData$LD_HaulNo)),
                                           "LD_GearID"]
        stationData@data$effort <- stationData@data$length * gearWidth[as.character(stationData@data$gear)]
    }else{                              # Point data: the effort is squared width of the gear.
        source(file.path(getOption("survey.scriptPath"), "2_Load_table_data_multiGear2.R"),
               encoding="latin1")
    }
}

head(catchData)
head(biologicalData)


if (exists("biologicalDataFilePrev") &&
    ! is.null(biologicalDataFilePrev))
{
    biologicalDataPrev <- read.csv(file = file.path(otherDataDir, biologicalDataFilePrev))

    if (ncol(biologicalDataPrev) == 2)
    {
        colnames(biologicalDataPrev) <- c("LD_SpeciesID", "Size")
    }else{
        biologicalDataPrev$LD_EventStartDate <-
            as.Date(as.character(biologicalDataPrev$LD_EventStartDate),
                    format = "%d/%m/%Y")
    }

    head(biologicalDataPrev)
}


## Analysed species:
species <- unique(as.character(catchData[is.element(catchData[ , "IsTargetY.N"],
                                                    c("Yes", "yes", "Y", "y")) |
                                         (! getOption("surveyTargetOnly")) ,
                                         "SpeciesID"]))

## ...and corresponding metrics:
speciesCatchUnit <- with(subset(catchData,
                                is.element(SpeciesID,
                                           species)),
                         table(as.character(SpeciesID), as.character(CatchUnitID))) > 0

speciesUnit <- speciesCatchUnit

## Data check here!
source(file.path(getOption("survey.scriptPath"), "2_Check_table_data.R"),
       encoding="latin1")

## Gear efficiency:
if (exists("gearEfficiency"))
{
    if ( ! is.null(gearEfficiency))
    {
        if (is.null(names(gearEfficiency)))
        {
            gearEfficiency <- rep(gearEfficiency,
                                  length.out = nrow(speciesUnit))

            names(gearEfficiency) <- row.names(speciesUnit)
        }else{}
    }else{}
}else{
    gearEfficiency <- NULL
}

## MLS
if (exists("MLS"))
{
    if ( ! is.null(MLS))
    {
        if (is.null(names(MLS)))
        {
            MLS <- rep(MLS,
                       length.out = nrow(speciesUnit))

            names(MLS) <- row.names(speciesUnit)
        }else{}
    }else{}
}else{
    MLS <- NULL
}


## Matching stations and kept data:
for (sp in row.names(speciesCatchUnit))
{
    for (unit in colnames(speciesCatchUnit))
    {
        if (speciesCatchUnit[sp, unit])
        {
            message("Adding metrics \"", unit, "\" for the species \"", sp, "\"")

            colName <- gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
                            paste(unit, sp))

            dataTmp <- subset(catchData,
                              SpeciesID == sp & CatchUnitID == unit)

            stationData@data[ , colName] <- dataTmp[match(stationData@data[ , "TRACKID"],
                                                          dataTmp[ , "LD_HaulNo"]) ,
                                                    "Caught"]

            ## Nothing caught interpreted as 0 (if station valid, i.e. in catchData):
            stationData@data[is.na(stationData@data[ , colName]) &
                             is.element(as.character(stationData@data[ , "TRACKID"]),
                                        as.character(catchData[ , "LD_HaulNo"])),
                             colName] <- 0


        }else{}
    }
}


## Corrections of metrics according to gear efficiency:
for (sp in row.names(speciesCatchUnit))
{
    for (unit in colnames(speciesCatchUnit))
    {
        if (speciesCatchUnit[sp, unit])
        {
            ## Efficiency correction:
            if (sp %in% names(gearEfficiency))
            {
                ## Original column:
                colName <- gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
                                paste(unit, sp))

                ## New metrics & column:
                unitCorr <- paste(unit, "Corr", sep="")

                message("Adding metrics \"", unitCorr,
                        "\" (corrected for gear efficiency) for the species \"", sp, "\"")

                colNamesEff <- gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
                                    paste(unitCorr , sp))

                stationData@data[ , colNamesEff] <- stationData@data[ , colName] / gearEfficiency[sp]

                ## Adds the metric to speciesUnit where relevant:
                if (! is.element(unitCorr,
                                 colnames(speciesUnit)))
                {
                    speciesUnit <- cbind(speciesUnit,
                                         matrix(data = FALSE,
                                                nrow = nrow(speciesUnit), ncol = 1,
                                                dimnames = list(row.names(speciesUnit),
                                                                unitCorr)))
                }else{}

                ## Turns the speciesUnit corresponding entry to TRUE:
                speciesUnit[sp, unitCorr] <- TRUE
            }else{}
        }else{}
    }
}


## Adds the density metrics:
for (sp in row.names(speciesUnit))
{
    for (unit in colnames(speciesUnit))
    {
        if (speciesUnit[sp, unit])
        {
            ## Original column:
            colName <- gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
                            paste(unit, sp))

            ## Density metrics:
            if (grepl("^(Number|KGs|NumberCorr|KGsCorr)", unit))
            {
                unitDens <- switch(unit,
                                   "Number" = {"AbundanceDensity"},
                                   "NumberCorr" = {"AbundanceDensityCorr"},
                                   "KGs" = {"WeightDensity"},
                                   "KGsCorr" = "WeightDensityCorr")

                colNameDens <- gsub("[[:blank:]+-]", "_",
                                    paste(unitDens, sp))

                message("Adding metrics \"", unitDens, "\" for the species \"", sp, "\"")

                ## Density / effort (usually m^2):
                stationData@data[ , colNameDens] <- stationData@data[ , colName] / stationData@data$effort

                ## The unit must be added to the unitSpecies table if does not exist yet:
                if ( ! is.element(unitDens, colnames(speciesUnit)))
                {
                    speciesUnit <- cbind(speciesUnit,
                                         matrix(data = FALSE,
                                                nrow = nrow(speciesUnit), ncol = 1,
                                                dimnames = list(row.names(speciesUnit),
                                                                unitDens)))
                }else{}

                ## Setting it to TRUE for that species:
                speciesUnit[sp, unitDens] <- TRUE
            }else{}
        }else{}
    }
}

## Identification of fields for which interpolations are calculated:
metrics <- sapply(colnames(speciesUnit),
                  function(unit, speciesUnit)
              {
                  sapply(row.names(speciesUnit),
                         function(sp, unit, speciesUnit)
                     {
                         if (speciesUnit[sp, unit])
                         {
                             return(data.frame("field" = gsub("[[:blank:]+-]", "_",
                                                              paste(unit, sp)),
                                               "species" = sp,
                                               "metrics" = unit,
                                               "sizeClass" = FALSE,
                                               "min" = NA,
                                               "max" = NA,
                                               stringsAsFactors = FALSE))
                         }else{
                             return(NULL)
                         }
                     },
                         speciesUnit = speciesUnit, unit = unit, simplify = FALSE)
              },
                  speciesUnit = speciesUnit, simplify = FALSE)

## Getting a data.frame from the result:
metrics <- do.call(rbind, unlist(metrics, recursive=FALSE))
row.names(metrics) <- metrics[ , "field"]






### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
