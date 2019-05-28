#-*- coding: latin-1 -*-

### File: 2_Size_class_metrics.R
### Time-stamp: <2019-02-01 12:33:50 yreecht>
###
### Created: 08/04/2016	10:28:36
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## Data by size selection:
sizeClasseData <- sapply(sizeSelections,
              function(x, speciesCatchUnit, biologicalData, catchData)
          {
              if (is.element(x[1], row.names(speciesCatchUnit)) &&
                  any(speciesCatchUnit[x[1], ]))
              {
                  resData <- catchData[is.element(catchData[ , "SpeciesID"],
                                                  x[1]),
                                       c("LD_HaulNo", "SpeciesID", "CatchUnitID", "Caught")]

                  resList <- list()

                  for (unit in colnames(speciesCatchUnit))
                  {
                      resDataUnit <- subset(resData,
                                            is.element(CatchUnitID, unit))

                      if (speciesCatchUnit[x[1], unit])
                      {
                          ## Biological data subsets:
                          biolSp <- subset(biologicalData,
                                           is.element(LD_SpeciesID, x[1]) & ! is.na(Size))

                          biolSpSz <- subset(biolSp,
                                             Size >= as.numeric(x[2]) & Size <= as.numeric(x[3]))

                          ## Total measured (at the individual level) by tow:
                          totMeasured <- tapply(X = biolSp$Weight,
                                                INDEX = as.character(biolSp$LD_HaulNo),
                                                FUN = switch(unit,
                                                             "Number" = length,
                                                             "KGs" = function(x){sum(x, na.rm = TRUE)},
                                                             function(x){NA}))

                          resDataUnit <- merge(resDataUnit,
                                               data.frame("LD_HaulNo" = names(totMeasured),
                                                          "Measured" = totMeasured),
                                               all.x = TRUE)

                          ## Measured for this size class (at the individual level) by tow:
                          sizeMeasured <- tapply(X = biolSpSz$Weight,
                                                 INDEX = as.character(biolSpSz$LD_HaulNo),
                                                 FUN = switch(unit,
                                                             "Number" = length,
                                                             "KGs" = function(x){sum(x, na.rm = TRUE)},
                                                             function(x){NA}))

                          resDataUnit <- merge(resDataUnit,
                                               data.frame("LD_HaulNo" = names(sizeMeasured),
                                                          "MeasuredSz" = sizeMeasured),
                                               all.x = TRUE)

                          ## Non measured are zeros:
                          resDataUnit$Measured[is.na(resDataUnit$Measured)] <- 0
                          resDataUnit$MeasuredSz[is.na(resDataUnit$MeasuredSz)] <- 0

                          ## Global ratio of the aggregated metrics for that size class
                          ## (used for extrapolation when no[t enough] individual[s] measured):
                          meanRatio <- sum(resDataUnit$MeasuredSz) / sum(resDataUnit$Measured)

                          ## Corrected catches:
                          resDataUnit$CaughtSz <-
                              ifelse(test = (resDataUnit$Measured == 0 |
                                             (resDataUnit$Measured < 5 &
                                              resDataUnit$Caught >= 10 &
                                              resDataUnit$MeasuredSz == 0)),
                                     yes = meanRatio * resDataUnit$Caught,
                                     no = (resDataUnit$MeasuredSz *
                                           resDataUnit$Caught / resDataUnit$Measured))

                          ## The metrics entry for this species/metrics/size class combination:
                          resMetrics <- data.frame(field = gsub("[[:blank:]+-]", "_",
                                                                paste(unit,
                                                                      paste(as.numeric(x[2:3]), collapse = " "),
                                                                      x[1])),
                                                   species = x[1],
                                                   metrics = unit,
                                                   sizeClass = TRUE,
                                                   min = as.numeric(x[2]),
                                                   max = as.numeric(x[3]))

                          row.names(resMetrics) <- resMetrics$field

                          ## We return a list of several results by species/metrics/size class combination:
                          if (all(is.na(resDataUnit$CaughtSz)))
                          {
                              ## In case no individual measurements were done, NaNs are produced
                              ## and the metrics must be discarded:
                              res <- NULL
                          }else{
                              res <- list(metrics = resMetrics,
                                          catchData = resDataUnit,
                                          biologicalData = biolSpSz)

                              ## List of results by unit:
                              resList <- c(resList, list(res))
                          }
                      }else{}
                  }
                  return(resList)
              }else{return(NULL)}
          },
              speciesCatchUnit = speciesCatchUnit,
              biologicalData = biologicalData,
              catchData = catchData,
              simplify = FALSE)

## Getting a list with one element by species/metrics/size class combination:
sizeClasseData <- unlist(x = sizeClasseData, recursive=FALSE, use.names=TRUE)
sizeClasseData <- sizeClasseData[! is.null(sizeClasseData)]
if (!is.null(sizeClasseData))
{
    names(sizeClasseData) <- sapply(sizeClasseData, function(x) x[["metrics"]]$field)
}
## Adding the corresponding rows in metrics:
metrics <- rbind(metrics,
                 do.call(rbind,
                         sapply(sizeClasseData,
                                function(x) x[["metrics"]],
                                simplify = FALSE)))

## sapply(metrics, class)
## summary(stationData@data)

## test <- stationData@data

addCol <- sapply(row.names(subset(metrics, sizeClass)),
                 function(i)
             {
                 x <- metrics[i, ]
                 res <- sizeClasseData[[i]]$catchData[ , c("LD_HaulNo", "CaughtSz")]
                 colnames(res) <- c("TRACKID", x[ , "field"])
                 res[ , "TRACKID"] <- as.character(res[ , "TRACKID"])
                 return(res)
             },
                 simplify = FALSE)


addCol <- Reduce(f = function(x, y){merge(x, y, all = TRUE)},
                 x = addCol)

## sapply(addCol, class)

if (! is.null(addCol))
{
    stationData@data <- cbind(stationData@data[ ,
                                               ! is.element(colnames(stationData@data),
                                                            tail(colnames(addCol), -1))],
                              addCol[match(as.character(stationData@data[ , "TRACKID"]),
                                           as.character(addCol[ , "TRACKID"])),
                                     -1,
                                     drop = FALSE])
}

## head(stationData@data)
## stationData@data[ , c("Track_ID", "TRACKID", "Number_Crassostrea_gigas", "Number_76_Inf_Crassostrea_gigas")]
## stationData@data[ , c("Track_ID", "TRACKID", "AbundanceDensity_Crassostrea_gigas", "AbundanceDensity_76_Inf_Crassostrea_gigas")]
## head(addCol)


## Corrections of size-class metrics according to gear efficiency:
for (i in row.names(subset(metrics, sizeClass)))
{
    sp <- metrics[i, "species"]
    unit <-  metrics[i, "metrics"]

    ## Efficiency correction:
    if (sp %in% names(gearEfficiency))
    {
        ## Original column:
        colName <- i

        ## New metrics & column:
        unitCorr <- paste(unit, "Corr", sep="")

        message("Adding metrics \"", unitCorr,
                "\" (corrected for gear efficiency) for the species \"", sp,
                "\" and size range [", metrics[i , "min"], ",", metrics[i , "max"], "]")

        colNamesEff <- gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
                            paste(unitCorr,
                                  metrics[i, "min"],
                                  metrics[i, "max"],
                                  sp))

        stationData@data[ , colNamesEff] <- stationData@data[ , colName] / gearEfficiency[sp]

        ## Updating the metrics table:
        metricsTmp <- metrics[i, ]
        metricsTmp[ , "metrics"] <- unitCorr
        metricsTmp[ , "field"] <- colNamesEff
        row.names(metricsTmp) <- colNamesEff

        metrics <- rbind(metrics, metricsTmp)
    }else{}
}
head(stationData@data)

## Adds the density metrics:
for (i in row.names(subset(metrics, sizeClass)))
{
    sp <- metrics[i, "species"]
    unit <-  metrics[i, "metrics"]

    ## Original column:
    colName <- i ## gsub("[[:blank:]+-]", "_", # field name by species, without blanks.
    ##    paste(unit, sp))

    ## Density metrics:
    if (grepl("^(Number|KGs|NumberCorr|KGsCorr)", unit))
    {
        unitDens <- switch(unit,
                           "Number" = {"AbundanceDensity"},
                           "NumberCorr" = {"AbundanceDensityCorr"},
                           "KGs" = {"WeightDensity"},
                           "KGsCorr" = "WeightDensityCorr")

        colNameDens <- gsub("[[:blank:]+-]", "_",
                            paste(unitDens,
                                  metrics[i, "min"],
                                  metrics[i, "max"],
                                  sp))

        message("Adding metrics \"", unitDens, "\" for the species \"", sp,
                "\" and size range [", metrics[i , "min"], ",", metrics[i , "max"], "]")

        ## Density / effort (usually m^2):
        stationData@data[ , colNameDens] <- stationData@data[ , colName] / stationData@data$effort

        ## Updating the metrics table:
        metricsTmp <- metrics[i, ]
        metricsTmp[ , "metrics"] <- unitDens
        metricsTmp[ , "field"] <- colNameDens
        row.names(metricsTmp) <- colNameDens

        metrics <- rbind(metrics, metricsTmp)
    }else{}
}

## Information:
message("\n## Preview of data:")
print(head(stationData, 3))

message("\n## Metrics / species table:")
print(speciesUnit)

densUnits <- colnames(speciesUnit)[grepl("^(AbundanceDensity|AbundanceDensityCorr|WeightDensity|WeightDensityCorr)", colnames(speciesUnit))]

message("\n## Breaks can be manually defined using",
        "\n\t options(surveyIDW.breaksList = list(\"<variable1>\" = c(<breaks>),",
        "\n\t                                     \"<variable2>\" = c(<breaks>),...))",
        "\n## for the following variables: \n\t\"",
        paste(## unlist(sapply(densUnits,
              ##               function(unit, speciesUnit)
              ##           {
              ##               gsub("[[:blank:]+-]", "_",
              ##                    paste(unit, row.names(speciesUnit)[speciesUnit[ , unit]]))
              ##           },
              ##               speciesUnit = speciesUnit))
              metrics[grepl("^(AbundanceDensity|AbundanceDensityCorr|WeightDensity|WeightDensityCorr)",
                            metrics[ , "metrics"]),
                      "field"],
              collapse = "\", \""),
        "\".")

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
