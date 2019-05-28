#-*- coding: latin-1 -*-

### File: 4_Density_IDW_interpolation.R
### Time-stamp: <2018-11-01 11:08:25 yreecht>
###
### Created: 16/03/2016	10:05:56
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Interpolation of abundance and biomass densities.
####################################################################################################

## Identification of fields for which interpolations are calculated:
interpMetrics <- metrics[is.element(metrics[ , "metrics"],
                                    c("AbundanceDensity", "WeightDensity",
                                      "AbundanceDensityCorr", "WeightDensityCorr")), ]


## Interpolations:
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

source(file.path(getOption("survey.scriptPath"), "./4_Zone_plot.R"),
       encoding="latin1")


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
