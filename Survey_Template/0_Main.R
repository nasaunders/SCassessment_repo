#-*- coding: latin-1 -*-

### File: 0_Main.R
### Time-stamp: <2018-11-01 15:52:12 yreecht>
###
### Template created: 14/03/2016	14:21:44
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Routines for survey analyses:
###   * density and biomass assessment by species.
###   * size distribution.
###
### ==================================
### Survey:
###
### <Add your survey description here>
####################################################################################################

rm(list=ls())

## =================================================================================================
## Files & pathes:
scriptDir <- "y:/Analyses/Surveys/Scripts"      # Change to fit your installed routines.

workingDir <- "y:/Analyses/Surveys/<your_working_directory>"
## Path for storing the results:
ResultsPath <- "./Results"              # Usually a relative path to workingDir.


## Station shapefile:
shpDataDir <- file.path("//Galwayfs03/FishData/INSHORE/Species/",
                        "<...>")
stationLayer <- "<shape_file_base_name>" # Shapefile name without extension.

maskLayer <- NULL                       # Selection (rough) mask: to differentiate several zones
                                        # and/or exclude stations (NULL or "" to ignore).
                                        # (Shapefile name without extension, in shpDataDir).

## ## If you want to use an alternative map and/or add a grid, specify the following:
## shpOtherDir <- ""
## coastlineLayer <- ""
## gridLayer <- ""
## unsuitableLayer <- ""    # Unsuitable grounds layer... to be added to the costline (default, if relevant).
## bedLayer <- ""           # Potential/forced (cf. options) bed for interpolations/extrapolations.

## Pathes for catch and life history traits files:
otherDataDir <-  file.path("//Galwayfs03/FishData/INSHORE/Species/",
                           "<...>")

catchDataFile <- "<Catch_data_file>.csv"           # Exported from the "BiologicalData" tab*
biologicalDataFile <- "<Biological_data_file>.csv" # Exported from the "Catch" tab*
                                        # *from the inshore database upload template.

biologicalDataFilePrev <- NULL          # "<Previous_Biological_data_file>.csv"
                                        # previous biological data for comparison
                                        # Ignored if NULL or not existent.

totalEffortPrev <- NULL                 # Total effort of the previous survey (in m^2)
                                        # Used for the standardized size distribution (prev. year).
                                        # Ignored if NULL

## Working directory:
setwd(workingDir) # Should be the one containing this file.

## =================================================================================================
## Gear information:
gearWidth <- <width in meters>          # Dredge/quadrat width in m.
                                        # for multi-gear use, e.g.:
                                        # c("Scientific quadrat" = 0.5,
                                        #   "Hand rakes" = sqrt(2))
                                        #
gearEfficiency <- NULL                  # Dredge efficiency (named vector by species,
                                        # e.g. c("Ostrea edulis" = 0.1753)
                                        #   applied to all species if unnamed,
                                        #   ignored if NULL or undefined for a species).

MLS <- NULL                             # Minimum Landing Sizes (by species).
                                        # e.g. c("Ostrea edulis" = 76)

## List of size class selections on which indices have to be calculeted as well:
sizeSelections <- list()                # e.g. list(c(species = "Ostrea edulis", min = 76, max = Inf))
                                        # Ignored if empty list

## =================================================================================================
## Growth parameters:

## From a data file (in the otherDataDir directory):
LWmeasurementsFile <- NULL              # "<LW_measurements_file>.csv"
                                        # In case you want to calculate the relationship from data.
                                        # (csv file with three columns "species", "size" and "weight",
                                        #  in that exact order). Ignored if NULL.

## Matrix of parameters of the relationship W = a * L^b, by species
##     / or a file name from where it is loaded
##     (overriden if a relevant measurement file is defined):
LWparam <- NULL                         # "<LW_params_file>.csv"
                                        # OR a matrix:
                                        # matrix(data = c(0.00001, 3.7232 # Params for Ostrea edulis
                                        #                 ),
                                        #        byrow = TRUE,
                                        #        ncol = 2,
                                        #        dimnames = list(c("Ostrea edulis"),
                                        #                        c("a", "b")))

## =================================================================================================
## Options:

## General options:
options(warn=TRUE                       # Printing warnings immediatly.
        , surveyTargetOnly = TRUE       # Produce indices for target species only.
        , surveyCoastMap = "Ireland"    # Any other value than Ireland => WorldMap
                                        #     (but overriden if an alternative map file is specified above).
        , surveyMapSource = "bing"      # "shapefile" (defined file or default) or any map type supported
                                        #     by the package "OpenStreetMap" ("osm", "osm-bw", "opencyclemap",...).
                                        #     "bing" by default. Then requires an active internet connection.
        , surveyBedForce = TRUE         # Whether to force the extend of the bed to bedLayer (if relevant).
        )

## Inverse distance interpolation options:
options(#
        surveyIDW.alpha = 500                 # Max distance (in projection unit, usually m)
                                              #     between stations to avoid concave parts
                                              #     in the survey zone definition.
        , surveyIDW.cellsize = 10             # Cell size (in projection unit, usually m).
        , surveyIDW.bufferWidth = 20          # Buffer width (in projection unit, usually m) around the
                                              #     stations for the survey zone definition.
        , surveyIDW.nmin = 0                  # Min number of neighbours.
        , surveyIDW.nmax = Inf                # Max number of neighbours.
        , surveyIDW.idp = 2                   # Inverse distance power (2 for an ArcMap like behaviour).
        , surveyIDW.n.categ = 3               # Number of breaks. Alternatively a breaks vector.
        , surveyIDW.break.method = "quantileStation" # Method for calculating breaks.
        , surveyIDW.log.scale = "Default"
        , surveyIDW.breaksList = NULL         # A named (by variable) list of breaks.
                                              #     Override breaks calculations for the variables present
                                              #     in the list.
        )

## Graphical options:
options(## surveyPlot.colIDW = function(n) # Color palette for the idw interpolations.
        ## {
        ##     rev(heat.colors(n))
        ## },
        ## surveyPlot.colGround = "#D7D79E", # Color for the ground.
        ## surveyPlot.colSea = "#d6edf8",    # Color for the sea.
        ## surveyPlot.main = FALSE,
        ## surveyPlot.addStations = TRUE,
        ## surveyPlot.colStations = "darkred",
        surveyPlot.dev = c("png") # Device(s) (in "X11", "pdf", "png", "jpg"))
        )

## #######################
## Geostatistics options
## (variogram parameters):

## options(surveyGeostat.vg.cutoff = 1000,
##         surveyGeostat.vg.width = 80)


## =================================================================================================
## Initialisation of the working environment:
source(file.path(scriptDir,
                 "0_Load_packages.R"),
       encoding="latin1")

source(file.path(scriptDir,
                 "0_Initialise_environment.R"),
       encoding="latin1",
       chdir = TRUE)

source(file.path(scriptDir,
                 "0_Paths.R"),
       encoding="latin1")

## =================================================================================================
## Loading and summarising spatial data (to notably fill the database upload template):
source(file.path(scriptDir,
                 "1_Load_spatial_data.R"),
       encoding="latin1")

## Here you can finish to prepare the catch and biological data for database upload
## (also used in these scripts) using the output file.

## =================================================================================================
## Loading catch and biological data and matching them with spatial data:
source(file.path(scriptDir,
                 "2_Load_table_data.R"),
       encoding="latin1")

source(file.path(scriptDir,
                 "2_Weights.R"),
       encoding="latin1")

source(file.path(scriptDir,
                 "2_Size_class_metrics.R"),
       encoding="latin1")

## =================================================================================================
## Analyses based on biological traits:

source(file.path(scriptDir,
                 "3_Size_distribution.R"),
       encoding="latin1")

## =================================================================================================
## Analyses of spatial data:
source(file.path(scriptDir,
                 "4_Density_IDW_interpolation.R"),
       encoding="latin1")

## At this stage you can change the break calculation/ force breaks and/or
##     modify the graphical parameters and re-run from here until you get
##     something satisfactory (some examples following):

## options(surveyIDW.breaksList = NULL)

## options(surveyIDW.breaksList = list("AbundanceDensity_Ostrea_edulis" = c(0, 0.1, 1, 2.5, Inf),
##                                     "AbundanceDensityCorr_Ostrea_edulis" = c(0, 0.5, 2.5, 5, Inf)))

## options(surveyIDW.break.method = "range", # Method for calculating breaks.
##         surveyIDW.log.scale = "Default")

source(file.path(scriptDir,
                 "4_Strata.R"),
       encoding="latin1")

## ## =================================================================================================
## ## Stratified frequentist estimate (deprecated, kept for comparison purpose):

## source(file.path(scriptDir,
##                  "4_Density_IDW_plot.R"),
##        encoding="latin1")

## source(file.path(scriptDir,
##                  "5_Raised_abundances.R"),
##        encoding="latin1")

## =================================================================================================
## Stratified Bayesian biomass assessment:

## ## Bayesian assessment options:
## options(surveyBayes.n.burnin = 2000,
##         surveyBayes.n.iter = 10000,
##         surveyBayes.n.chains = 3,
##         surveyBayes.models = c("N", "LN", "G"))

source(file.path(scriptDir,
                 "5_Raised_abundances_Bayes.R"),
       encoding="latin1")

## =================================================================================================
## Geostatistics assessment:

## ## Variogram fitting options:
## options(surveyGeostat.vg.cutoff = 3000,
##         surveyGeostat.vg.width = 80,
##         surveyGeostat.vg.autobreaks = FALSE,
##         surveyGeostat.vg.likelihood = TRUE,
##         surveyGeostat.vg.lik.method = "ML",
##         surveyGeostat.CGS.nsim = 500,
##         surveyGeostat.CGS.nmax = 30)

source(file.path(scriptDir,
                 "6_Kriging_biomass.R"),
       encoding="latin1")


## =================================================================================================
## Add your own statistics here:

## summary(stationData@data$AbundanceDensityCorr_Ostrea_edulis)
## summary(subset(biologicalData,
##                LD_SpeciesID == "Ostrea edulis"))
## sd(subset(biologicalData,
##           LD_SpeciesID == "Ostrea edulis")$Size)

## sum(subset(catchData,
##            SpeciesID == "Ostrea edulis")$Caught)

## sum(! is.na(subset(biologicalData,
##                    LD_SpeciesID == "Ostrea edulis")$Size))

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
