#-*- coding: latin-1 -*-

### File: 0_Main_data_preparation.R
### Time-stamp: <2019-01-18 11:47:48 yreecht>
###
### Created: 06/09/2017	09:04:49
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## ##################################################
## Paths:

## Routines:
scriptPath <- "y:/Analyses/Surveys/Scripts"

options("FilePrefix" = "",
        "FileSuffix" = "")

## Data:
measurementsFile <- file.path("//Galwayfs03/FishData/INSHORE/Species",
                              "<...>")  # measurements, see template.

weighLengthFile <- NULL                 # optional, external weight-length
                                        #    relationship.

samplesFile <-  file.path("//Galwayfs03/FishData/INSHORE/Species",
                          "<...>")      # Samples, see template.


## Output:
outputDir <- file.path("//Galwayfs03/FishData/INSHORE/Species",
                       "<...>/")

## ##################################################
## Running the routines:
source(file.path(scriptPath, "0_Data_preparation.R"), chdir = TRUE)



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
