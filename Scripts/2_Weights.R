#-*- coding: latin-1 -*-

### File: 3_weights.R
### Time-stamp: <2018-09-12 13:17:08 yreecht>
###
### Created: 23/03/2016	17:27:36
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################
.glmRes <- list()

if ((exists("LWmeasurementsFile") &&
     ! is.null(LWmeasurementsFile)) ||
    sum( ! is.na(biologicalData$Size) & ! is.na(biologicalData$Weight)) > 30)
{
    if (exists("LWmeasurementsFile") &&
        ! is.null(LWmeasurementsFile))
    {
        LWdata <- read.csv(file = file.path(otherDataDir, LWmeasurementsFile))
    }else{
        if (isTRUE(getOption("surveyWL.all.obs")))
        { ## Using all observations:
            LWdata <- biologicalDataTot[! is.na(biologicalDataTot$Size) &
                                        ! is.na(biologicalDataTot$Weight),
                                        c("LD_SpeciesID", "Size", "Weight")]
        }else{ ## Using active stations only:
            LWdata <- biologicalData[! is.na(biologicalData$Size) &
                                     ! is.na(biologicalData$Weight),
                                     c("LD_SpeciesID", "Size", "Weight")]
        }
    }

    colnames(LWdata) <- c("species", "size", "weight")

    ## Calculate the weight-length relationships:
    source(file.path(scriptDir, "2_Weight-length_relationships.R"))

    if (exists("LWparam") &&
        ! is.null(LWparam))
    {
        LWparam <- rbind(LWparamCalc,
                         LWparam[ ! is.element(row.names(LWparam),         # Overwritting former parameters
                                               row.names(LWparamCalc)), ]) # if calculated from data.
    }else{
        LWparam <- LWparamCalc
    }
}else{}


## lines(1:250, LWparamCalc[1, "a"] * (1:250)^LWparamCalc[1, "b"], lty = 2)


if (exists("LWparam") && ! is.null(LWparam) && any(is.na(biologicalData$Weight)))
{
    tmpBData <- sapply(unique(as.character(biologicalData$LD_SpeciesID)),
                       function(sp, data, params)
                   {
                       if (sp == "") return(NULL) # blank rows to remove.

                       subData <- subset(data, LD_SpeciesID == sp)

                       if (is.element(sp, row.names(params)))
                       {
                           subData[is.na(subData[ , "Weight"]), "Weight"] <-
                                    params[sp, "a"] * (subData[is.na(subData[ , "Weight"]),
                                                               "Size"] ^ params[sp, "b"])
                       }else{}

                       return(subData)
                   },
                       data = biologicalData,
                       params = LWparam,
                       simplify = FALSE)

    biologicalData <- do.call(rbind, tmpBData)
}else{
    message("\n## Nothing to calculate!")
}






### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
