#-*- coding: latin-1 -*-

### File: 1_Data_preparation_checks.R
### Time-stamp: <2019-01-18 15:01:26 yreecht>
###
### Created: 07/09/2017	14:07:51
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

errors <- c("brokenMissing" = FALSE,
            "stSamplesMissing" = FALSE,
            "sampleDate" = FALSE,
            "measDate" = FALSE,
            "dateMismatch" = FALSE,
            "gradeMismatch" = FALSE,
            "spMismatch" = FALSE)

if (any(is.na(samples$broken)))
{
    message("\nError: missing values in broken")

    errors["brokenMissing"] <- TRUE
}

## Basic checks on unicity of date/station pairs:
stationDate <- table(samples$stationid, samples$date)
stationDate[stationDate > 0] <- 1

stationDate2 <- table(measurements$stationid, measurements$date)
stationDate2[stationDate2 > 0] <- 1

## Station names mismatch (all those in measurements must be in samples):
stSamples <- row.names(stationDate)
stMeas <- row.names(stationDate2)
stJoined <- intersect(stSamples, stMeas)
stMeasOnly <- setdiff(stMeas, stSamples)

if (length(stMeasOnly))
{
    message("\nError: station(s) in measurements not present in the samples:\n    * ",
            paste(stMeasOnly,
                  collapse = "\n    * "))

    errors["stSamplesMissing"] <- TRUE
}

## Ambiguities might persist if within one vessel/team data:
if ( ! isTRUE(all.equal(range(apply(X = stationDate,
                                    MARGIN = 1,
                                    FUN = sum)),
                        c(1, 1))))
{
    message("\nError: sample station(s) over several dates:\n    * ",
            paste(row.names(stationDate)[apply(X = stationDate,
                                               MARGIN = 1,
                                               FUN = sum) > 1],
                  collapse = "\n    * "))

    errors["sampleDate"] <- TRUE
}

if ( ! isTRUE(all.equal(range(apply(X = stationDate2,
                                    MARGIN = 1,
                                    FUN = sum)),
                        c(1, 1))))
{
    message("\nError: measurement station(s) over several dates:\n    * ",
            paste(row.names(stationDate2)[apply(X = stationDate2,
                                                MARGIN = 1,
                                                FUN = sum) > 1],
                  collapse = "\n    * "))

    errors["measDate"] <- TRUE
}

## ...and entry errors may exist that lead to Samples - Measurements date mismatch:
if ( ! all(stationDate[stJoined, ] ==
           stationDate2[stJoined, ]))
{
    stationDateDiff <- stationDate[stJoined, ] == stationDate2[stJoined, ]
    ## For convenience, the station list is printed:
    message("\nError: station(s) with sample-measurement date mismatch:\n    * ",
            paste(row.names(stationDateDiff)[! apply(stationDateDiff, 1, all)],
                  collapse = "\n    * "))

    errors["dateMismatch"] <- TRUE
}

## Station - Grade mismatch:
stationGrade <- with(subset(samples,
                            species %in% measurements$species),
                     table(stationid, grade))
stationGrade[stationGrade > 0] <- 1

stationGrade2 <- table(measurements$stationid, measurements$grade)
stationGrade2[stationGrade2 > 0] <- 1



if ( any( ! colnames(stationGrade2) %in% colnames(stationGrade)) ||
     ! all(stationGrade[stJoined, colnames(stationGrade2)] ==
           stationGrade2[stJoined, ]))
{
    if (all(colnames(stationGrade2) %in% colnames(stationGrade)))
    {

        stationGradeDiff <- stationGrade[stJoined, ] == stationGrade2[stJoined, ]
        sampleMissing <- stationGrade2[stJoined, ] - stationGrade[stJoined, ] > 0
        measurementMissing <- stationGrade2[stJoined, ] - stationGrade[stJoined, ] < 0
        if (any(sampleMissing))
        {
            ## For convenience, the station list is printed:
            message("\nError: station(s) with sample-measurement grade mismatch:\n    * ",
                    paste(row.names(sampleMissing)[! apply( ! sampleMissing, 1, all)],
                          collapse = "\n    * "))

            errors["gradeMismatch"] <- TRUE
        }else{}

        if (any(measurementMissing))
        {
            ## For convenience, the station list is printed:
            warning("\nGrade(s) with no measurement in station(s):\n    * ",
                    paste(row.names(measurementMissing)[! apply( ! measurementMissing, 1, all)],
                          collapse = "\n    * "))
        }else{}
    }else{
        extraGrades <- colnames(stationGrade2)[! colnames(stationGrade2) %in% colnames(stationGrade)]
        ## For convenience, extra grades in measurements are printed:
        message("\nError: grades in the measurements, not present in the samples:\n    * ",
                paste(extraGrades, collapse = "\n    * "))

        errors["gradeMismatch"] <- TRUE
    }
}

## Station - Species mismatch:
stationSpecies <- table(samples$stationid, samples$species)
stationSpecies[stationSpecies > 0] <- 1

stationSpecies2 <- table(measurements$stationid, measurements$species)
stationSpecies2[stationSpecies2 > 0] <- 1



## if ( ! all(stationSpecies[row.names(stationSpecies2), ] ==
##            stationSpecies2))
## {
##     stationSpeciesDiff <- stationSpecies[row.names(stationSpecies2), ] == stationSpecies2
##     ## For convenience, the station list is printed:
##     message("\nError: station(s) with sample-measurement species mismatch:\n    * ",
##             paste(row.names(stationSpeciesDiff)[! apply(stationSpeciesDiff, 1, all)],
##                   collapse = "\n    * "))
##
##     errors["spMismatch"] <- TRUE
## }

## Station-species-grade mismatch

samplesSSG <- unique(subset(samples, select = c("species", "stationid", "grade")))
measurementsSSG <- unique(subset(measurements, select = c("species", "stationid", "grade")))

idxNoMatch <- which(is.na(match(x = apply(measurementsSSG, 1, paste, collapse = "+-+"),
                                table = apply(samplesSSG, 1, paste, collapse = "+-+"),
                                nomatch=NA)))


if (length(idxNoMatch))
{
    extraSSG <- apply(measurementsSSG[idxNoMatch, ], 1,
                      function(x)
                      {
                          paste(paste0(c("Sp:", "St:", "Gd:"),
                                       x),
                                collapse = ", ")
                      })
    message("\nError: Species-station-grade combination(s) in the measurements, not present in the samples:\n    * ",
            paste(extraSSG, collapse = "\n    * "),
            "\n")

    errors["gradeMismatch"] <- TRUE
}


## ##################################################
## Apply errors:
if (any(errors))
{
    if (errors["brokenMissing"]) stop("Samples: missing broken values!")
    if (errors["stSamplesMissing"]) stop("Samples: missing station present in measurements!")
    if (errors["sampleDate"]) stop("Samples: stations over several dates!")
    if (errors["measDate"]) stop("Measurements: stations over several dates!")
    if (errors["dateMismatch"]) stop("Samples - Measurements date mismatch(es)!")
    if (errors["gradeMismatch"]) stop("Samples - Measurements grade mismatch(es)!")
    if (errors["spMismatch"]) stop("Samples - Measurements species mismatch(es)!")
}
### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
