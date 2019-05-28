#-*- coding: latin-1 -*-

### File: 5_Raised_abundances.R
### Time-stamp: <2018-11-01 11:08:24 yreecht>
###
### Created: 22/03/2016	17:32:52
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


alph <- 0.05


stats <- sapply(names(densInterpPolyFusion),
                function(i, SP, stationData)
            {
                SP <- SP[[i]]

                ## Density category of every station (represented by its centroid):
                categs <- gCovers(spgeom1 = SP,
                                  spgeom2 = gCentroid(spgeom = stationData, byid = TRUE),
                                  byid=TRUE)

                ## Working on only stations with a category (some may be excluded if you used a mask):
                idx <- apply(categs, 1, sum) == 1

                stationsTmp <- stationData[idx, ]

                stationCateg <- data.frame("stationID"=stationsTmp$TRACKID,
                                           "categDens"=factor(colnames(categs)[apply(categs[idx, ,
                                                                                            drop = FALSE],
                                                                                     1, which)],
                                                              levels=colnames(categs))) ## [!!!] check that still works.

                ## Append categories to the temporary station data:
                stationsTmp@data <- cbind(stationsTmp@data[ ,
                                                           !is.element(colnames(stationsTmp@data),
                                                                       colnames(stationCateg))],
                                          stationCateg)

                head(stationsTmp@data)

                ## Statistics by strata:
                statDens <- data.frame("Area" = gArea(SP,
                                                      byid = TRUE)[levels(stationsTmp$categDens)],
                                       "D" = tapply(stationsTmp[[i]],
                                                    stationsTmp$categDens,
                                                    mean, na.rm=TRUE),
                                       "sdD" = tapply(stationsTmp[[i]],
                                                      stationsTmp$categDens,
                                                      sd, na.rm=TRUE),
                                       "nbD" = tapply(stationsTmp[[i]],
                                                      stationsTmp$categDens,
                                                      function(x)sum( ! is.na(x))),
                                       "seD" = tapply(stationsTmp[[i]],
                                                      stationsTmp$categDens,
                                                      sd, na.rm=TRUE) /
                                           sqrt(tapply(stationsTmp[[i]],
                                                       stationsTmp$categDens,
                                                       function(x)sum( ! is.na(x)))))

                ## Aggregated statistics:
                totalD <- c("Area" = sum(statDens[ , "Area"], na.rm = TRUE),
                            "D" = weighted.mean(x = statDens[ , "D"],
                                                w = statDens[ , "Area"], na.rm=TRUE),
                            "dsD" = NA,
                            "nbD" = sum(statDens[ , "nbD"], na.rm = TRUE),
                            ## sd(var_st) (i.e. SE) according to Cochran (1953) and Fieller (1954):
                            "seD" = sqrt(sum(statDens[ , "Area"]^2 * statDens[ , "sdD"]^2 /
                                         statDens[ , "nbD"], na.rm = TRUE) /
                                         (sum(statDens[ , "Area"], na.rm = TRUE) ^ 2)))
                statDens <- rbind(statDens,
                                  Total = totalD)

                ## CI distance to the mean (given the value of alpha):
                statDens$dAlpha <- qnorm(1 - getOption("surveyStat.alpha") /
                                         2) * statDens[ , "seD"]
                statDens$dAlpha[is.na(statDens$dAlpha)] <- 0

                statDens <- cbind(statDens,
                                  "Ab" = (statDens$D  * statDens[ , "Area"] *
                                          ifelse(is.element(interpMetrics[i, "metrics"],
                                                            c("WeightDensity", "WeightDensityCorr")),
                                                 1e-3, 1)),
                                  "abAlpha" = (statDens$dAlpha  * statDens[ , "Area"] *
                                               ifelse(is.element(interpMetrics[i, "metrics"],
                                                            c("WeightDensity", "WeightDensityCorr")),
                                                      1e-3, 1)))

                ## Calculation of weight statistics for combining with densities if relevant:
                if (is.element(interpMetrics[i, "metrics"], c("AbundanceDensity", "AbundanceDensityCorr")) &&
                    ! all(is.na(biologicalData$Weight)))
                {
                    ## Mean weight:
                    tmpBiologicalData <- subset(biologicalData,
                                                (is.element(LD_HaulNo,
                                                            stationsTmp$TRACKID) &
                                                 is.element(LD_SpeciesID,
                                                            interpMetrics[i, "species"])))


                    if (interpMetrics[i, "sizeClass"])
                    {
                        tmpBiologicalData <- subset(tmpBiologicalData,
                                                    Size >= interpMetrics[i, "min"] &
                                                    Size <= interpMetrics[i, "max"] &
                                                    ! is.na(Size))
                    }else{}

                    tmpBiolData <- cbind(tmpBiologicalData[ ,
                                                           ! is.element(colnames(tmpBiologicalData),
                                                                        colnames(stationCateg))],
                                         stationCateg[match(as.character(tmpBiologicalData[ , "LD_HaulNo"]),
                                                            as.character(stationCateg[ , "stationID"])),
                                                      ])


                    tail(tmpBiolData)
                    dim(tmpBiolData)

                    statWeight <- data.frame("Area" = gArea(SP,
                                                            byid = TRUE)[levels(stationsTmp$categDens)],
                                             "W"=tapply(tmpBiolData$Weight,
                                                        as.list(tmpBiolData[ , c("categDens"), drop=FALSE]),
                                                        mean,
                                                        na.rm=TRUE),
                                             "sdW"=tapply(tmpBiolData$Weight,
                                                          as.list(tmpBiolData[ , c("categDens"), drop=FALSE]),
                                                          sd,
                                                          na.rm=TRUE),
                                             "nbW"=tapply(tmpBiolData$Weight,
                                                          as.list(tmpBiolData[ , c("categDens"), drop=FALSE]),
                                                          function(x)sum( ! is.na(x))),
                                             "seW" = tapply(tmpBiolData$Weight,
                                                            as.list(tmpBiolData[ , c("categDens"), drop=FALSE]),
                                                            sd,
                                                            na.rm=TRUE) /
                                                 sqrt(tapply(tmpBiolData$Weight,
                                                             as.list(tmpBiolData[ , c("categDens"), drop=FALSE]),
                                                             function(x)sum( ! is.na(x)))))

                    ## Aggregated statistics:
                    totalW <- c("Area" = sum(statWeight[ , "Area"], na.rm = TRUE),
                                "W" = weighted.mean(x = statWeight[ , "W"],
                                                    w = statWeight[ , "Area"], na.rm=TRUE),
                                "dsW" = NA,
                                "nbW" = sum(statWeight[ , "nbW"], na.rm = TRUE),
                                ## sd(var_st) (i.e. SE) according to Cochran (1953) and Saville (1977):
                                "seW" = sqrt(sum(statWeight[ , "Area"]^2 * statWeight[ , "sdW"]^2 /
                                                 statWeight[ , "nbW"], na.rm = TRUE) /
                                             (sum(statWeight[ , "Area"], na.rm = TRUE) ^ 2)))

                    statWeight <- rbind(statWeight,
                                        Total = totalW)[ , -1] # No need anymore for replicating the Area.

                    ## CI on weight:
                    statWeight$wAlpha <- qnorm(1 - getOption("surveyStat.alpha") /
                                               2) * statWeight[ , "seW"]
                    statWeight$wAlpha[is.na(statWeight$wAlpha) & ! is.na(statWeight$W)] <- 0

                    ## Estimated biomass density:
                    stats <- cbind(statDens, statWeight,
                                   "BD" = statDens$D * statWeight$W,
                                   "seBD" = statDens$D * statWeight$W *
                                       sqrt(statDens$seD^2 / statDens$D^2 + statWeight$seW^2 / statWeight$W^2),
                                   "bdAlpha" = statDens$D * statWeight$W *
                                       sqrt(statDens$dAlpha^2 / statDens$D^2 + statWeight$wAlpha^2 / statWeight$W^2))

                    stats$BD[is.na(stats$BD) & stats$D == 0] <- 0
                    stats$seBD[is.na(stats$seBD) & stats$D == 0] <- 0
                    stats$bdAlpha[is.na(stats$bdAlpha) & stats$D == 0] <- 0


                    stats["Total", "BD"] <- weighted.mean(x = head(stats[ , "BD"], -1),
                                                          w = head(stats[ , "Area"], -1), na.rm=TRUE)

                    stats["Total", "seBD"] <- sqrt(sum(head(stats[ , "Area"]^2 * stats[ , "seBD"]^2, -1)) /
                                                   (sum(head(stats[!is.na(stats[ , "seBD"]) , "Area"],
                                                             -1), na.rm = TRUE) ^ 2))

                    stats["Total", "bdAlpha"] <- stats["Total", "seBD"] *
                        qnorm(1 - getOption("surveyStat.alpha") / 2)

                    ## Estimated biomasse:
                    stats <- cbind(stats,
                                   "B" =stats$BD  * stats[ , "Area"] * 1e-6, # In tons.
                                   "bAlpha" = stats$bdAlpha  * stats[ , "Area"] * 1e-6)


                }else{
                    ## if (is.element(interpMetrics[i, "metrics"], c("AbundanceDensity", "AbundanceDensityCorr")) &&
                    ##     all(is.na(biologicalData$Weight))) return(NULL)
                    stats <- statDens
                }

                ## Renaming the fields properly:
                Dname <- switch(interpMetrics[i, "metrics"],
                                "AbundanceDensity" =,
                                "AbundanceDensityCorr" = "Density",
                                "WeightDensity" =,
                                "WeightDensityCorr" = "BiomassDensity",
                                "UnknownDensity")

                AbName <- switch(interpMetrics[i, "metrics"],
                                 "AbundanceDensity" =,
                                 "AbundanceDensityCorr" = "Abundance",
                                 "WeightDensity" =,
                                 "WeightDensityCorr" = "Biomass",
                                 "UnknownDensity")

                colNames <- c(paste(c(D = "mean",
                                      sdD = "sd",
                                      nbD = "nb",
                                      seD = "se",
                                      dAlpha = "d"),
                                    Dname, sep=""),
                              paste(c(Ab = "",
                                      abAlpha = "d"),
                                    AbName, sep=""),
                              c(W = "meanWeight",
                                sdW = "sdWeight",
                                nbW = "nbWeight",
                                seW = "seWeight",
                                wAlpha = "dWeight",
                                BD = "meanBiomassDensity",
                                seBD = "seBiomassDensity",
                                bdAlpha = "dBiomassDensity",
                                B = "Biomass2",
                                bAlpha = "dBiomass2",
                                Area = "Area"))
                names(colNames)[1:7] <- c("D", "sdD", "nbD", "seD", "dAlpha", "Ab", "abAlpha")

                colnames(stats) <- colNames[colnames(stats)]

                return(stats)
            },
                SP = densInterpPolyFusion,
                stationData = stationData,
                simplify = FALSE)

invisible(sapply(names(stats),
                 function(i, stats)
             {
                 write.csv(stats[[i]],
                           file = file.path(ResultsPath, "Stratified_(deprecated)",
                                            paste("statsAbundance_", i, ".csv", sep = "")))

                 message("Statistics on abundance from the metric \"",
                         interpMetrics[i, "metrics"],
                         "\" and species \"",
                         interpMetrics[i, "species"],
                         "\"\n\texported in the file ",
                         paste("statsAbundance_", i, ".csv", sep = ""))
             },
                 stats = stats))

## densInterpPolyFusion[[2]]
## gArea(densInterpPolyFusion[[2]], byid = TRUE)
## names(densInterpPolyFusion)

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
