#-*- coding: latin-1 -*-

### File: 5_Raised_abundances_Bayes.R
### Time-stamp: <2018-11-01 12:13:18 yreecht>
###
### Created: 30/08/2016	15:00:04
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


stationsCat <- sapply(names(densInterpPolyFusion),
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
                                                                 levels=levels(SP$categ)))

                   ## Append categories to the temporary station data:
                   stationsTmp@data <- cbind(stationsTmp@data[ ,
                                                              !is.element(colnames(stationsTmp@data),
                                                                          colnames(stationCateg))],
                                             stationCateg)

                   head(stationsTmp@data)

                   return(stationsTmp)
               },
               SP = densInterpPolyFusion,
               stationData = stationData,
               simplify = FALSE)

## stationsCat[[3]]@data


biolCat <- sapply(names(densInterpPolyFusion),
                  function(i, SP, stationData)
           {
               SP <- SP[[i]]

               ## Density category of every station (represented by its centroid):
               categs <- gCovers(spgeom1 = SP,
                                 spgeom2 = gCentroid(spgeom = stationData, byid = TRUE),
                                 byid=TRUE)

               ## Working on only stations with a category (some may be excluded if you used a mask):
               idx <- apply(categs, 1, sum) == 1

               stationsTmp <- stationData[idx, , drop = FALSE]

               stationCateg <- data.frame("stationID"=stationsTmp$TRACKID,
                                          "categDens"=factor(colnames(categs)[apply(categs[idx, ,
                                                                                           drop = FALSE],
                                                                                    1, which)],
                                                             levels=levels(SP$categ)))

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

               return(tmpBiolData)
           },
           SP = densInterpPolyFusion,
           stationData = stationData,
           simplify = FALSE)

## debuggingState(on=FALSE)
biomassCat <- sapply(X = names(densInterpPolyFusion),
                     FUN = function(i, stationsCat, biolCat, metrics)
              {
                  spC <- metrics[i, "species"]
                  statC <- stationsCat[[i]]@data[ , c(i, "stationID", "categDens")]
                  biolC <- subset(biolCat[[i]],
                                  LD_SpeciesID == spC)[ , c("Weight", "stationID", "categDens")]

                  mType <- metrics[i, "metrics"]

                  if (is.element(mType,
                                 c("AbundanceDensity", "AbundanceDensityCorr")))
                  {
                      if (all(is.na(biolC$Weight))) return(NULL)

                      meanweight <- by(biolC,
                                       as.list(biolC[ , c("categDens"), drop = FALSE]),
                                       function(x)
                                    {
                                        data.frame("categDens" = as.character(unique(x$categDens)),
                                                   "stationID" = tapply(as.character(x$stationID),
                                                                        as.character(x$stationID),
                                                                        unique),
                                                   "mWeight" = tapply(x$Weight,
                                                                      as.character(x$stationID),
                                                                      mean, na.rm = TRUE))
                                    }, simplify = FALSE)

                      meanweight <- do.call(rbind, meanweight)

                      biomDens <- merge(statC, meanweight, all = TRUE)
                      biomDens$biomDensity <- biomDens[ , i] * biomDens$mWeight / 1000 # Kgs

                      biomDens$biomDensity[is.na(biomDens$biomDensity) &
                                           biomDens[ , i] == 0] <- 0

                      if (any(is.na(biomDens$biomDensity)))
                      {
                          ##browser()
                          warning("## metrics \"", i, "\": possible data inconsistencies:",
                                  "\n  ##   density > 0 but no individual measured (weight/size+WL)",
                                  "\n  ##   in stations ",
                                  paste(as.character(biomDens$stationID[is.na(biomDens$biomDensity)]),
                                        collapse = ", "))
                      }
                  }else{
                      biomDens <- statC
                      biomDens$biomDensity <- biomDens[ , i]
                      biomDens[ , i] <- NULL
                  }

                  head(biomDens)
                  dim(biomDens)

                  return(biomDens)
              },
                  stationsCat = stationsCat, biolCat = biolCat, metrics = metrics,
                  simplify=FALSE, USE.NAMES=TRUE)

biomassCat <- biomassCat[ ! sapply(biomassCat, is.null)]


statBayes <- sapply(X = names(biomassCat),
                    FUN = function(i, biomassCat, densInterpPolyFusion)
             {
                 message("\n## ##################################################",
                         "\n## ", i ,
                         ":\n")

                 gc()

                 biomassData <- biomassCat[[i]]

                 areas <- densInterpPolyFusion[[i]]@data

                 listBdens <- tapply(biomassData$biomDensity,
                                     biomassData$categDens,
                                     function(x) x[! is.na(x)], simplify = FALSE)

                 ## Max number of stations by strata:
                 maxStat <- max(sapply(listBdens, length))

                 ## Matrix of biomass densities:
                 matBdens <- do.call(cbind,
                                     sapply(listBdens,
                                            function(x, maxL)
                                     {
                                         if (maxL > length(x)) x <- c(x, rep(NA, maxL - length(x)))
                                         return(x)
                                     },
                                            maxL = maxStat,
                                            simplify = FALSE))

                 ## Number of stations by strata:
                 Nbc <- apply(matBdens, 2, function(x) sum(! is.na(x)))
                 head(matBdens, 29)

                 ## Surface area
                 areascl <- areas$area

                 ## require(rjags)
                 ## require(R2jags)

                 ## if (i == "AbundanceDensity_78_Inf_Ostrea_edulis") browser()

                 dataList <- list(matBdens = matBdens +
                                      ifelse(any(na.omit(as.vector(matBdens)) == 0),
                                             10^(floor(log10(min(matBdens[matBdens != 0], na.rm = TRUE))) - 2),
                                             0),
                                  ## minDens = apply(matBdens, 2, min, na.rm = TRUE) +
                                  ##     ifelse(any(na.omit(as.vector(matBdens)) == 0),
                                  ##            0.001, 0),
                                  ## maxDens = apply(matBdens, 2, max, na.rm = TRUE) +
                                  ##     ifelse(any(na.omit(as.vector(matBdens)) == 0),
                                  ##            0.001, 0),
                                  meanDens = apply(matBdens, 2, mean, na.rm = TRUE) +
                                      ifelse(any(na.omit(as.vector(matBdens)) == 0),
                                             10^(floor(log10(min(matBdens[matBdens != 0], na.rm = TRUE))) - 2),
                                             0),
                                  Nbc = Nbc,
                                  areascl = areascl)

                 sd0 <- which(apply(dataList$matBdens, 2, sd, na.rm = TRUE) == 0)

                 if (length(sd0))
                 {
                     matDensTmp <- dataList$matBdens

                     matDensTmp[ , sd0] <- matDensTmp[ , sd0, drop = FALSE] +
                         sapply(sd0,
                                function(x, data)
                         {
                             res <- rnorm(n = nrow(data),
                                          mean = 0,
                                          sd = data[1, x] * 1e-03)
                         },
                         data = matDensTmp)

                     dataList$matBdens <- matDensTmp
                 }

                 funcInit <- function()
                 {
                     res <- list("mu" = runif(n = ncol(dataList$matBdens),
                                              min = apply(dataList$matBdens, 2, min, na.rm = TRUE) * 1,
                                              max = apply(dataList$matBdens, 2, max, na.rm = TRUE) * 1),
                                 "sigma" = runif(n = ncol(dataList$matBdens),
                                                 min = apply(dataList$matBdens, 2, sd, na.rm = TRUE) * 0.8,
                                                 max = apply(dataList$matBdens, 2, sd, na.rm = TRUE) * 1.2))

                     res$sigma[is.na(res$sigma)] <- mean(res$sigma, na.rm = TRUE)
                     ## res$mu[res$sigma == 0] <- runif(sum(res$sigma == 0),
                     ##                                 res$mu[res$sigma == 0] * 0.8,
                     ##                                 res$mu[res$sigma == 0] * 1.2)

                     ## res$sigma[res$sigma == 0] <- runif(sum(res$sigma == 0),
                     ##                                    10^(floor(log10(min(res$sigma[res$sigma != 0]))) - 2) * 0.8,
                     ##                                    10^(floor(log10(min(res$sigma[res$sigma != 0]))) - 2) * 1.2)

                     return(res)
                 }

                 modelScripts <- c("LN" = file.path(scriptDir,
                                                    "5_Biomass_LN.jag"),
                                   "G" = file.path(scriptDir,
                                                   "5_Biomass_G.jag"),
                                   "N" = file.path(scriptDir,
                                                   "5_Biomass_N.jag"))[getOption("surveyBayes.models")]

                 models <- sapply(modelScripts,
                                  function(x, ...)
                                  {
                                      message("\n## Model ", x, ":")
                                      jags(model.file = x, ...)
                                  },
                                  data = dataList,
                                  n.chains = getOption("surveyBayes.n.chains"),
                                  n.iter = getOption("surveyBayes.n.iter"),
                                  n.burnin = getOption("surveyBayes.n.burnin"),
                                  inits = funcInit,
                                  n.thin = 1,
                                  parameters.to.save = c("mu", "sigma", "Bcl",
                                                         "BDtot", "Btot"),
                                  simplify = FALSE)


                 dic <- sapply(models,
                               function(x) x$BUGSoutput$DIC)

                 samples <- sapply(models,
                                   function(x) as.mcmc.list(x$BUGSoutput),
                                   simplify = FALSE)

                 bestModel <- which.min(dic)
                 message("\n## ", i , " best model for biomass assessment is: ",
                         c("Log Normal", "Normal", "Gamma")[bestModel],
                         "\n\n")

                 res <- list("metrics" = i,
                             "models" = models,
                             "samples" = samples,
                             "dic" = dic,
                             "bestModel" = bestModel,
                             "Nstations" = Nbc,
                             "selectedModel" = models[[bestModel]],
                             "selectedSample" = samples[[bestModel]],
                             "selectedDic" = dic[[bestModel]])

                 return(res)
             },
             biomassCat = biomassCat, densInterpPolyFusion = densInterpPolyFusion,
             simplify=FALSE, USE.NAMES=TRUE)

## #################################################################################################
##

distributionSummary <- sapply(statBayes,
       function(x)
       {##browser()
           res <- MCMCdistribution.summary(sampleVec=x$selectedSample[ , "Btot"],
                                           credMass=1 - getOption("surveyStat.alpha"),
                                           central=c("mean", "median"), print=FALSE)

           ## Convergence and CI diagnostics:
           Neff <- unname(coda::effectiveSize(as.mcmc(do.call(rbind,
                                                              x$selectedSample))[ , "Btot"]))

           N <- length(do.call(rbind, x$selectedSample)[ , "Btot"])

           ##
           Nmin.i <- raftery.diag(data = as.mcmc(do.call(rbind, x$selectedSample))[ , "Btot"],
                                  q = getOption("surveyStat.alpha") / 2)$resmatrix[ , "N"]

           Nmin.u <- raftery.diag(data = as.mcmc(do.call(rbind, x$selectedSample))[ , "Btot"],
                                  q = 1 - getOption("surveyStat.alpha") / 2)$resmatrix[ , "N"]

           Nmin <- max(Nmin.i, Nmin.u)

           PSR <- gelman.diag(x$selectedSample[ , "Btot"],
                              autoburnin = FALSE)$psrf[1, ]
           names(PSR) <- c("PSR", "PSR.CI.up")


           cat("\n## Biomass (t) from \"", x$metrics, "\", model ", names(x$bestModel), ":\n", sep = "")
           print(cbind(matrix(round(res,
                                    ifelse(ceiling(log10(res)) > 3, 0, 3 - ceiling(log10(res)))),
                              nrow = 1,
                              dimnames = list(" ", names(res))),
                       data.frame("convergence" = ifelse(unname(PSR["PSR.CI.up"]) < 2,
                                                         "OK", "Check"),
                                  "CI.diag" = ifelse(N >= Nmin,
                                                     "OK", "Increase N"))))

           res <- c(res, "model" = names(x$bestModel),
                    round(PSR["PSR.CI.up"], 2),
                    "convergence" = ifelse(unname(PSR["PSR.CI.up"]) < 2,
                                           "OK", "Check"),
                    "N" = N, "N.eff" = round(Neff),
                    "N.req" = Nmin,
                    "CI.diag" = ifelse(N >= Nmin,
                                       "OK", "Increase N"))

           return(res)
       })

colnames(distributionSummary) <- sub("^Abundance", "Abundance+L-W",
                                     sub("^([[:alpha:]]+)Density", "\\1",
                                         sub("Corr_", "_(corrected)_", colnames(distributionSummary))))

write.csv(t(distributionSummary),
          file = file.path(ResultsPath, "Bayesian",
                           "TotalBiomass_distributions_summary.csv"))

## #################################################################################################
## Detailed stat:

distSummTot <- sapply(statBayes,
                      function(x, SP)
               {
                   tmpBreaksNames <- row.names(SP[[x$metrics]])

                   ## Organising the distributions information in a data.frame:
                   res <- sapply(seq_along(tmpBreaksNames),
                                 function(i, x, nbreaks){
                                     ## Mu and Bcl names without "[]" if only one stratum:
                                     mu <- ifelse(nbreaks > 1,
                                                  paste0("mu[", i , "]"),
                                                  "mu")

                                     Bcl <- ifelse(nbreaks > 1,
                                                   paste0("Bcl[", i , "]"),
                                                   "Bcl")

                                     ## Summary of the biomass density distribution:
                                     BD <- MCMCdistribution.summary(sampleVec = x$selectedSample[ , mu],
                                                                    credMass = 1 - getOption("surveyStat.alpha"),
                                                                    central = c("mean"),
                                                                    print=FALSE)

                                     names(BD) <- c("meanBiomassDensity",
                                                    paste0("CI",
                                                           round(100 * (1 - getOption("surveyStat.alpha"))),
                                                           "% ", c("low", "upp"), " BiomassDensity"))

                                     ## Summary of the biomass distribution:
                                     B <- MCMCdistribution.summary(sampleVec = x$selectedSample[ , Bcl],
                                                                   credMass = 1 - getOption("surveyStat.alpha"),
                                                                   central = c("mean"),
                                                                   print=FALSE)

                                     names(B) <- c("meanBiomass (t)",
                                                   paste0("CI",
                                                          round(100 * (1 - getOption("surveyStat.alpha"))),
                                                          "% ", c("low", "upp"), " Biomass"))

                                     resi <- data.frame(bestModel = "")

                                     resi <- cbind(resi,
                                                   matrix(BD, nrow = 1, dimnames = list(NULL , names(BD))),
                                                   matrix(B, nrow = 1, dimnames = list(NULL , names(B))))
                                 },
                                 x = x, nbreaks = length(tmpBreaksNames),
                                 simplify = FALSE)

                   res <- do.call(rbind, res)

                   res <- cbind(Area = gArea(spgeom = SP[[x$metrics]], byid=TRUE),
                                StationNumber = x$N,
                                res)

                   totRes <- data.frame(gArea(spgeom = SP[[x$metrics]], byid=FALSE),
                                        sum(x$N),
                                        names(x$bestModel),
                                        matrix(MCMCdistribution.summary(sampleVec = x$selectedSample[ , "BDtot"],
                                                                        credMass = 1 - getOption("surveyStat.alpha"),
                                                                        central = c("mean"),
                                                                        print=FALSE),
                                               nrow = 1),
                                        matrix(MCMCdistribution.summary(sampleVec = x$selectedSample[ , "Btot"],
                                                                        credMass = 1 - getOption("surveyStat.alpha"),
                                                                        central = c("mean"), print=FALSE),
                                               nrow = 1))
                   colnames(totRes) <- colnames(res)

                   res <- rbind(res,
                                Total = totRes)

                   write.csv(res,
                             file = file.path(ResultsPath, "Bayesian",
                                              paste(x$metrics, "_Biomass_assessment.csv", sep = "")))

                   return(res)
               }, SP = densInterpPolyFusion, simplify = FALSE)


## #################################################################################################
## Plots :
invisible(sapply(statBayes,
                 function(x)
{
    ## browser()
    fileNameTrace <- paste("TracePlot", x$metrics, names(x$bestModel), sep = "_")
    fileNameDensity <- paste("DensityPlot", x$metrics, names(x$bestModel), sep = "_")

    devNoX11 <- getOption("surveyPlot.dev")
    devNoX11 <- devNoX11[ ! is.element(devNoX11, c("X11", "x11"))]
    ## Graphical representation:
    for (devType in devNoX11)
    {
        gc()

        if (devType %in% c("pdf", "jpg", "jpeg", "png"))
        {
            width <- c(png = 1000, jpg = 1000, jpeg = 1000, pdf = 6)
            dev <- openDev(device = devType,
                           directory = file.path(ResultsPath, "Bayesian/Diagnostics"),
                           filename = fileNameTrace,
                           width = width,
                           height = width ,
                           pointsize = c(png = 24, jpg = 24, jpeg = 24, pdf = 10),
                           verbose=FALSE)
        }else{
            warning("Device \"", devType, "\" not supported")
            next()
        }

        params <- grep("^(mu|sigma)", colnames(x$selectedSample[[1]]), value = TRUE)

        par(mfrow = c(2, 2))

        traceplot(x$selectedSample[ , params])

        dev.off()

        if (devType %in% c("pdf", "jpg", "jpeg", "png"))
        {
            width <- c(png = 1000, jpg = 1000, jpeg = 1000, pdf = 6)
            dev <- openDev(device = devType,
                           directory = file.path(ResultsPath, "Bayesian/Diagnostics"),
                           filename = fileNameDensity,
                           width = width,
                           height = width ,
                           pointsize = c(png = 24, jpg = 24, jpeg = 24, pdf = 10),
                           verbose=FALSE)
        }else{
            warning("Device \"", devType, "\" not supported")
            next()
        }

        params <- grep("^B[^D]", colnames(x$selectedSample[[1]]), value = TRUE)

        par(mfrow = c(2, 2))

        densplot.hdi(x$selectedSample[ , params], show.obs=FALSE, trimxThr = 1.5)


        ## str(test)
        ## hdi <- MCMCdistribution.summary(sampleVec=x$selectedSample[ , "Btot"],
        ##                                 credMass=0.95,
        ##                                 central=c("none"), print=FALSE)


        ## densplot.shade(x$selectedSample[ , "Btot"], hdi = hdi, show.obs=FALSE)
        dev.off()

    }
}))


## par(mfcol = c(3, 5))
## test <- densplot.hdi(statBayes[[1]]$selectedSample)

## X11()
## test2 <- densplot.hdi(statBayes[[1]]$selectedSample[ , "Btot"], mai = "test", xlab = "bla")
## attr(test2, "HDIs")

## dim(statBayes[[1]]$selectedSample[ , "Btot"])

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
