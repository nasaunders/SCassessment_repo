#-*- coding: utf-8 -*-

### File: 1_test_krigging_bt.R
### Time-stamp: <2019-04-10 18:17:23 yreecht>
###
### Created: 02/08/2017	12:41:51
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

if (grepl("windows", .Platform$OS.type))
{
    memory.limit(size = ifelse(.Platform$r_arch == "i386",
                               4095, 10240))
}

boxcox.transf <- function(x, lambda = 0, offset = 0,
                          backtransformation = FALSE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 10 Jul 2018, 16:13

    if (isTRUE(backtransformation))
    {
        if (lambda == 0)
        {
            res <- exp(x)
        }else{
            res <- (lambda * x + 1) ^ (1 / lambda)
        }

        res <- res - offset
    }else{
        x <- x + offset

        stopifnot(all(x[ ! is.na(x)] > 0))

        if (lambda == 0)
        {
            res <- log(x)
        }else{
            res <- (x ^ lambda - 1) / lambda
        }
    }
    return(res)
}

boxcox.mean.backtransf <- function(mean, var, lambda = 0, offset = 0)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 10 Jul 2018, 16:41

    if (lambda == 0)
    {
        mean <- exp(mean) * ( 1 + var / 2)
    }else{
        mean <- (lambda * mean + 1) ^ (1 / lambda) *
            (1 + (var * (1 - lambda)) / (2 * (lambda * mean + 1) ^ 2))
    }

    mean <- mean - offset
    return(mean)
}

source(file.path(scriptDir,
                 "6_Kriging_biomass_functions.R"),
       encoding="latin1")


resGeostats <- list()

## Raster grid from the concave hull:
concaveHull <- concaveHull.survey(bufferedTracks = gBuffer(stationData,
                                                           byid=TRUE,
                                                           width=getOption("surveyIDW.bufferWidth"),
                                                           capStyle="ROUND"),
                                  SelectionMask = mask,
                                  unsuitableArea = unsuitableGround,
                                  suitableArea = suitableGround,
                                  suitableAreaForce = getOption("surveyBedForce"),
                                  alpha = getOption("surveyIDW.alpha"),
                                  bufferWidth =  getOption("surveyIDW.bufferWidth"))

if (is.null(getOption("surveyGeostat.predict.raster")) ||
    ! exists(getOption("surveyGeostat.predict.raster")) ||
    ! class(get(getOption("surveyGeostat.predict.raster"))) %in% c("SpatialPixelsDataFrame"))
{
    data.pix <- zone2SpatialPixelsDataFrame(spgeom = concaveHull,
                                            cellsize = getOption("surveyIDW.cellsize"))
}else{
    data.pix <- get(getOption("surveyGeostat.predict.raster"))
}


## Zone as a possible covariate if mask defined:
if (! is.null(mask))
{
    data.pix[["zone"]] <- factor(over(data.pix, mask)[ , 1])
}


if (is.null(getOption("surveyGeostat.vg.cutoff")))
{
    options(surveyGeostat.vg.cutoff = round(as.vector(dist(t(bbox(data.pix)))) / 2))
}else{}

if (is.null(getOption("surveyGeostat.vg.width")))
{
    options(surveyGeostat.vg.width = round(getOption("surveyGeostat.vg.cutoff") / 15))
}else{}

if (is.null(getOption("krigging.model")))
{
    options(krigging.model = "~ 1")
}


## summary(data.pix)
## proj4string(data.pix)

## %in% densUnits
## metricsName <- row.names(metrics[metrics$metrics %in% densUnits, ])[1]
## metricsName <- row.names(metrics[metrics$metrics %in% "AbundanceDensity", ])[2]
for (metricsName in row.names(metrics[metrics$metrics %in% densUnits, ]))
{
    message("\n## ", paste0(rep("#", 80), collapse = ""))
    message("## ", metricsName, ":")

    meanWeightName <- sub("^AbundanceDensity", "meanWeight", metricsName)
    densBiomassName <- sub("^AbundanceDensity", "AbundanceDensityLW", metricsName)

    ## Preparation of data (calculate the biomass density if needed):
    if (metrics[metricsName, "metrics"] %in% c("AbundanceDensity", "AbundanceDensityCorr"))
    {
        ## Mean individual weight for the species and size class (if relevant):
        if (metrics[metricsName, "sizeClass"])
        {
            meanWeight <- with(subset(biologicalData,
                                      Size < metrics[metricsName, "max"] &
                                      Size >= metrics[metricsName, "min"] &
                                      LD_SpeciesID ==  metrics[metricsName, "species"]),
                               tapply(Weight, LD_HaulNo, mean, na.rm = TRUE))
        }else{
            meanWeight <- with(subset(biologicalData,
                                      LD_SpeciesID ==  metrics[metricsName, "species"]),
                               tapply(Weight, LD_HaulNo, mean, na.rm = TRUE))
        }

        ## Checking if weights available for calculating biomass from density:
        if (all(is.na(meanWeight[stationData@data[match(names(meanWeight),
                                                        stationData@data[ , "TRACKID"]),
                                                  metricsName] > 0])))
        {
            message("\n\t## No weight data: skipping the metrics.")
            next()
        }else{
            if (any(is.na(meanWeight[stationData@data[match(names(meanWeight),
                                                            stationData@data[ , "TRACKID"]),
                                                      metricsName] > 0])))
            {
                warning("Some missing mean weight by station")
            }
        }

        ## Saving mean weight for the size class (for convenience):
        stationData@data[ , meanWeightName] <- meanWeight[match(toupper(stationData$TRACKID),
                                                                toupper(names(meanWeight)))]

        ## Mean weight when no individual measured assumed to be the overall mean for the size class
        ##   ... might be idealy replaced by idw on mean weight.
        stationData@data[(! is.na(stationData@data[ , metricsName]) &
                          stationData@data[ , metricsName] > 0 &
                          is.na(stationData@data[ , meanWeightName])),
                         meanWeightName] <- mean(stationData@data[ , meanWeightName], na.rm = TRUE)

        ## Biomass density, in *Kg* (from g):
        stationData@data[ , densBiomassName] <-
            stationData@data[ , metricsName] * stationData@data[ , meanWeightName] * 1e-3

        ## Some biomass are NA as no individual for estimating mean weight... but should be zero:
        stationData@data[is.na(stationData@data[ , densBiomassName]) &
                         ! is.na(stationData@data[ , metricsName]) &
                         stationData@data[ , metricsName] == 0,
                         densBiomassName] <- 0
    }else{
        ## densBiomassName <- metricsName
    }

    ## ##################################################
    ## Geostat:

    require(gstat)


    ## Need for point data:
    if ("SpatialLinesDataFrame" %in% class(stationData))
    {
        terms <- as.character(as.list(as.formula(getOption("krigging.model"))))[-1]
        terms <- unlist(strsplit(terms, split = "[[:blank:]]*[*+-/][[:blank:]]*"))

        addColZone <- names(stationData@data)[names(stationData@data) %in%
                                              c("zone",
                                                terms[suppressWarnings(is.na(as.numeric(terms)))])]

        ## names(stationData@data)
        stationPosData <- SpatialPointsDataFrame(coords = gCentroid(stationData, byid = TRUE),
                                                 data = stationData@data[ , c("TRACKID",
                                                                              "TIMESTART",
                                                                              "length",
                                                                              "effort",
                                                                              metricsName,
                                                                              meanWeightName,
                                                                              densBiomassName,
                                                                              addColZone)],
                                                 proj4string = CRS(proj4string(stationData)),
                                                 match.ID = TRUE)
    }else{
        stationPosData <- stationData
    }


    ## Choice of modelling approach:
    if (sum(stationData@data[ , densBiomassName] == 0, na.rm = TRUE) /
        sum( ! is.na(stationData@data[ , densBiomassName])) >= getOption("surveyGeostat.propZero.max"))
    {
        ## ####################################################################################################
        ## Two-steps model:
        ##   1) presence/absence.
        ##   2) abundances for presence only.

        formulaPres <- formula(paste0("pres ", getOption("krigging.model")))

        message("\n# Two-step model (modelling presence)!")

        ## #################
        ## Presence/absence:
        stationPosData$pres <- stationPosData@data[ , metricsName] > 0 # metricsName instead of
                                        # densBiomassName to avoid intermediate issues with
                                        # biomass calculation.
    }else{
        ## ####################################################################################################
        ## One-step model:

        formulaPres <- NULL

        message("\n# One-step model (abundance only)!")
    }

    ## Estimate the cst to add in boxcox transformation:
    nbZero <- sum(stationData@data[ , densBiomassName] == 0, na.rm = TRUE)

    if (as.logical(nbZero) && is.null(formulaPres) || isTRUE(getOption("surveyGeostat.boxcox.cst")))
    {
        if (!is.null(getOption("surveyGeostat.boxcox.cst")) &&
            is.numeric(getOption("surveyGeostat.boxcox.cst")))
        {
            cst <- getOption("surveyGeostat.boxcox.cst")
        }else{
            cst <- min(stationData@data[stationData@data[ , densBiomassName] != 0,
                                        densBiomassName], na.rm = TRUE) * 0.1
        }
    }else{
        cst <- 0
    }

    formulaDens <- formula(paste0("boxcox.transf(x = ", densBiomassName,
                                  ", lambda = lambda, offset = cst) ", getOption("krigging.model")))

    ## BoxCox transformation:
    lambda <- boxCoxDiag(stationPosData = stationPosData,
                         varName = densBiomassName,
                         presenceOnly = TRUE,
                         ResultsPath = ResultsPath,
                         offset.bc = cst)

    ## Transformation and removing non-finit values:
    idxPres <- (! is.na(stationPosData@data[ , densBiomassName]) &
                stationPosData@data[ , densBiomassName] > ifelse(is.null(formulaPres), -1, 0))

    suppressWarnings(stationPosData@data$densTransfo[idxPres] <-
                         boxcox.transf(x = stationPosData@data[idxPres , densBiomassName],
                                       lambda = lambda, offset = cst))
    stationPosData@data$densTransfo[! is.finite(stationPosData@data$densTransfo)] <- NA

    ## head( stationPosData@data)

    ## ##################################################
    ## Variogram fitting:
    varioFit <- tryCatch(surveyKrigingVariograms(formulaPres = formulaPres,
                                                 formulaDens = formulaDens,
                                                 varName = densBiomassName,
                                                 stationPosData = stationPosData,
                                                 autoBreak = getOption("surveyGeostat.vg.autobreaks"),
                                                 cutoff = getOption("surveyGeostat.vg.cutoff"),
                                                 width = getOption("surveyGeostat.vg.width"),
                                                 use.lik = getOption("surveyGeostat.vg.likelihood"),
                                                 lik.method = getOption("surveyGeostat.vg.lik.method"),
                                                 krigging.model = getOption("krigging.model"),
                                                 lambda = lambda, offset.bc = cst,
                                                 offset.name = "cst",
                                                 verbose=TRUE),
                         error = function(e)
                         {
                             warning("Could not fit the variogram for variable \"", densBiomassName,
                                     "\"\n\n", e)
                             return(NULL)
                         })

    names(varioFit)

    ## ##################################################
    ## Ordinary/Universal kriging (two-steps model):
    krigRes <- tryCatch(surveyKriging(formulaDens = formulaDens, vmfDens = varioFit$vmfDensAuto,
                                      formulaPres = formulaPres, vifPres = varioFit$vifPresAuto,
                                      stationPosData = stationPosData,
                                      data.pix = data.pix,
                                      varName = densBiomassName,
                                      correctPres = TRUE,
                                      lambda = lambda, offset.bc = cst,
                                      offset.name = "cst"),
                        error = function(e)
                        {
                            warning("Could not run kriging for variable \"", densBiomassName,
                                    "\"\n\n", e)
                            return(NULL)
                        })

    ## names(krigRes)

    ## ###################################################
    ## Conditional Gaussian simulations (two-steps model):
    krigResSim <- tryCatch(
                           surveyGaussianSimu(formulaDens = formulaDens, vmfDens = varioFit$vmfDensAuto,
                                              formulaPres = formulaPres, vifPres = varioFit$vifPresAuto,
                                              stationPosData = stationPosData,
                                              data.pix = data.pix,
                                              varName = densBiomassName,
                                              lambda = lambda, offset.bc = cst,
                                              offset.name = "cst",
                                              nmax = getOption("surveyGeostat.CGS.nmax"),
                                              nsim = getOption("surveyGeostat.CGS.nsim")),
                           error = function(e)
                           {
                               warning("Could not run conditional gaussian simulation for variable \"",
                                       densBiomassName,
                                       "\"\n\n", e)
                               return(NULL)
                           })

    ## names(krigResSim)

    resGeostats[[metricsName]] <- c(varioFit,
                                    krigRes,
                                    krigResSim,
                                    list(metrics = densBiomassName,
                                         models = list("presence" = formulaPres,
                                                       "density" = formulaDens),
                                         stationPosData = stationPosData,
                                         boxcoxParams = c(lambda = lambda, offset.bc = cst)))
}

## #################################################################################################
## Graphics and statistics:
## names(resGeostats)
resGeostats <-
    sapply(resGeostats,
           function(x)
           {
               if (!is.null(x$krigingPres))
               {
                   ## ##############################################################################
                   ## Graphs:
                   PresPlotVar <- spplot(x$krigingPres, "var1.var", asp=1, col.regions=cm.colors(64),
                                         main = "Presence/Absence variability")
                   PresPlotPred <- spplot(x$krigingPres, "var1.pred", asp=1, col.regions=bpy.colors(64),
                                          main = "Presence probability")

                   biomDplotVar <- spplot(x$krigingDens, "var1.var", asp=1, col.regions=cm.colors(64),
                                          main = "BoxCox-transformed biomass density variability | presence")
                   biomDplotPred <- spplot(x$krigingDens, "var1.pred", asp=1, col.regions=bpy.colors(64),
                                           main = "BoxCox-transformed biomass density | presence")

                   biomDplotBT <- spplot(x$krigingDens, "DensPos", asp=1, col.regions=bpy.colors(64),
                                         main = expression("Biomass density | presence "~(kg.m^{-2})))

                   biomDplotComb <- spplot(x$krigingDens, "Dens", asp=1, col.regions=bpy.colors(64),
                                           main = expression("Biomass density"~(kg.m^{-2})))

                   graphs <- list("presence.var" = PresPlotVar, "presence.pred" = PresPlotPred,
                                  "density.var" = biomDplotVar, "density.pred" =  biomDplotPred,
                                  "density|pres.back.transformed" = biomDplotBT,
                                  "density.back.transformed" = biomDplotComb)
               }else{
                   if (!is.null(x$krigingDens))
                   {
                       biomDplot1 <- spplot(x$krigingDens, "var1.var", asp=1, col.regions=cm.colors(64),
                                            main = "BoxCox-transformed biomass density variability")
                       biomDplot2 <- spplot(x$krigingDens, "var1.pred", asp=1, col.regions=bpy.colors(64),
                                            main = "BoxCox-transformed biomass density")

                       biomDplot3 <- spplot(x$krigingDens, "Dens", asp=1, col.regions=bpy.colors(64),
                                            main = expression("Biomass density "~(kg.m^{-2})))

                       graphs  <-  list("density.pred" = biomDplot2,
                                        "density.var" = biomDplot1,
                                        "density.back.transformed" = biomDplot3)
                   }else{
                       graphs <- list(NULL)
                   }
               }

               ## if (metricsName == "AbundanceDensity_Ensis_magnus") browser()

               meanB <- if (!is.null(x$krigingDens))
                        {
                            sum(x$krigingDens[["Dens"]], na.rm = TRUE) * prod(data.pix@grid@cellsize) / 1e3
                        }else{
                            NA
                        }

               if (! is.null(x$kDensSim))
               {
                   simuBiomass <- x$abdSimu

                   ## message("\n# Result for ", metricsName, ":")
                   resHDI <- MCMCdistribution.summary(sampleVec=simuBiomass,
                                                      credMass=0.95,
                                                      central=c("mean", "median"), print=FALSE)
               }else{
                   simuBiomass <- NA
                   resHDI <- c("mean" = NA,
                               "median" = NA,
                               "95% HDI inf" = NA,
                               "95% HDI sup" = NA)
               }


               res <- c(x,
                        list(meanB = meanB,
                             graphs = graphs,
                             resHDI = resHDI))

               return(res)
           }, simplify = FALSE)

## ####################################################################################################
## Simulations summary:
distributionSummaryGS <- sapply(resGeostats,
                                function(x)
                         {
                             res <- x$resHDI
                             names(res)[1:2] <- paste0("simu_", names(res)[1:2])

                             res <- c(krig_mean = x$meanB,
                                      "delta%" = unname(round(100 * (res[1] - x$meanB) / x$meanB,
                                                              2)),
                                      res)

                             res <- res[c(1, 3:2, 4:length(res))]

                             cat("\n## Biomass (t) from \"",
                                 ## sub("^AbundanceDensityLW", "AbundanceDensity+WL",
                                 ##     sub("Corr_", "-Corr_",
                                 x$metrics,
                                 ## )),
                                 ":\n", sep = "")
                             print(round(res,
                                         ifelse(test = ceiling(log10(abs(res))) > 3,
                                                yes = 0,
                                                no =ifelse(3 - ceiling(log10(abs(res))) > 2,
                                                           2,
                                                           3 - ceiling(log10(abs(res)))))))

                             res <- c(res,
                                      "model" = ifelse("viPres" %in% names(x),
                                                       "Two-steps", "One-step"))

                             return(res)
                         })

colnames(distributionSummaryGS) <- sub("^Abundance", "Abundance+L-W",
                                       sub("^([[:alpha:]]+)Density", "\\1",
                                           sub("Corr_", "_(corrected)_", colnames(distributionSummaryGS))))

write.csv(t(distributionSummaryGS),
          file = file.path(ResultsPath, "Geostatistics",
                           paste0("TotalBiomass_distributions_summaryGS_",
                                  "(", gsub("[[:blank:]]", "", getOption("krigging.model")),
                                  ").csv")))


## ####################################################################################################
## Graphical results:

## Graphics:
sapply(resGeostats,
       function(x)
       {
           if (! is.null(x$krigingDens))
           {
               ## Maps (several formats possible):
               for (devType in getOption("surveyPlot.dev"))
               {

                   if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
                   {
                       width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6) * 5 / 6
                       dev <- openDev(device = devType,
                                      directory = file.path(ResultsPath,
                                                            "Geostatistics"),
                                      filename = paste0("BiomassMap_from_", x$metrics,
                                                        "(",
                                                        gsub("[[:blank:]]", "",
                                                             getOption("krigging.model")),
                                                        ")"),
                                      width = width,
                                      height = width *
                                          diff(bboxMapZone[2, ]) / diff(bboxMapZone[1, ]),
                                      pointsize = c(png = 20, jpg = 20, jpeg = 20, X11 = 9, pdf = 9),
                                      counter = FALSE,
                                      verbose=FALSE)
                   }else{
                       warning("Device \"", devType, "\" not supported")
                       next()
                   }

                   plotDensityMap(raster=x$krigingDens, variable="Dens",
                                  palpers=colorRampPalette(rev(brewer.pal(9, "RdYlGn"))),
                                  legendTitle=expression("Biomass density\n ", (kg.m^{-2})))

                   if (devType != "X11") dev.off()
               }
           }
       })

resGeostats[[1]]$vmfDensAuto
resGeostats[[1]]$vifPresAuto

## Detailed graphs for diagnostics:
sapply(resGeostats,
       function(x)
       {
           ## Variogram (empirical + fit):
           png(file = file.path(ResultsPath,
                                "Geostatistics/Diagnostics",
                                paste("Variograms_", x$metrics,
                                      "(", gsub("[[:blank:]]", "", getOption("krigging.model")), ")",
                                      "_%01d.png", sep = "")),
               width = 800)

           if (is.null(x$models[["presence"]])) ## ("vLB" %in% names(x))
           {
               print(plot(x$vmDens, model = x$vmfDensAuto, plot.numbers = TRUE,
                          main = paste("Variogram for log(biomass) from \"",
                                       x$metrics, "\"",
                                       sep = ""),
                          ylim = c(0, 1.06 * max(c(x$vmDens$gamma, sum(x$vmfDensAuto$psill))))))
           }else{

               print(plot(x$viPres, model = x$vifPresAuto, plot.numbers = TRUE,
                          main = paste("Variogram for presence from \"",
                                       x$metrics, "\"",
                                       sep = ""),
                          ylim = c(0, 1.06 * max(c(x$viPres$gamma, sum(x$vifPresAuto$psill))))))

               print(plot(x$vmDens, model = x$vmfDensAuto, plot.numbers = TRUE,
                          main = paste("Variogram for log(biomass | presence) from \"",
                                       x$metrics, "\"",
                                       sep = ""),
                          ylim = c(0, 1.06 * max(c(x$vmDens$gamma, sum(x$vmfDensAuto$psill))))))
           }

           dev.off()

           ## Detailed predictions and prediction variance maps:
           png(file = file.path(ResultsPath,
                                "Geostatistics/Diagnostics",
                                paste("Untransformed_maps_", x$metrics,
                                      "(", gsub("[[:blank:]]", "", getOption("krigging.model")), ")",
                                      "_%01d.png", sep = "")),
               width = 800)
           ##
           if (is.null(x$models[["presence"]]))
           {
               plot(x$graphs[[2]])
               plot(x$graphs[[1]])
           }else{
               plot(x$graphs[[2]])
               plot(x$graphs[[1]])
               plot(x$graphs[[4]])
               plot(x$graphs[[3]])
           }
           dev.off()
       })

## Tiff outputs:
sapply(resGeostats,
       function(x)
       {
           writeRaster(x = raster(x$krigingDens["Dens"]),
                       filename = file.path(ResultsPath, "Geostatistics/Tiff",
                                            paste0("BiomassMap_from_", x$metrics,
                                                   "(",
                                                   gsub("[[:blank:]]", "",
                                                        getOption("krigging.model")),
                                                   ").tif")),
                       overwrite=TRUE)
       })


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
