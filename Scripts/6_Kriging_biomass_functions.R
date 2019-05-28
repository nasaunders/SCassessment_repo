#-*- coding: utf-8 -*-

### File: 6_Kriging_biomass_functions.R
### Time-stamp: <2019-05-20 16:02:32 yreecht>
###
### Created: 05/04/2019	08:43:30
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


surveyKrigingVariograms <- function(formulaPres = NULL,
                                    formulaDens = NULL,
                                    varName = NULL,
                                    stationPosData,
                                    autoBreak = getOption("surveyGeostat.vg.autobreaks"),
                                    cutoff = getOption("surveyGeostat.vg.cutoff"),
                                    width = getOption("surveyGeostat.vg.width"),
                                    use.lik = getOption("surveyGeostat.vg.likelihood"),
                                    lik.method = getOption("surveyGeostat.vg.lik.method"),
                                    krigging.model = getOption("krigging.model"),
                                    lambda = NULL, offset.bc = NULL, offset.name = "cst",
                                    verbose = TRUE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  5 Apr 2019, 08:46

    assign(x = offset.name, value = offset.bc)

    if (! is.null(formulaPres))
    {
        if (isTRUE(autoBreak))
        {
            ## Calculating breaks for the variograms:
            stationDistances <- as.vector(dist(coordinates(stationPosData)))
            breaksPres <- c(0, tail(quantile(stationDistances[stationDistances <
                                                              cutoff],
                                             probs = seq(0, 1, by = max(0.05, 30 / length(stationDistances)))
                                             ##c(0, 0.05, 0.10,
                                                     ##   seq(0.15, 0.2, 0.05), seq(0.25, 0.85, 0.1), 1)
                                             ),
                                    -1))

            ## Empirical indicator variogram:
            vi <- variogram(formulaPres,
                            location = subset(stationPosData, ! is.na(pres)),
                            ## width = width,
                            cutoff = cutoff,
                            boundaries = breaksPres)
        }else{
            ## Empirical indicator variogram:
            vi <- variogram(formulaPres,
                            location = subset(stationPosData, ! is.na(pres)),
                            width = width,
                            cutoff = cutoff)
        }

        ## Observation likelihood or empirical least-squares fit?:
        if (isTRUE(use.lik))
        {
            ## Max likelihood fitting of variograms and selection of the best model (BIC):
            vifAuto <- fitVariogramLike(SpObj = stationPosData,
                                        formula = formulaPres,
                                        model = c("Sph", "Exp", "Mat"),
                                        lik.method = lik.method)
        }else{
            ## Choose automatically the best LS fit of the empirical variogram:
            vifAuto <- fit.variogram(vi,
                                     vgm(c("Exp", "Sph",
                                           "Ste", "Mat")))
        }

        if (verbose)
        {
            print(vifAuto)
            message("")
        }
    }else{
        vi  <- NULL
        vifAuto  <-  NULL
    }

    ## Density:
    if(! is.null(formulaDens))
    {
        ## Variogram:
        if (isTRUE(autoBreak))
        {
            ## Calculating breaks for the variograms:
            idxP <- ## is.null(formulaPres) |
                if(! is.null(formulaPres))
                {
                    (as.logical(stationPosData$pres) &
                     ! is.na(as.logical(stationPosData$pres)))
                }else{
                    !is.na(stationPosData@data[ , varName])
                }

            if (is.null(idxP)) idxP <- rep(TRUE, nrow(stationPosData))

            stationDistancesPres <-
                as.vector(dist(coordinates(stationPosData[idxP, ])))
            breaks <- c(0, tail(quantile(stationDistancesPres[stationDistancesPres <
                                                              cutoff],
                                         probs = seq(0, 1, by = max(0.05, 30 / length(stationDistancesPres)))
                                         ## c(0, 0.025, 0.05,
                                         ##           seq(0.10, 0.2, 0.05), seq(0.25, 0.85, 0.1), 1)
                                         ),
                                -1))

            vmDens <- variogram(formulaDens,
                                subset(stationPosData,
                                       ! is.na(stationPosData@data[ , varName]) &
                                       stationPosData@data[ , varName] > 0),
                                ## width = width,
                                cutoff = cutoff,
                                boundaries = breaks)
        }else{
            vmDens <- variogram(formulaDens,
                                subset(stationPosData,
                                       ! is.na(stationPosData@data[ , varName]) &
                                       stationPosData@data[ , varName] > 0),
                                width = width,
                                cutoff = cutoff)
        }

        ## Observation likelihood or empirical least-squares fit?:
        if (isTRUE(use.lik))
        {
            ## Max likelihood fitting of variograms and selection of the best model (BIC):
            subDataVarioFit <- subset(stationPosData,
                                      ! is.na(stationPosData@data[ , varName]) &
                                      stationPosData@data[ , varName] > 0)

            ## Drop levels for possibly absent covariates (e.g. zone with no presence):
            subDataVarioFit@data <- droplevels(subDataVarioFit@data)

            vmfDensAuto <- fitVariogramLike(SpObj = subDataVarioFit,
                                            formula = formula(paste0("densTransfo ",
                                                                     krigging.model)),
                                            model = c("Sph", "Exp", "Mat"),
                                            lik.method = lik.method)
        }else{
            ## Choose automatically the best LS fit of the empirical variogram:
            vmfDensAuto <- fit.variogram(vmDens,
                                         vgm(c("Exp", "Sph", "Ste", "Mat")))
        }

        if (verbose)
        {
            print(vmfDensAuto)
            message("")
        }
    }else{
        vmDens <- NULL
        vmfDensAuto <- NULL
    }

    return(list(vifPresAuto = vifAuto,
                vmfDensAuto = vmfDensAuto,
                viPres = vi,
                vmDens = vmDens))
}

boxCoxDiag <- function(stationPosData,
                       varName = densBiomassName,
                       presenceOnly,
                       ResultsPath,
                       offset.bc = 0)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  5 Apr 2019, 11:07

    if (isTRUE(presenceOnly))
    {
        dataPres <- stationPosData@data[(! is.na(stationPosData@data[ , varName]) &
                                         stationPosData@data[ , varName] > 0),
                                        varName]

        idxPres <- (! is.na(stationPosData@data[ , varName]) &
                    stationPosData@data[ , varName] > 0)
    }else{
        dataPres <- stationPosData@data

        idxPres <- (! is.na(stationPosData@data[ , varName]))
    }

    ## Dist:
    ## X11(height = 3.5)
    png(file = file.path(ResultsPath,
                         "Geostatistics/Diagnostics",
                         paste0("Normality_", metricsName, ".png")),
        width = 500, height = 500)
    par(mfcol = c(2, 2))
    hist(dataPres,
         main = metricsName, xlab = "Dens. biomass", breaks = 20)

    ## Fitting the lambda for the
    lambdas <- boxcox(dataPres + offset.bc ~ 1,
                      plotit = TRUE, lambda = seq(-5, 5, .005))
    mtext(text = expression(BoxCox~lambda~max~likelihood~fit), cex = 1.2)
    lambda <- lambdas$x[which.max(lambdas$y)]

    if (abs(lambda) > 2) warnings("Lambda = ", round(lambda, 3),
                                  ": BoxCox transformation might not be relevant!")

    ## Normality diagnostics after transformation:
    hist(boxcox.transf(x = dataPres, lambda = lambda, offset = offset.bc),
         main = "BoxCox transformed (| presence)",
         xlab = "BoxCox(Dens. biomass)", breaks = 20, probability = TRUE)
    legend(x = "topright", legend = parse(text = paste0("lambda==", round(lambda, 2))),
           bty = "n", cex = 1.5)
    lines(density(boxcox.transf(x = dataPres, lambda = lambda, offset = offset.bc),
                  na.rm = TRUE, adjust = 2),
          col = "blue", lwd = 2)
    qqnorm(boxcox.transf(x = dataPres, lambda = lambda, offset = offset.bc))
    qqline(boxcox.transf(x = dataPres, lambda = lambda, offset = offset.bc),
           col = "red", lwd = 2)
    dev.off()

    return(lambda)
}

surveyKriging <- function(formulaDens, vmfDens,
                          formulaPres = NULL, vifPres = NULL,
                          stationPosData,
                          data.pix,
                          varName,
                          correctPres = TRUE,
                          lambda, offset.bc, offset.name = "cst",
                          block = NULL, ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  5 Apr 2019, 11:36

    assign(x = offset.name, value = offset.bc)

    if (is.null(block))
    {
        block <- data.pix@grid@cellsize
    }else{}

    ## Optional presence model first:

    if (!is.null(formulaPres))
    {
        ## Indicator krigging:
        kPres <- krige(formula = formulaPres,
                       locations = subset(stationPosData, ! is.na(pres)),
                       newdata = data.pix, model = vifPres,
                       indicators = TRUE,
                       block = block,
                       ...)

        ## Correction of presence probabilities out of the [0,1] range:
        if (isTRUE(correctPres))
        {
            kPres[["var1.pred"]][!is.na(kPres[["var1.pred"]]) &
                                 kPres[["var1.pred"]] < 0] <- 0

            kPres[["var1.pred"]][!is.na(kPres[["var1.pred"]]) &
                                 kPres[["var1.pred"]] > 1] <- 1
        }

        ## Index of kept data for presence only kriging:
        idx <- (! is.na(stationPosData@data[ , varName]) &
                stationPosData@data[ , varName] > 0)


    }else{
        kPres <- NULL

        ## Index of kept data for kriging:
        idx <- ! is.na(stationPosData@data[ , varName])
    }

    ## Subset of data for kriging variable of interest:
    subDataPred <- subset(stationPosData,
                          idx)

    ## Drop levels for possibly absent covariates (e.g. zone with no presence):
    subDataPred@data <- droplevels(subDataPred@data)

    suppressWarnings(kDensPos <- krige(formula = formulaDens,
                                       locations = subDataPred,
                                       newdata = data.pix, model = vmfDens,
                                       nsim = 0, indicators = FALSE,
                                       block = block,
                                       ...))

    ## Cox-box back transformation with bias correction:
    kDensPos[["DensPos"]] <- boxcox.mean.backtransf(mean = kDensPos[["var1.pred"]],
                                                    var = kDensPos[["var1.var"]],
                                                    lambda = lambda, offset = offset.bc)

    if (!is.null(formulaPres))
    {
        ## Sometimes the backtransformation using mean and var returns NaN
        ##   for very low values => replaced by zeros):
        kDensPos[["DensPos"]][is.nan(kDensPos[["DensPos"]]) &
                              ! is.na(kPres[["var1.pred"]])] <- 0

        ## summary(kDensPos[["DensPos"]])

        ## Combination of presence and biomass density predictions (with correction for negative values):
        kDensPos[["Dens"]] <- ifelse(kDensPos[["DensPos"]] * kPres[["var1.pred"]] >= 0,
                                     kDensPos[["DensPos"]] * kPres[["var1.pred"]],
                                     0) #Correction for possible negative predictions.

    }else{
        ## Correction for possible negative predicted values:
        kDensPos[["Dens"]] <- ifelse(kDensPos[["DensPos"]] >= 0,
                                     kDensPos[["DensPos"]],
                                     0)
    }

    return(list(krigingPres = kPres,
                krigingDens = kDensPos))
}


surveyGaussianSimu <- function(formulaDens, vmfDens,
                               formulaPres = NULL, vifPres = NULL,
                               stationPosData,
                               data.pix,
                               varName,
                               correctPres = TRUE,
                               lambda, offset.bc, offset.name = "cst",
                               nmax = getOption("surveyGeostat.CGS.nmax"),
                               nsim = getOption("surveyGeostat.CGS.nsim"),
                               debug.level = -1, ...)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  5 Apr 2019, 15:13

    assign(x = offset.name, value = offset.bc)

    ## Optional presence model first:
    if (!is.null(formulaPres))
    {
        dataPres <- subset(stationPosData, ! is.na(pres))
        ## gc()

        ## Presence/absence simulation:
        kPresSim <- tryCatch(krige(formula = formulaPres,
                                   locations = dataPres,
                                   newdata = data.pix,
                                   nsim = nsim,
                                   indicators = TRUE,
                                   model = vifPres,
                                   nmax = nmax,
                                   debug.level = 0), # For some reason, first run of indicator
                                        # kriging fails when previous run was not indicator
                             error = function(e)
                             {
                                 krige(formula = formulaPres,
                                       locations = dataPres,
                                       newdata = data.pix,
                                       nsim = nsim,
                                       indicators = TRUE,
                                       model = vifPres,
                                       nmax = nmax,
                                       debug.level = debug.level,...)
                             })

        summary(subset(stationPosData, ! is.na(pres)))
        head(data.pix)
        ## Index of kept data for presence only kriging:
        idx <- (! is.na(stationPosData@data[ , varName]) &
                stationPosData@data[ , varName] > 0)
    }else{
        kPresSim <- NULL

        ## Index of kept data for kriging:
        idx <- ! is.na(stationPosData@data[ , varName])
    }

    ## Subset of data for kriging variable of interest:
    subDataPred <- subset(stationPosData,
                          idx)

    ## Drop levels for possibly absent covariates (e.g. zone with no presence):
    subDataPred@data <- droplevels(subDataPred@data)

    ## Conditional gaussian simulations of density:
    kDensPosSim <- krige(formula = formulaDens,
                         locations = subDataPred,
                         newdata = data.pix,
                         model = vmfDens,
                         nmax = nmax,
                         nsim = nsim,
                         block = numeric(0),
                         indicators = FALSE,
                         debug.level = debug.level,...)

    kDensSim <- kDensPosSim

    if (!is.null(formulaPres))
    {
        kDensSim@data <- boxcox.transf(x = kDensSim@data,
                                       lambda = lambda, offset = offset.bc,
                                       backtransformation = TRUE) * kPresSim@data

        kDensSim@data[(is.na(kDensSim@data) &
                       ! is.na(kPresSim@data)) |
                      kDensSim@data < 0] <- 0
    }else{
        kDensSim@data <- boxcox.transf(x = kDensSim@data[ , ],
                                       lambda = lambda, offset = offset.bc,
                                       backtransformation = TRUE)

        kDensSim@data[kDensSim@data < 0] <- 0
    }

    ## Total simulated abundance by km^2
    abdSimu <- tryCatch(sapply(kDensSim@data, sum, na.rm = TRUE) *
                        prod(data.pix@grid@cellsize) / 1e3,
                        error=function(e) return(NULL))

    ## hist(abdSimu, breaks = 25)

    return(list(kDensTransfSim = kDensPosSim,
                kDensSim = kDensSim,
                abdSimu = abdSimu,
                kPresSim = kPresSim))
}

plotDensityMap <- function(raster, variable, palpers, legendTitle)
{
    library(RColorBrewer)

    oldpar <- par("mar", no.readonly=TRUE)
    par(mar=c(0, 0, 0, 0)+0.1, xpd=FALSE,
        oma = c(0, 0,
                ifelse(getOption("surveyPlot.main"), 1, 0),
                0) + 0.2)

    ## Base map:
    if (is.element("SpatialGridDataFrame",
                   class(coastline)))
    {
        image(coastline, red = 1, green = 2, blue = 3,
              xlim=bboxMapZone[1, ],
              ylim=bboxMapZone[2, ]## ,
              ## setParUsrBB = TRUE
              )
    }else{
        plot(coastline,
             col = getOption("surveyPlot.colGround"),
             bg = getOption("surveyPlot.colSea"),
             xlim=bboxMapZone[1, ],
             ylim=bboxMapZone[2, ])
    }



    ## breaksCommon <- prettyQuantileBreaks(raster[[variable]], n.categ = 100, n.digits=4)
    if (diff(range(raster[[variable]], na.rm = TRUE)) == 0)
    {
        ## value +/- 1%:
        rangeRast <- mean(raster[[variable]], na.rm = TRUE) * c(0.99, 1.01)
        breaksRange <- prettyQuantileBreaks(rangeRast, n.categ = 100, n.digits=4)

        pl <- plotRasterColor(raster = raster, variable = variable,
                              add = TRUE, persPalette = palpers, Quantile.type=2,
                              breaks =  breaksRange)
    }else{
        ## Most cases:
        ## breaksRange <- prettyQuantileBreaks(raster[[variable]], n.categ = 100, n.digits=4)

        ## breaksRange[breaksRange < 0] <- 0

        pl <- plotRasterColor(raster = raster, variable = variable,
                              add = TRUE, persPalette = palpers, Quantile.type=2,
                              ## breaks =  breaksRange)
                              n.categ = 100)
    }

    plbreaks <- pl$breaks
    plbreaks[plbreaks < 0] <- 0

    do.call(what = colorBarSmall,
            args = c(list(breaks = plbreaks, col = pl$col,
                          title = legendTitle),
                     getOption("surveyPlot.colorBarSmall.opt")))

    scale.bar(x = "bottomleft", convert = "Nm", inset = c(0.02, 0.07))

    return(invisible(pl))
}


surveyKriging.CrossValidation <- function(formulaDens, vmfDens,
                                          formulaPres = NULL, vifPres = NULL,
                                          stationPosData,
                                          varName,
                                          correctPres = TRUE,
                                          simuVal = TRUE,
                                          krigVal = TRUE,
                                          lambda, offset.bc, offset.name = "cst")
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  8 Apr 2019, 18:32

    res <- sapply(1:nrow(stationPosData),
                  function(i)
                  {
                      resPred <- tryCatch(surveyKriging(formulaDens = formulaDens,
                                                        vmfDens = vmfDens,
                                                        formulaPres = formulaPres,
                                                        vifPres = vifPres,
                                                        stationPosData = stationPosData[-i, ],
                                                        data.pix = stationPosData[i, ],
                                                        varName = varName,
                                                        correctPres = correctPres,
                                                        lambda = lambda,
                                                        offset.bc = offset.bc,
                                                        offset.name = offset.name,
                                                        block=numeric(0),
                                                        debug.level = 0),
                                          error = function(e)return(NULL))

                      if (is.null(resPred)) return(NULL)

                      if (! is.null(formulaPres))
                      {
                          ZscPres <- (stationPosData@data[i, "pres"] -
                                      resPred$krigingPres@data[ , "var1.pred"]) /
                              sqrt(resPred$krigingPres@data[ , "var1.var"])
                      }else{
                          ZscPres <- NA
                      }

                      ZscDens <- (stationPosData@data[i, "densTransfo"] -
                                  resPred$krigingDens@data[ , "var1.pred"]) /
                          sqrt(resPred$krigingDens@data[ , "var1.var"])

                      biasDens <- (stationPosData@data[i, varName] -
                                   resPred$krigingDens@data[ , "Dens"])

                      if (krigVal)
                      {
                          resKrig <- krigValidation(stations = stationPosData[i, ],
                                                    densKrig = resPred$krigingDens,
                                                    presKrig = resPred$krigingPres,
                                                    varName = varName,
                                                    boxcox.params= c(lambda, offset.bc),
                                                    prop=0.95, returnData=TRUE)

                          resKrig <- resKrig$data$stations@data
                      }else{
                          resKrig <- NULL
                      }

                      return(cbind(data.frame(ZscDens = ZscDens, biasDens = biasDens, ZscPres = ZscPres),
                                   resKrig))

                  }, simplify = FALSE)

    res <- do.call(rbind, res)


    if (isTRUE(krigVal))
    {
        ME.Krig <- mean(res$test1, na.rm = TRUE)
        RMSE.Krig <- sqrt(mean(res$test1^2, na.rm = TRUE))

        Perc.Krig <- mean(res$test2, na.rm = TRUE)
    }else{
        ME.Krig <- NULL
        RMSE.Krig <- NULL

        Perc.Krig <- NULL
    }

    if (isTRUE(simuVal))
    {
        simuRes <-  sapply(1:nrow(stationPosData),
                           function(i)
                           {
                               resSim <- tryCatch(surveyGaussianSimu(formulaDens = formulaDens,
                                                                     vmfDens = vmfDens,
                                                                     formulaPres = formulaPres,
                                                                     vifPres = vifPres,
                                                                     stationPosData = stationPosData[-i, ],
                                                                     data.pix = stationPosData[i, ],
                                                                     varName = varName,
                                                                     correctPres = correctPres,
                                                                     lambda = lambda,
                                                                     offset.bc = offset.bc,
                                                                     offset.name = offset.name,
                                                                     debug.level = 0),
                                                  error = function(e)return(NULL))

                               if (is.null(resSim)) return(NULL)

                               test <- simValidation(stations = stationPosData[i, ],
                                                     densSimu = resSim$kDensSim,
                                                     varName=varName,
                                                     prop=0.95,
                                                     returnData=TRUE)

                               return(test$data$stations@data)
                           }, simplify = FALSE)

        simuRes <- do.call(rbind, simuRes)

        ME.Sim <- mean(simuRes$test1, na.rm = TRUE)
        RMSE.Sim <- sqrt(mean(simuRes$test1^2, na.rm = TRUE))

        Perc.Sim <- mean(simuRes$test2, na.rm = TRUE)
    }else{
        simuRes <- NULL
        ME.Sim <- NULL
        RMSE.Sim <- NULL

        Perc.Sim <- NULL
    }

    ## summary(res)
    ## sapply(res, mean, na.rm = TRUE)
    ## sapply(res, sd, na.rm = TRUE)

    ## plot(seq(-2.5, 2.5, 0.01), dnorm(seq(-2.5, 2.5, 0.01)),
    ##      lwd = 2, col = "red", type = "l", lty = 2)
    ## lines(density(res$ZscDens, na.rm = TRUE, bw = 0.5), col = "darkblue")
    ## lines(density(res$ZscPres, na.rm = TRUE, bw = 0.8), col = "darkgreen")

    return(list(resDF = res,
                resDFsimu = simuRes,
                meanZsc = sapply(res[ , 1:3], mean, na.rm = TRUE),
                sdZsc = sapply(res[ , 1:3], sd, na.rm = TRUE),
                simuValid = c(ME.Sim = ME.Sim,
                              RMSE.Sim = RMSE.Sim,
                              Perc.Sim = Perc.Sim),
                krigingValid = c(ME.Krig = ME.Krig,
                                 RMSE.Krig = RMSE.Krig,
                                 Perc.Krig = Perc.Krig)))

}

krigValidation <- function(stations, densKrig, presKrig = NULL,
                           varName, presVarName = "pres",
                           boxcox.params,
                           prop = 0.95, returnData = FALSE)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  9 Apr 2019, 12:48

    densKrigStat <- densKrig["Dens"]

    probLow <- (1 - prop^(1 / ifelse(is.null(presKrig), 1, 2))) / 2
    probUpp <- 1 - (1 - prop^(1 / ifelse(is.null(presKrig), 1, 2))) / 2

    ## densKrigStat$lwr <- qnorm(p = probLow, RainSqrt$var1.pred, sqrt(RainSqrt$var1.var))

    ## transfDens <- boxcox.transf(stations@data[ , c(varName)],
    ##                             lambda = boxcox.params[1],
    ##                             offset = boxcox.params[2]+
    ##                                 min(stations@data[stations@data[ , c(varName)] > 0 , c(varName)],
    ##                                     na.rm = TRUE) / 10)

    ## covDP <- cov(transfDens,
    ##              stations@data[ , c(presVarName)])

    ## covSqDP <- cov(transfDens^2,
    ##                stations@data[ , c(presVarName)]^2)

    ## ## covDP <- cov(densKrig$var1.pred,
    ## ##              presKrig$var1.pred)

    ## ## covSqDP <- cov(densKrig$var1.pred^2,
    ## ##                presKrig$var1.pred^2)

    ## densKrigStat$varDens <- (covSqDP + (densKrig$var1.var + densKrig$var1.pred ^ 2) *
    ##                          (presKrig$var1.var + presKrig$var1.pred ^ 2)) -
    ##     (covDP + densKrig$var1.pred * presKrig$var1.pred) ^ 2

    ## densKrigStat$predDens <- densKrig$var1.pred *
    ##     sapply(presKrig$var1.pred,
    ##            function(x) min(max(x, 0), 1))

    ## summary(densKrigStat$varDens)

    ## densKrigStat$densBT <- boxcox.mean.backtransf(mean = densKrigStat$predDens,
    ##                                               var = densKrigStat$varDens,
    ##                                               lambda = boxcox.params[1],
    ##                                               offset = boxcox.params[2])

    densKrigStat$dens.lwr <-  qnorm(p = probLow,
                                    mean = densKrig$var1.pred,
                                    sd = sqrt(densKrig$var1.var))

    densKrigStat$dens.upr <-  qnorm(p = probUpp,
                                    mean = densKrig$var1.pred,
                                    sd = sqrt(densKrig$var1.var))

    densKrigStat$densBT.lwr <- boxcox.transf(x = densKrigStat$dens.lwr,
                                             lambda = boxcox.params[1],
                                             offset = boxcox.params[2],
                                             backtransformation = TRUE)


    densKrigStat$densBT.upr <- boxcox.transf(x = densKrigStat$dens.upr,
                                             lambda = boxcox.params[1],
                                             offset = boxcox.params[2],
                                             backtransformation = TRUE)

    if (!is.null(presKrig))
    {
        densKrigStat$pres.lwr <-  qnorm(p = probLow,
                                        mean = presKrig$var1.pred,
                                        sd = sqrt(presKrig$var1.var))

        densKrigStat$pres.upr <-  qnorm(p = probUpp,
                                        mean = presKrig$var1.pred,
                                        sd = sqrt(presKrig$var1.var))

        densKrigStat$Density.lwr <- densKrigStat$densBT.lwr * densKrigStat$pres.lwr
        densKrigStat$Density.upr <- densKrigStat$densBT.upr * densKrigStat$pres.upr
    }else{
        densKrigStat$Density.lwr <- densKrigStat$densBT.lwr
        densKrigStat$Density.upr <- densKrigStat$densBT.upr
    }
    stations@data[ , c("lwr.Krig",
                       "prd.Krig",
                       "upr.Krig")] <- over(stations,
                                            densKrigStat)[ , c("Density.lwr", "Dens", "Density.upr")]

    ## Estimation of Mean Error (ME) and Root Mean Square Error (RMSE):
    stations$test1 <- stations[[varName]] - stations$prd.Krig

    me.Krig <- mean(stations$test1, na.rm = TRUE)
    rmse.Krig <- sqrt(mean(stations$test1^2, na.rm = TRUE))

    ## Percentile validation...
    ##   What proportion of observations fall in the "prop" prediction interval:
    stations$test2 <- apply(stations@data[ , c(varName, "lwr.Krig", "upr.Krig")],
                            1,
                            function(x) x[1] > x[2] & x[1] < x[3])

    perc.Krig <- mean(stations$test2, na.rm = TRUE)

    return(list(ME.Krig = me.Krig, RMSE.Krig = rmse.Krig,
                prop.Intvl.Krig = perc.Krig,
                params = c(varName = varName,
                           prop = prop),
                data = if (isTRUE(returnData))
                       {
                           list(stations = stations[c(varName, "prd.Krig",
                                                      "lwr.Krig", "upr.Krig", "test1", "test2")],
                                densKrigStat = densKrigStat)
                       }else{NULL}))
}

simValidation <- function(stations, densSimu, varName,
                          prop = 0.95, returnData = FALSE,
                          centralFun = mean, simuStat = NULL)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  9 Apr 2019, 11:41

    probLow <- (1 - prop) / 2
    probUpp <- 1 - (1 - prop) / 2

    if ((missing(densSimu) || is.null(densSimu)) && is.null(simuStat))
        stop("One of the parameters \"densSimu\" and \"simuStat\" must be provided")

    if (missing(simuStat) || is.null(simuStat))
    {
        simuStat <- densSimu[0]

        simuStat$central <- apply(densSimu@data, 1, centralFun, na.rm = TRUE)
        simuStat$sd <- apply(densSimu@data, 1, sd, na.rm = TRUE)
        simuStat$lwr <- apply(densSimu@data, 1, quantile, probs = probLow, na.rm = TRUE)
        simuStat$upr <- apply(densSimu@data, 1, quantile, probs = probUpp, na.rm = TRUE)
    }

    stations@data[ , c("lwr.Sim",
                       "prd.Sim",
                       "upr.Sim")] <- over(stations,
                                           simuStat)[ , c("lwr", "central", "upr")]

    ## Estimation of Mean Error (ME) and Root Mean Square Error (RMSE):
    stations$test1 <- stations[[varName]] - stations$prd.Sim

    me.Sim <- mean(stations$test1, na.rm = TRUE)
    rmse.Sim <- sqrt(mean(stations$test1^2, na.rm = TRUE))

    ## Percentile validation...
    ##   What proportion of observations fall in the "prop" prediction interval:
    stations$test2 <- apply(stations@data[ , c(varName, "lwr.Sim", "upr.Sim")],
                            1,
                            function(x) x[1] > x[2] & x[1] < x[3])

    perc.Sim <- mean(stations$test2, na.rm = TRUE)


    return(list(ME.Sim = me.Sim, RMSE.Sim = rmse.Sim,
                prop.Intvl.Sim = perc.Sim,
                params = c(varName = varName,
                           prop = prop,
                           centralFun = centralFun),
                data = if (isTRUE(returnData))
                       {
                           list(stations = stations[c(varName, "prd.Sim",
                                                      "lwr.Sim", "upr.Sim", "test1", "test2")],
                                simuStat = simuStat)
                       }else{NULL}))
}




### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
