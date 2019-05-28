#-*- coding: latin-1 -*-

### File: 3_Size_distribution.R
### Time-stamp: <2018-10-04 16:13:17 yreecht>
###
### Created: 06/04/2016	11:42:13
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################


for (sp in row.names(speciesUnit))
{
    prev <- FALSE

    bdfCurrSp <- subset(biologicalData,
                        is.element(LD_SpeciesID, sp))

    if (nrow(bdfCurrSp) == 0)
    {
        message("\n## No size information for species \"", sp, "\"")
        next()
    }

    sizeRange <- range(bdfCurrSp$Size, na.rm = TRUE)
    densCurr <- density(na.omit(bdfCurrSp$Size),
                        bw = getOption("surveyPlot.density.bw"))
    densRange <- range(c(0, densCurr$y))

    ## Legend (Month Year):
    leg <- format(##
                  if(all(is.na(bdfCurrSp$LD_EventStartDate)))
                  {
                      head(stationData@data[(! is.na(stationData@data$TIMESTART)) &
                                            (as.character(stationData@data$TRACKID) ==
                                             as.character(head(bdfCurrSp$LD_HaulNo, 1))) ,
                                            "TIMESTART"],
                           1)
                  }else{
                      head(na.omit(bdfCurrSp$LD_EventStartDate),
                           1)
                  },
                  format = "%b %Y")

    colLeg <- c("dodgerblue4")

    if (exists("biologicalDataPrev") &&
        ! is.null(biologicalDataPrev) &&
        is.element(sp, biologicalDataPrev$LD_SpeciesID)) # Checks presence of species in previous survey.
    {
        bdfPrevSp <- subset(biologicalDataPrev,
                            is.element(LD_SpeciesID, sp))

        prev <- as.logical(nrow(bdfPrevSp))

        sizeRange <- range(c(sizeRange,
                             range(bdfPrevSp$Size, na.rm = TRUE)))

        densPrev <- density(na.omit(bdfPrevSp$Size),
                            bw = getOption("surveyPlot.density.bw"))
        densRange <- range(c(densRange, densPrev$y))

        freqRange <- range(densPrev$y) * sum(! is.na(bdfPrevSp$Size))

        if (exists("totalEffortPrev") && ! is.null(totalEffortPrev))
        {
            abDensRange <- freqRange / totalEffortPrev
        }else{
            abDensRange <- c(0, 0)
        }

        leg <- c(leg,
                 ifelse(ncol(bdfPrevSp) > 2,
                        format(bdfPrevSp$LD_EventStartDate[1], format = "%b %Y"),
                        "Previous"))

        colLeg <- c(colLeg, "deepskyblue")
    }else{
        freqRange <- NULL
        abDensRange <- NULL
    }

    histCurr <- hist(na.omit(bdfCurrSp$Size),
                     breaks = seq(from = floor(sizeRange[1]) - 0.5,
                                  to = ceiling(sizeRange[2]) + 0.5),
                     plot = FALSE)

    ## ylims for density and frequency
    densRange <- range(c(densRange, histCurr$density))

    freqRange <- floor(range(c(range(c(0, densCurr$y)) * sum(histCurr$counts),
                               freqRange,
                               histCurr$counts))) + c(0, 1)

    ## ##################################################
    ## Ground density histogram:

    ## Total effort in squared m:
    totalEffort <- sum(stationData@data$effort[as.character(stationData@data$TRACKID) %in%
                                               as.character(catchData$LD_HaulNo)],
                       na.rm = TRUE)

    histCurrMod <- histCurr
    histCurrMod$density <- histCurrMod$counts / totalEffort

    abDensRange <-range(c(abDensRange, histCurrMod$density))

    ##
    MLSadd <- TRUE
    ## Graphical representation:
    for (devType in getOption("surveyPlot.dev"))
    {
        if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
        {
            width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6.5)
            dev <- openDev(device = devType,
                           directory = ResultsPath,
                           filename = paste0("Size_distribution_",
                                             gsub("[[:blank:]]+", "_", sp)),
                           width = width,
                           height = width * (5 / 9),
                           pointsize = c(png = 24, jpg = 24, jpeg = 24, X11 = 12, pdf = 12),
                           counter = FALSE,
                           verbose=FALSE)
        }else{
            warning("Device \"", devType, "\" not supported")
            next()
        }

        par(mgp = c(1.9, 0.7, 0), mar = c(3, 3, 2.5, 0.5) + 0.1)
        plot(histCurr, freq = FALSE, border = "grey", ylim = densRange,
             xlab = "Size",
             main = if(getOption("surveyPlot.main")){bquote(italic(.(sp)))}else{""})


        if (! is.null(MLS) &&
            is.element(sp, names(MLS)))
        {
            colLeg <- c(colLeg, "firebrick3")

            if (MLSadd) leg <- c(leg, "MLS") ## paste("MLS (=", MLS[sp], " mm)", sep = ""))
            MLSadd <- FALSE

            abline(v = MLS[sp], col = tail(colLeg, 1), lwd = 3)
        }

        if (prev) lines(densPrev, lwd = 3, col = colLeg[2])
        lines(densCurr, lwd = 3, col = colLeg[1])

        legend(x = "topright", legend = leg, col = colLeg, lwd = 3,
               bg = "white")

        if (devType != "X11") dev.off()
    }

    for (devType in getOption("surveyPlot.dev"))
    {
        if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
        {
            width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6.5)
            dev <- openDev(device = devType,
                           directory = ResultsPath,
                           filename = paste0("Size_distribution_freq_",
                                             gsub("[[:blank:]]+", "_", sp)),
                           width = width,
                           height = width * (5 / 9),
                           pointsize = c(png = 24, jpg = 24, jpeg = 24, X11 = 12, pdf = 12),
                           counter = FALSE,
                           verbose=FALSE)
        }else{
            warning("Device \"", devType, "\" not supported")
            next()
        }

        par(mgp = c(1.9, 0.7, 0), mar = c(3, 3, 2.5, 0.5) + 0.1)
        plot(histCurr, freq = TRUE, border = "grey", ylim = freqRange,
             xlab = "Size",
             main = if(getOption("surveyPlot.main")){bquote(italic(.(sp)))}else{""})


        if (! is.null(MLS) &&
            is.element(sp, names(MLS)))
        {
            colLeg <- c(colLeg, "firebrick3")

            if (MLSadd) leg <- c(leg, "MLS") ## paste("MLS (=", MLS[sp], " mm)", sep = ""))
            MLSadd <- FALSE

            abline(v = MLS[sp], col = tail(colLeg, 1), lwd = 3)
        }

        if (prev)
        {
            freqPrev <- densPrev
            freqPrev$y <- freqPrev$y * sum(! is.na(bdfPrevSp$Size))
            lines(freqPrev, lwd = 3, col = colLeg[2])
        }
        freqCurr <- densCurr
        freqCurr$y <- freqCurr$y * sum(histCurr$counts)
        lines(freqCurr, lwd = 3, col = colLeg[1])

        legend(x = "topright", legend = leg, col = colLeg, lwd = 3,
               bg = "white")

        if (devType != "X11") dev.off()
    }

    for (devType in getOption("surveyPlot.dev"))
    {
        if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
        {
            width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6.5)
            dev <- openDev(device = devType,
                           directory = ResultsPath,
                           filename = paste0("Size_distribution_abDensity_",
                                             gsub("[[:blank:]]+", "_", sp)),
                           width = width,
                           height = width * (5 / 9),
                           pointsize = c(png = 24, jpg = 24, jpeg = 24, X11 = 12, pdf = 12),
                           counter = FALSE,
                           verbose=FALSE)
        }else{
            warning("Device \"", devType, "\" not supported")
            next()
        }

        par(mgp = c(1.9, 0.7, 0), mar = c(3, 3.2, 2.5, 0.5) + 0.1)
        plot(histCurrMod, freq = FALSE, border = "grey", ylim = abDensRange,
             xlab = "Size",
             ylab = expression("Density"~(indiv %.% m^{-2})),
             main = if(getOption("surveyPlot.main")){bquote(italic(.(sp)))}else{""})


        if (! is.null(MLS) &&
            is.element(sp, names(MLS)))
        {
            colLeg <- c(colLeg, "firebrick3")

            if (MLSadd) leg <- c(leg, "MLS") ## paste("MLS (=", MLS[sp], " mm)", sep = ""))
            MLSadd <- FALSE

            abline(v = MLS[sp], col = tail(colLeg, 1), lwd = 3)
        }

        ## Legend for the effort:
        legEffort <- parse(text = paste0("bold(\"", leg[1], ":\")~",
                                         prettyNum(totalEffort, format = "g", digits = 3),
                                         "~m^2"))


        if (prev && exists("totalEffortPrev") && ! is.null(totalEffortPrev))
        {
            abDensPrev <- densPrev
            abDensPrev$y <- abDensPrev$y * sum(! is.na(bdfPrevSp$Size)) / totalEffortPrev
            lines(abDensPrev, lwd = 3, col = colLeg[2])

            legEffort <- c(legEffort,
                           parse(text = paste0("\"", leg[2], ":\"~",
                                               prettyNum(totalEffortPrev, format = "g", digits = 3),
                                               "~m^2")))
        }

        abDensCurr <- freqCurr
        abDensCurr$y <- abDensCurr$y / totalEffort
        lines(abDensCurr, lwd = 3, col = colLeg[1])

        legend(x = "topright", legend = leg, col = colLeg, lwd = 3,
               bg = "white")

        legend(x = "right", legend = legEffort,
               title = expression(bold("Total sampling effort:")),
               bty = "n", inset = 0.01)

        if (devType != "X11") dev.off()
    }
}





### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
