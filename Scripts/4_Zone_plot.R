#-*- coding: latin-1 -*-

### File: 4_Zone_plot.R
### Time-stamp: <2018-11-19 12:26:35 yreecht>
###
### Created: 08/09/2016	12:28:15
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

breaksFake <- range(do.call(c,
                            sapply(densInterp[[1]],
                                   function(x) x@data[ , paste(names(densInterp)[1], "pred", sep=".")],
                                   simplify = FALSE)),
                    na.rm = TRUE)

## Conversion to polygones:
zones <- SPixDFlist2SPolyDF(densInterp[[1]],
                            variable=paste(names(densInterp)[1], "pred", sep="."),
                            breaks=breaksFake)


for (devType in getOption("surveyPlot.dev"))
{

    if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
    {
        width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6)
        dev <- openDev(device = devType,
                       directory = ResultsPath,
                       filename = c("Survey_zones"),
                       counter = FALSE,
                       width = width,
                       height = width *
                           diff(bboxMapZone[2, ]) / diff(bboxMapZone[1, ]),
                       pointsize = c(png = 36, jpg = 36, jpeg = 36, X11 = 11, pdf = 12),
                       verbose=FALSE)
    }else{
        warning("Device \"", devType, "\" not supported")
        next()
    }

    par(mar=c(0, 0, 0, 0)+0.1, xpd=FALSE,
        oma = c(0, 0,
                ifelse(getOption("surveyPlot.main"), 1, 0),
                0.0) + 0.2)

    ## Base map:
    if (is.element("SpatialGridDataFrame",
                   class(coastline)))
    {
        image(coastline, red = 1, green = 2, blue = 3,
              xlim=bboxMapZone[1, ],
              ylim=bboxMapZone[2, ],
              setParUsrBB = TRUE)
    }else{
        plot(coastline,
             col = getOption("surveyPlot.colGround"),
             bg = getOption("surveyPlot.colSea"),
             xlim=bboxMapZone[1, ],
             ylim=bboxMapZone[2, ])
    }

    plot(zones, col = "#C1CDF5", add = TRUE)

    tmpLegend <- c("Stations", ifelse(nrow(zones) - 1,
                                      "Survey zones",
                                      "Survey zone"))

    if (exists("surveyGrid") &&
        is.element("SpatialPolygonsDataFrame",
                   class(surveyGrid)))
    {
        plot(surveyGrid,
             col = NA,
             border = getOption("surveyPlot.colGrid"),
             ## bg = NA,
             add = TRUE)

        tmpLegend <- c(tmpLegend, "Survey grid")
    }

    ## if (isTRUE(getOption("surveyPlot.addStations")))
    ## {
    plot(stationData, add = TRUE,
         col = getOption("surveyPlot.colStations"),
         lwd = 2)
    ## }

    legend(x = getOption("surveyPlot.legendPos"),
           legend = tmpLegend,
           col = c(getOption("surveyPlot.colStations"), NA),
           lwd = c(2, NA, NA),
           fill = c(NA, "#C1CDF5", NA),
           border = c(NA, "black", getOption("surveyPlot.colGrid")),
           bg = "white", inset=0.02, cex = 0.8)

    if (getOption("surveyPlot.scaleBar"))
    {
        do.call(what = scale.bar, args = getOption("surveyPlot.scaleBar.opt"))
    }

    if (getOption("surveyPlot.northArrow"))
    {
        do.call(what = north.arrow, args = getOption("surveyPlot.northArrow.opt"))
    }

    if (devType != "X11") dev.off()
}





### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
