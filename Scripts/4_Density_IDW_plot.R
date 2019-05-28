#-*- coding: latin-1 -*-

### File: 4_Density_IDW_plot.R
### Time-stamp: <2018-10-03 15:32:22 yreecht>
###
### Created: 22/03/2016	14:23:01
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

sapply(names(densInterp),
       function(i)
   {
       ## Graphical representation:
       for (devType in getOption("surveyPlot.dev"))
       {
           if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
           {
               width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6)
               dev <- openDev(device = devType,
                              directory = file.path(ResultsPath, "Stratified_(deprecated)"),
                              filename = c(paste(i, "_IDW", sep="")),
                              width = width,
                              height = width *
                                  (5 / 6) * diff(bboxMapZone[2, ]) / diff(bboxMapZone[1, ]),
                              pointsize = c(png = 24, jpg = 24, jpeg = 24, X11 = 11, pdf = 10),
                              counter = FALSE,
                              verbose=FALSE)
           }else{
               warning("Device \"", devType, "\" not supported")
               next()
           }

           layout(matrix(1:2, ncol=2), width=c(5, 1))
           oldpar <- par("mar", no.readonly=TRUE)
           par(mar=c(0, 0, 0, 0)+0.1, xpd=FALSE, oma = c(0, 0,
                                                         ifelse(getOption("surveyPlot.main"), 1, 0),
                                                         0.8) + 0.2)

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

           breaks <- breaksList[[i]]

           ## Plotting the interpolations:
           invisible(sapply(densInterp[[i]],
                            function(x)
                        {
                            image(x[paste(i, "pred", sep=".")],
                                  add=TRUE,
                                  col=getOption("surveyPlot.colIDW")(length(breaks) - 1),
                                  breaks=breaks)
                        }))

           logScale <- ((tolower(getOption("surveyIDW.log.scale")) == "default" &&
                        (is.element(i, names(getOption("surveyIDW.breaksList"))) ||
                         is.element(match.arg(getOption("surveyIDW.break.method"),
                                              c("quantile", "range", "logrange", "quantileStation")),
                                    c("logrange", "quantile")))) ||
                        isTRUE(getOption("surveyIDW.log.scale")))

           if (logScale &&
               any(breaks <= 0))
           {
               breakspos <- breaks + min(breaks, na.rm = TRUE) + 0.1 * diff(breaks[1:2])
           }else{
               breakspos <- breaks
           }

           if (isTRUE(getOption("surveyPlot.addStations")))
           {
               plot(stationData, add = TRUE,
                    col = getOption("surveyPlot.colStations"))
           }

           if (is.element(interpMetrics[i, "metrics"],
                          names(getOption("surveyPlot.legendTitle"))))
           {
               legendTitle <- getOption("surveyPlot.legendTitle")[interpMetrics[i, "metrics"]]
           }else{
               legendTitle <- interpMetrics[i, "metrics"]
           }

           if (getOption("surveyPlot.scaleBar"))
           {
               do.call(what = scale.bar, args = getOption("surveyPlot.scaleBar.opt"))
           }

           if (getOption("surveyPlot.northArrow"))
           {
               do.call(what = north.arrow, args = getOption("surveyPlot.northArrow.opt"))
           }

           par(mar=c(3, 0.2, 6, 3.8) + 0.1, las=1, xpd = NA, cex.axis = 0.8)
           image.scale(col = getOption("surveyPlot.colIDW")(length(breaks) - 1),
                       breaks = breakspos,
                       axis.pos=4, add.axis=FALSE,
                       title = legendTitle,
                       cex.title = 0.75, adj.title = 0,
                       log = ifelse(logScale, "y", ""))
           axis(side = 4, at = breakspos, labels = breaks, cex = 0.8)

           if (getOption("surveyPlot.main"))
           {
               mtext(text = bquote(italic(.(interpMetrics[i, "species"]))),
                     line = 0,
                     side = 3, outer = TRUE, adj = 0.5)
           }

           if (devType != "X11") dev.off()
       }

   })


sapply(names(densInterpPoly),
       function(i)
   {
       ## Graphical representation:
       for (devType in getOption("surveyPlot.dev"))
       {
           if (devType %in% c("X11", "pdf", "jpg", "jpeg", "png"))
           {
               width <- c(png = 1000, jpg = 1000, jpeg = 1000, X11 = 6.5, pdf = 6)
               dev <- openDev(device = devType,
                              directory = file.path(ResultsPath, "Stratified_(deprecated)"),
                              filename = c(paste(i, "_IDWcontours", sep="")),
                              width = width * 5 / 6,
                              height = width *
                                  diff(bboxMapZone[2, ]) / diff(bboxMapZone[1, ]),
                              pointsize = c(png = 36, jpg = 36, jpeg = 36, X11 = 11, pdf = 12),
                              counter = FALSE,
                              verbose=FALSE)
           }else{
               warning("Device \"", devType, "\" not supported")
               next()
           }

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

           breaks <- breaksList[[i]]

           ## Plotting the interpolations:
           plot(densInterpPoly[[i]],
                 add=TRUE,
                 col=getOption("surveyPlot.colCateg")(length(breaks) - 1)[densInterpPoly[[i]]$categ])

           if (is.element(interpMetrics[i, "metrics"],
                          names(getOption("surveyPlot.legendTitle"))))
           {
               legendTitle <- getOption("surveyPlot.legendTitle")[interpMetrics[i, "metrics"]]
           }else{
               legendTitle <- interpMetrics[i, "metrics"]
           }

           legend(x=getOption("surveyPlot.legendPos"),
                  legend=rev(levels(densInterpPoly[[i]]$categ)),
                  fill=rev(getOption("surveyPlot.colCateg")(length(breaks) - 1)),
                  title = legendTitle,
                  bg="white", inset=0.02, cex = 0.8)


           if (getOption("surveyPlot.scaleBar"))
           {
               do.call(what = scale.bar, args = getOption("surveyPlot.scaleBar.opt"))
           }

           if (getOption("surveyPlot.northArrow"))
           {
               do.call(what = north.arrow, args = getOption("surveyPlot.northArrow.opt"))
           }

           if (getOption("surveyPlot.main"))
           {
               mtext(text = bquote(italic(.(interpMetrics[i, "species"]))),
                     line = 0,
                     side = 3, outer = TRUE)
           }

           if (isTRUE(getOption("surveyPlot.addStations")))
           {
               plot(stationData, add = TRUE,
                    col = getOption("surveyPlot.colStations"))
           }

           if (devType != "X11") dev.off()
       }
   })



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
