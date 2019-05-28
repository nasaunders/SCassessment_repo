#-*- coding: latin-1 -*-

### File: 4_Strata.R
### Time-stamp: <2018-11-07 15:32:01 yreecht>
###
### Created: 13/04/2016	11:06:58
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## Named list of breaks (calculated from quantiles) by interpolated field:
breaksList <- sapply(names(densInterp),
                     function(i, densInterp, interpMetrics)
                 {
                     ## ...Otherwise they are calculated from quantiles of interpolated values:
                     values <- do.call(c,
                                       sapply(densInterp[[i]],
                                              function(x)x[[paste(i, "pred", sep=".")]],
                                              simplify=FALSE))

                     if ( ! is.null(getOption("surveyIDW.breaksList")) &&
                         i %in% names(getOption("surveyIDW.breaksList")))
                     {
                         breaks <- getOption("surveyIDW.breaksList")[[i]]

                         if (is.infinite(tail(breaks, 1)))
                         {
                             breaks[length(breaks)] <- max(values, na.rm = TRUE)

                             if (diff(tail(breaks, 2)) <= 0) breaks <- head(breaks, -1)

                             n.digits <- ceiling(log10(1 / min(diff(tail(breaks, 2))))) + 1

                             breaks[length(breaks)] <- round(breaks[length(breaks)], digits = n.digits)

                             breaks[length(breaks)] <- ifelse(breaks[length(breaks)] < max(values, na.rm=TRUE),
                                                              breaks[length(breaks)] + 10^(-n.digits[length(n.digits)]),
                                                              breaks[length(breaks)])
                         }else{}
                     }else{
                         if (length(getOption("surveyIDW.n.categ")) - 1)
                         {
                             ## In case the breaks are forced by the user:
                             breaks <- getOption("surveyIDW.n.categ")
                         }else{
                             type <- match.arg(getOption("surveyIDW.break.method"),
                                               c("quantile", "range", "logrange", "quantileStation"))

                             if (type == "range")
                             {
                                 breaks <- prettyRangeBreaks(x = values,
                                                             n.categ = getOption("surveyIDW.n.categ"),
                                                             log = FALSE)
                             }else{}

                             if (type == "logrange")
                             {
                                 breaks <- prettyRangeBreaks(x = values,
                                                             n.categ = getOption("surveyIDW.n.categ"),
                                                             log = TRUE)
                             }else{}

                             if (type == "quantile")
                             {
                                 breaks <- prettyQuantileBreaks(x = values,
                                                                n.categ = getOption("surveyIDW.n.categ"))
                             }else{}

                             if (type == "quantileStation")
                             {
                                 values <- stationData@data[ , i]
                                 breaks <- prettyQuantileBreaks(x = values,
                                                                n.categ = getOption("surveyIDW.n.categ"))
                             }else{}

                             if (any(duplicated(breaks))) breaks <- unique(breaks)
                         }
                     }
                     return(breaks)
                 },
                     densInterp = densInterp, interpMetrics = interpMetrics,
                     simplify = FALSE)


## Calculate the polygons for every interpolated field:
densInterpPoly <- sapply(names(densInterp),
                         function(i, densInterp, breaksList)
                     {
                         ## Sometimes need for memory cleaning before conversion:
                         gc()

                         ## Conversion to polygones:
                         SP <- SPixDFlist2SPolyDF(densInterp[[i]],
                                                  variable=paste(i, "pred", sep="."),
                                                  breaks=breaksList[[i]])

                         ## Add polygon area:
                         SP@data <- cbind(SP@data, area=gArea(SP, byid=TRUE))

                         return(SP)
                     },
                         densInterp = densInterp, breaksList = breaksList,
                         simplify = FALSE)

## Save the Polygons in shapefiles:
invisible(sapply(names(densInterpPoly),
                 function(i, densInterpPoly, interpMetrics)
             {##browser()
                 tmpSP <- densInterpPoly[[i]]

                 ## Some fields may be arrays and need to be converted in vectors:
                 tmpSP@data[ , sapply(tmpSP@data, class) == "array"] <-
                     sapply(tmpSP@data[ , sapply(tmpSP@data, class) == "array"],
                            as.vector)

                 ## Shorter field names for export:
                 colnames(tmpSP@data) <- gsub(".", "_",
                                              sub(paste(interpMetrics[i, "field"], ".", sep=""),
                                                  "",
                                                  colnames(tmpSP@data)),
                                              fixed = TRUE)

                 tmpSP <- spTransform(x = tmpSP, CRSobj = CRS(projargs = projITM))

                 tryCatch(writeOGR(obj = tmpSP,
                                   dsn = normalizePath(file.path(ResultsPath, "shapeFiles"),
                                                       winslash = "/"),
                                   layer = paste(i, "_IDWcontours", sep=""),
                                   driver="ESRI Shapefile", overwrite_layer=TRUE, morphToESRI = TRUE),
                          error = function(e)
                 {
                     warning("Not able to write the .shp layer \"",
                             paste(i, "_IDWcontours", sep=""),
                             "\":\n\t",
                             e)
                 })

                 message("## File ",
                         file.path(ResultsPath, "shapeFiles",
                                   paste(i, "_IDWcontours.shp", sep="")),
                         " created.")
             },
                 densInterpPoly = densInterpPoly, interpMetrics = interpMetrics))

## Surveyed area:
surveyedArea <- spTransform(x = gBuffer(densInterpPoly[[1]],
                                       byid = FALSE,
                                       width = 0),
                            CRSobj = CRS(projITM))

surveyedArea <- SpatialPolygonsDataFrame(Sr = surveyedArea,
                                         data = data.frame(id = seq_len(length(surveyedArea)),
                                                           row.names = row.names(surveyedArea)),
                                         match.ID = TRUE)

tryCatch(writeOGR(obj = surveyedArea,
                  dsn = file.path(ResultsPath, "shapeFiles"),
                  layer = "Surveyed_area",
                  driver = "ESRI Shapefile", overwrite_layer = TRUE),
         error = function(e)
         {
             warning("Could not create ",
                     file.path(ResultsPath, "shapeFiles",
                               "Surveyed_area.shp:\n\t"),
                     e)
         })

## Contour polygons with the averaged metrics (area weighted):
densInterpPolyFusion <- sapply(X = names(densInterpPoly),
                               FUN = function(i, SP)
                           {
                               SP2 <- aggregate(SP[[i]][paste(i, "pred", sep = ".")],
                                                by=list(categ=SP[[i]]$categ),
                                                FUN = mean,
                                                dissolve=TRUE,
                                                areaWeighted=TRUE)

                               SP2@data <- cbind(SP2@data, "area"=gArea(SP2, byid=TRUE))

                               spChFIDs(SP2) <- SP2$categ

                               return(SP2)
                           },
                               SP = densInterpPoly,
                               simplify = FALSE)



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
