#-*- coding: latin-1 -*-

### File: 1_Load_spatial_data.R
### Time-stamp: <2019-01-31 12:40:33 yreecht>
###
### Created: 14/03/2016	15:37:04
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Generic script for loading data (survey specific calculation in a separate script).
####################################################################################################

## Normalize, remove final "/":
if (exists("shpOtherDir"))
{
    shpOtherDir <- file.path(shpOtherDir)
}
shpDataDir <- file.path(shpDataDir)

if (is.null(getOption("surveyCoastMap"))) options(surveyCoastMap = "Ireland")

## -------------------------------------------------------------------------------------------------
## Usual projections:
.EPSG <- make_EPSG()

projWGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
projITM <- .EPSG[.EPSG$code == 2157 & ! is.na(.EPSG$code), "prj4"] # Irish Transverse Mercator.
## -------------------------------------------------------------------------------------------------

## Loading shapefiles:
stationData <- readOGR(dsn=shpDataDir, layer=stationLayer, encoding="utf8")

## Renaming fields if necessary:
if ( ! any(grepl("TRACKID", colnames(stationData@data))))
{
    colnames(stationData@data)[which(toupper(colnames(stationData@data)) %in%
                                     c("TRACKID", "TRACK_ID"))] <- "TRACKID"
}

if ( ! any(grepl("TIMESTART", colnames(stationData@data))))
{
    colnames(stationData@data)[which(toupper(colnames(stationData@data)) %in%
                                     c("TIMESTART", "TIME_START",
                                       "START_TIME", "STARTTIME"))] <- "TIMESTART"
}

source(file.path(getOption("survey.scriptPath"), "1_Check_spatial_data.R"),
       encoding="latin1")

source(file.path(getOption("survey.scriptPath"), "1_Summary_spatial_data.R"),
       encoding="latin1")

stationDataTot <- stationData           # For station summary (all tracks).

## Selection mask (if any):
if ( ! exists("maskLayer") || is.null(maskLayer) || maskLayer == "")
{
    mask <- NULL
}else{
    mask <- readOGR(dsn=shpDataDir, layer=maskLayer, encoding="utf8")

    mask <- spTransform(x = mask, CRSobj = CRS(proj4string(stationData)))

    stationData <- stationData[as.vector(gIntersects(stationData,
                                                     gBuffer(mask, width = 0),
                                                     byid = TRUE)) , ]

    ## Saving the zone as potential krigging covariate:
    stationData[["zone"]] <- factor(over(stationData, mask)[ , 1])
}

## Zones for plotting maps by survey:
surveyMapZones <- readOGR(dsn = file.path(getOption("survey.scriptPath"), ".geoData"),
                          layer = "SurveyMapZones", encoding = "utf8")

surveyMapZones <- spTransform(surveyMapZones, CRSobj=CRS(proj4string(stationData)))

## Identification of the proper map zone:
mapZone <- findPlotArea(stationData = stationData,
                        surveyMapZones = surveyMapZones,
                        addStationExtent=TRUE)

bboxMapZone <- bbox(obj = mapZone)

## Base map:
if ((! exists("coastlineLayer")) || is.null(coastlineLayer))  #
{
    mapZoneWGS84 <- spTransform(mapZone, CRSobj=projWGS84)

    ## Extending  the bounding box by 10% in every direction:
    bboxArea <- bbox(mapZoneWGS84) +
        sweep(x = 0.1 * matrix(rep(c(-1, 1), each = 2), ncol = 2),
              MARGIN = 1,
              STAT = apply(X = bbox(mapZoneWGS84),
                           MARGIN = 1,
                           FUN = function(x) diff(range(x))),
              FUN = "*")

    ## bboxArea <- sweep(x = bbox(mapZoneWGS84), MARGIN = 2, STATS = c(-0.1, +0.1), FUN="+")

    if (getOption("surveyMapSource")[1] == "shapefile")
    {
        if (getOption("surveyCoastMap") %in% c("Ireland", "ireland"))
        {
            coastline <- readOGR(dsn = file.path(getOption("survey.scriptPath"), ".geoData"),
                                 layer = "Ireland", encoding = "utf8")


            ## Cropping areas for the map (very large):
            cropArea <- as(extent(bboxArea["x", ],
                                  bboxArea["y", ]),
                           "SpatialPolygons")

            proj4string(cropArea) <- proj4string(mapZoneWGS84)

            cropAreaProj <- spTransform(x = cropArea, CRSobj = CRS(proj4string(coastline)))

            ## Crop the map:
            coastline <- gIntersection(coastline, cropAreaProj, byid = TRUE)

        }else{
            baseMap <- map(database="worldHires",
                           regions=".",                    # select any polygon within the long/lat range.
                           fill = TRUE, col = gray(c(0.6)),
                           xlim=bboxArea["x", ],
                           ylim=bboxArea["y", ],
                           interior=TRUE, plot=FALSE)

            coastline <- map2SpatialPolygons(baseMap,
                                             proj4string = projWGS84,
                                             IDs = sapply(strsplit(baseMap$names, ":"),
                                                          function(x) x[1]))

            coastline <- SpatialPolygonsDataFrame(Sr=coastline,
                                                  data=data.frame(ID=sapply(coastline@polygons,
                                                                            slot, name="ID"),
                                                                  row.names=sapply(coastline@polygons,
                                                                                   slot, name="ID")))
        }

        coastline <- spTransform(coastline, CRSobj=CRS(proj4string(stationData)))
    }else{
        ## require("OpenStreetMap")

        avTypes <- c("osm", "osm-bw", "maptoolkit-topo", "waze", "bing", "stamen-toner",
                     "stamen-terrain", "stamen-watercolor", "osm-german", "osm-wanderreitkarte",
                     "mapbox", "esri", "esri-topo", "nps", "apple-iphoto", "skobbler",
                     "hillshade", "opencyclemap", "osm-transport", "osm-public-transport",
                     "osm-bbike", "osm-bbike-german")


        mapType <- match.arg(getOption("surveyMapSource")[1],
                             avTypes)

        message("## Downloading the basemap... it may take a while!")
        coastline <-  tryCatch(openmap(upperLeft = c(bboxArea[2, 2], bboxArea[1, 1]),
                                       lowerRight = c(bboxArea[2, 1], bboxArea[1, 2]),
                                       type=mapType, minNumTiles=12, mergeTiles=FALSE),
                               error = function(e)
                      {
                          tryCatch(openmap(upperLeft = c(bboxArea[2, 2], bboxArea[1, 1]),
                                           lowerRight = c(bboxArea[2, 1], bboxArea[1, 2]),
                                           type=mapType, minNumTiles=10, mergeTiles=FALSE),
                                   error = function(e)
                          {
                              tryCatch(openmap(upperLeft = c(bboxArea[2, 2], bboxArea[1, 1]),
                                               lowerRight = c(bboxArea[2, 1], bboxArea[1, 2]),
                                               type=mapType, minNumTiles=9, mergeTiles=FALSE),
                                       error = function(e)
                              {
                                  tryCatch(openmap(upperLeft = c(bboxArea[2, 2], bboxArea[1, 1]),
                                                   lowerRight = c(bboxArea[2, 1], bboxArea[1, 2]),
                                                   type=mapType, minNumTiles=4, mergeTiles=FALSE),
                                           error = function(e)
                                  {
                                      openmap(upperLeft = c(bboxArea[2, 2], bboxArea[1, 1]),
                                              lowerRight = c(bboxArea[2, 1], bboxArea[1, 2]),
                                              type=mapType, minNumTiles=2, mergeTiles=FALSE)
                                  })
                              })
                          })
                      })

        coastline <- openproj(x = coastline, projection = proj4string(stationData))

        coastline <- raster(x = coastline)
        ## plot(coastline)
        ## plot.OpenStreetMap(coastline)

        coastline <- as(coastline, "SpatialGridDataFrame")

        ## require(ggmap)
        ## coastline <- get_map(location = c(left = bboxArea[1, 1], bottom = bboxArea[2, 1],
        ##                                   right = bboxArea[1, 2], top = bboxArea[2, 2]),
        ##                      maptype = "satellite",
        ##                      source = "google")

        ## plot(coastline)
        ## class(coastline)
        ## str(coastline)



        ## test <- coastline
        ## attr(test, "class") <- "raster"
        ## test <- as(test, "SpatialGridDataFrame")
    }

}else{
    if (grepl("\\.tiff?$", coastlineLayer))
    {
        ## coastline <- raster(file.path(shpOtherDir, coastlineLayer))
        coastline <- readGDAL(file.path(shpOtherDir, coastlineLayer))
        ## coastline <- spTransform(coastline, CRSobj=CRS(proj4string(stationData)))
    }else{
        coastline <- readOGR(dsn=shpOtherDir, layer=coastlineLayer, encoding="utf8")
        coastline <- spTransform(coastline, CRSobj=CRS(proj4string(stationData)))
    }
}

if ((exists("gridLayer")) && (! is.null(gridLayer)))
{
    surveyGrid <- readOGR(dsn=shpOtherDir, layer=gridLayer, encoding="utf8")
    surveyGrid <- spTransform(surveyGrid, CRSobj=CRS(proj4string(stationData)))
}

if (exists("unsuitableLayer") && (! is.null(unsuitableLayer)))
{
    unsuitableGround <- readOGR(dsn=shpOtherDir, layer=unsuitableLayer, encoding="utf8")
    unsuitableGround <- gBuffer(spTransform(unsuitableGround, CRSobj=CRS(proj4string(stationData))),
                                width = 0)

    ## High water line can be added to unsuitable grounds:
    if (is.element(class(coastline),
                   c("SpatialPolygons", "SpatialPolygonsDataFrame")))
    {
        unsuitableGround <- tryCatch(gUnion(spgeom1 = gBuffer(unsuitableGround, width = 0),
                                            spgeom2 = gBuffer(coastline, width = 0),
                                            byid=FALSE, drop_lower_td=TRUE,
                                            unaryUnion_if_byid_false=TRUE, checkValidity=TRUE),
                                     error = function(e)
                            {
                                warning("Impossible to add the land to unsuitable grounds.")
                                return(unsuitableGround)
                            })
    }
}else{
    unsuitableGround <- NULL
}

if (exists("bedLayer") && (! is.null(bedLayer)))
{
    suitableGround <- readOGR(dsn=shpOtherDir, layer=bedLayer, encoding="utf8")
    suitableGround <- gBuffer(spTransform(suitableGround, CRSobj=CRS(proj4string(stationData))),
                              byid = TRUE,
                              width = 0)
}else{
    suitableGround <- NULL
}

class(suitableGround)

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
