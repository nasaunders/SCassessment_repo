#-*- coding: latin-1 -*-

### File: 0_Initialise_environment.R
### Time-stamp: <2019-04-08 18:06:10 yreecht>
###
### Created: 15/03/2016	17:22:56
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Load the functions and set defaults options if required.
####################################################################################################

## Storing the path of the scripts (for reliably find .geoData later notably):
options("survey.scriptPath" = getwd())

## Libraries directory (formerly functionsDir) stored in a
##   new var to avoid override by old analysis config scripts:
libDir <- normalizePath(file.path(scriptDir, "lib"))

options("survey.functionsPath" = libDir)

## Functions:
source(file.path(libDir, "Base_functions", "1_Functions_base.R"), encoding="latin1")
source(file.path(libDir, "Bayesian_functions", "1_Functions_Bayesian.R"), encoding="latin1")
source(file.path(libDir, "Geostat_functions", "1_Functions_geostats.R"), encoding="latin1")
source(file.path(libDir, "Geostat_functions", "1_Functions_hulls.R"), encoding="latin1")
source(file.path(libDir, "Graphical_functions", "1_Functions_graphics.R"), encoding="latin1")
source(file.path(libDir, "Graphical_functions", "1_Functions_maps.R"), encoding="latin1")
source(file.path(libDir, "Spatial_functions", "1_Functions_spatial.R"), encoding="latin1")



YGrBlPalette <- function(n)
{
    tmpYGrBlPalette <- colorRampPalette(colors = c("#FFFF80", "#8FF041", "#40AE3F",
                                                   "#2FA190", "#215896", "#0D1178"),
                                        interpolate = "spline")
    if (n < 4)
    {
        cols <- tmpYGrBlPalette(4)[seq_len(n)]
    }else{
        cols <- tmpYGrBlPalette(n)
    }

    return(cols)
}

## Defaults options:
defaults <-
    list(surveyTargetOnly = TRUE,       # Produce indices for target species only.
         surveyCoastMap = "Ireland",    # Any other value than Ireland => WorldMap
                                        #     (except if an alternative map is specified below).
         surveyMapSource = "shapefile", # Shapefile or a type for openmap (e.g. "bing" for satellite,
                                        #      see the help of the OpenStreetMap package).
         surveyBedForce = TRUE,         # Whether to force the extend of the bed to bedLayer.

         ## =============================
         ## Weight-length relationship:
         surveyWL.all.obs = TRUE,       # Use all observations for the weight-length relationship.
                                        # If FALSE, only stations in the assessed zone are used.

         ## =============================
         ## Survey zone:
         surveyIDW.alpha = 500,                # Max distance (in projection unit, usually m)
                                               #     between stations to avoid concave parts
                                               #     in the survey zone definition.
         surveyIDW.cellsize = 10,              # Cell size (in projection unit, usually m).
         surveyIDW.bufferWidth = 20,           # Buffer width (in projection unit, usually m) around the
                                               #     stations for the survey zone definition.

         ## =============================
         ## Inverse distance weighting:
         surveyIDW.nmin = 0,                   # Min number of neighbours. || (0 ->
         surveyIDW.nmax = Inf,                 # Max number of neighbours. ||  Inf => automatically defined).
         surveyIDW.idp = 2,                    # Inverse distance power (2 for an ArcMap like behaviour).
         surveyIDW.n.categ = 5,                # Number of breaks for the contour identification.
                                               #     Alternatively the breaks themselves can be given
                                               #     (only one metric recommended).
         surveyIDW.break.method = "quantileStation",
                                               # Method for calculating breaks.
                                               #     Choose among "quantileStation", "quantile",
                                               #     "range", "logrange"
         surveyIDW.breaksList = NULL,          # A named (by species) list of breaks.
                                               #     Override breaks calculations for the species present
                                               #     in the list.
         surveyIDW.log.scale = "default",      # Has the density values to be on a log scale?
                                               #     "default" will choose depending on the break method.

         ## =============================
         ## Plots:
         surveyPlot.colIDW = function(n) # Color palette for the idw interpolations.
         {
             rev(heat.colors(n))
         },
         surveyPlot.colCateg = YGrBlPalette,## function(n) # Color palette for the contours.
         ## {
         ##     bpy.colors(n, cutoff.tails=0.3)
         ## },
         surveyPlot.colGround = gray(0.96),  # Color for the ground.
         surveyPlot.colSea = "#d6edf8",      # Color for the sea.
         surveyPlot.colGrid = "grey",        # Color for the grid.
         surveyPlot.main = FALSE,            # Is the main title plotted.
         surveyPlot.addStations = TRUE,      # Are stations plotted?
         surveyPlot.colStations = "darkred", # Color for plotting stations.
         surveyPlot.dev = c("X11", "pdf"),   # Device(s) (in "X11", "pdf", "png", "jpg")
         surveyPlot.density.bw = 2,          # Band width for size distribution density representation.
         ## Density legend title by metrics:
         surveyPlot.legendTitle = c("AbundanceDensity" = expression("Density"~scriptstyle((italic(indiv.m^{-2})))),
                                    "AbundanceDensityCorr" = expression("Density"~scriptstyle((italic(indiv.m^{-2})))),
                                    "WeightDensity" = expression("Biomass dens."~scriptstyle((italic(kg.m^{-2})))),
                                    "WeightDensityCorr" = expression("Biomass dens."~scriptstyle((italic(kg.m^{-2}))))),
         surveyPlot.legendPos = "bottomright",            # Legend position.
         surveyPlot.scaleBar = TRUE,                      # Has a scalebar to be plotted?
         surveyPlot.scaleBar.opt = list(x = "bottomleft", # Scalebar parameters.
                                        inset = 0.02,
                                        ratio = 0.15),
         surveyPlot.northArrow = FALSE,                      # Has a north arrow to be plotted?
         surveyPlot.northArrow.opt = list(x = "topleft",    # North arrow parameters.
                                          inset = 0.02,
                                          ratio = 0.09),
         surveyPlot.colorBarSmall.opt = list(xspan = 0.05,
                                             ymin = 0.05,
                                             xmax = 0.88,
                                             legend.type = "value",
                                             nsign = 2),
         ## =============================
         ## Stats:
         surveyStat.alpha = 0.05,       # Alpha risk for IC calc.

         ## ==============================
         ## Bayesian stratified estimate:
         surveyBayes.n.burnin = 2000,
         surveyBayes.n.iter = 15000,
         surveyBayes.n.chains = 3,
         surveyBayes.models = c("N", "LN", "G"), # among c("N", "LN", "G")

         ## ==============================
         ## Geostatistical estimate options:
         surveyGeostat.propZero.max = 0.05,  # Max prop. of zeros to do one-step kriging.
         surveyGeostat.CGS.nsim = 500,       # Number of Conditional Gaussian Simulations.
         surveyGeostat.CGS.nmax = 30,        # Max number of neighbours for simulations (CGS).
                                        # Can be lowered to reduce computation time and memory.
         krigging.model = "~ 1",             # Covariate model (none by default).
         surveyGeostat.vg.cutoff = NULL,     # Cutoff for empirical variogram fitting/plots.
                                        # Automatically estimated from data if NULL.
         surveyGeostat.vg.width = NULL,      # break width for empirical variogram fitting/plots.
                                        # Automatically estimated from data if NULL.
         surveyGeostat.vg.autobreaks = TRUE, # Automatic calculation of breaks for the variogram
                                        # fitting/plot (based on irregular quantiles)
         surveyGeostat.vg.likelihood = TRUE, # Likelihood fitting
         surveyGeostat.vg.lik.method = "ML"  # One of "ML" (max likelihood), "REML" (restricted ML)
         )

## Index of undefined options:
idxNullOpt <- sapply(do.call(options, as.list(names(defaults))), is.null)

## Set the default for undefined options (if any):
if (any(idxNullOpt))
{
    do.call(options, defaults[idxNullOpt])

    message("\nSetting the following option(s) to default:")

    invisible(sapply(names(idxNullOpt)[idxNullOpt],
              function(x)
                 {
                     message(paste("    * ", x, " = ", defaults[x], sep = ""))
                 }))
}


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
