#-*- coding: latin-1 -*-

### File: 5_tmp_jags_ext.R
### Time-stamp: <2018-11-01 11:02:06 yreecht>
###
### Created: 07/12/2016	17:34:18
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## #####################################################################
## dens jags

## debuggingState(on=FALSE)
weightDensCat <- sapply(X = names(densInterpPolyFusion),
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

                      ## browser()

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
                                                                      mean, na.rm = TRUE),
                                                   "nWeight" = tapply(x$Weight,
                                                                      as.character(x$stationID),
                                                                      function(x) sum(! is.na(x))))
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
                      return(NULL)
                  }

                  head(biomDens)
                  dim(biomDens)

                  return(biomDens)
              },
                  stationsCat = stationsCat, biolCat = biolCat, metrics = metrics,
                  simplify=FALSE, USE.NAMES=TRUE)

weightDensCat <- weightDensCat[ ! sapply(weightDensCat, is.null)]

tapply(weightDensCat[[1]]$mWeight, weightDensCat[[1]]$categDens, mean, na.rm = TRUE)

head(weightDensCat[[1]])

i <- names(weightDensCat)[1]

biomassData <- weightDensCat[[i]]
head(biomassData)

areas <- densInterpPolyFusion[[i]]@data

listDens <- tapply(biomassData[ , i],
                   biomassData$categDens,
                   function(x) x[! is.na(x)], simplify = FALSE)

## Max number of stations by strata:
maxStat <- max(sapply(listDens, length))

## Matrix of biomass densities:
matDens <- do.call(cbind,
                   sapply(listDens,
                          function(x, maxL)
                   {
                       if (maxL > length(x)) x <- c(x, rep(NA, maxL - length(x)))
                       return(x)
                   },
                   maxL = maxStat,
                   simplify = FALSE))

## meanDens <- apply(matDens, 2, mean, na.rm = TRUE)

listWbar <- tapply(biomassData$mWeight,
                   biomassData$categDens,
                   function(x) x / 1e3, simplify = FALSE) # kg

## Matrix of biomass densities:
matWbar <- do.call(cbind,
                   sapply(listWbar,
                          function(x, maxL)
                   {
                       x <- x
                       if (maxL > length(x)) x <- c(x, rep(NA, maxL - length(x)))
                       return(x)
                   },
                   maxL = maxStat,
                   simplify = FALSE))

listNWbar <- tapply(biomassData$nWeight,
                    biomassData$categDens,
                    function(x)
             {
                 x[is.na(x)] <- 100
                 return(x)
             }, simplify = FALSE)

## Matrix of biomass densities:
matNWbar <- do.call(cbind,
                    sapply(listNWbar,
                           function(x, maxL)
                    {
                        if (maxL > length(x)) x <- c(x, rep(NA, maxL - length(x)))
                        return(x)
                    },
                    maxL = maxStat,
                    simplify = FALSE))



## Number of stations by strata:
Nbc <- apply(matDens, 2, function(x) sum(! is.na(x)))
head(matDens, 29)

## Surface area
areascl <- areas$area

require(rjags)
require(R2jags)

dataList <- list(matDens = matDens +
                     ifelse(any(na.omit(as.vector(matDens)) == 0),
                            10^(floor(log10(min(matDens[matDens != 0], na.rm = TRUE))) - 2),
                            0),
                 matWbar = matWbar +
                     ifelse(any(na.omit(as.vector(matWbar)) == 0),
                            10^(floor(log10(min(matWbar[matWbar != 0], na.rm = TRUE))) - 2),
                            0),
                 nw = matNWbar,
                 ## minDens = apply(matDens, 2, min, na.rm = TRUE) +
                 ##     ifelse(any(na.omit(as.vector(matDens)) == 0),
                 ##            0.001, 0),
                 ## maxDens = apply(matDens, 2, max, na.rm = TRUE) +
                 ##     ifelse(any(na.omit(as.vector(matDens)) == 0),
                 ##            0.001, 0),
                 meanDens = apply(matDens, 2, mean, na.rm = TRUE) +
                     ifelse(any(na.omit(as.vector(matDens)) == 0),
                            10^(floor(log10(min(matDens[matDens != 0], na.rm = TRUE))) - 2),
                            0),
                 meanW = apply(matWbar, 2, mean, na.rm = TRUE) +
                     ifelse(any(na.omit(as.vector(matWbar)) == 0),
                            10^(floor(log10(min(matWbar[matWbar != 0], na.rm = TRUE))) - 2),
                            0),
                 Nbc = Nbc,
                 areascl = areascl)

sd0 <- which(apply(dataList$matDens, 2, sd, na.rm = TRUE) == 0)

if (length(sd0))
{
    matDensTmp <- dataList$matDens

    matDensTmp[ , sd0] <- matDensTmp[ , sd0, drop = FALSE] +
        sapply(sd0,
               function(x, data)
        {
            res <- rnorm(n = nrow(data),
                         mean = 0,
                         sd = data[1, x] * 1e-03)
        },
        data = matDensTmp)

    dataList$matDens <- matDensTmp
}

funcInit <- function()
{
    res <- list("muDens" = runif(n = ncol(dataList$matDens),
                                 min = apply(dataList$matDens, 2, min, na.rm = TRUE) * 1,
                                 max = apply(dataList$matDens, 2, max, na.rm = TRUE) * 1),
                "sigmaDens" = runif(n = ncol(dataList$matDens),
                                    min = apply(dataList$matDens, 2, sd, na.rm = TRUE) * 0.8,
                                    max = apply(dataList$matDens, 2, sd, na.rm = TRUE) * 1.2),
                "muWbar" = runif(n = ncol(dataList$matWbar),
                                 min = apply(dataList$matWbar, 2, min, na.rm = TRUE) * 1,
                                 max = apply(dataList$matWbar, 2, max, na.rm = TRUE) * 1),
                "sigmaWbar" = runif(n = ncol(dataList$matWbar),
                                    min = apply(dataList$matWbar, 2, sd, na.rm = TRUE) * 0.8,
                                    max = apply(dataList$matWbar, 2, sd, na.rm = TRUE) * 1.2))

    res$sigmaDens[is.na(res$sigmaDens)] <- mean(res$sigmaDens, na.rm = TRUE)
    res$sigmaWbar[is.na(res$sigmaWbar)] <- mean(res$sigmaWbar, na.rm = TRUE)

    return(res)
}

funcInit()

testN <- jags(model.file = "../5_Biomass_N_tmp.jag",
              data = dataList, n.chains = 3, n.iter = 45000,
              n.burnin = 10000,
              inits = funcInit,
              n.thin = 1,
              parameters.to.save = c("muDens", "sigmaDens",
                                     "muWbar", "sigmaWbar",
                                     ## "Dens",
                                     "matWbar",
                                     "Bcl", "Bdens",
                                     "BDtot", "Btot"))

testN

testLN <- jags(model.file = "../5_Biomass_LN_tmp.jag",
               data = dataList, n.chains = 3, n.iter = 45000,
               n.burnin = 15000,
               inits = funcInit,
               n.thin = 1,
               parameters.to.save = c("muDens", "sigmaDens",
                                      "muWbar", "sigmaWbar",
                                      ## "Dens",
                                      "matWbar",
                                      "Bcl", "Bdens",
                                      "BDtot", "Btot"))

testLN








### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
