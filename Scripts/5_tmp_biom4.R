#-*- coding: latin-1 -*-

### File: 5_tmp.R
### Time-stamp: <2018-11-01 11:08:24 yreecht>
###
### Created: 20/07/2016	15:45:38
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

               stationsTmp <- stationData[idx, ]

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




head(stationsCat[[1]]@data)
head(biolCat[[1]])

statmp <- stationsCat[[1]]@data[ , c("WeightDensityCorr_Ensis_arcuatus", "stationID", "categDens")]
head(statmp)

bioltmp <- subset(biolCat[[1]],
                  ## LD_SpeciesID == "Cerastoderma edule")[ , c("Weight", "stationID", "categDens")]
                  LD_SpeciesID == "Ensis arcuatus")[ , c("Weight", "stationID", "categDens")]
head(bioltmp)


areas <- densInterpPolyFusion[[1]]@data

names(densInterpPolyFusion)

sum(areas$area)

X11()
par(mfcol = c(2, 2))

tapply(statmp$WeightDensityCorr_Ensis_arcuatus,
       statmp$categDens,
       hist, breaks = 10)

X11()
hist(statmp$WeightDensityCorr_Ensis_arcuatus, breaks = 10)

## meanweight <- by(bioltmp,
##                  bioltmp[ , c("categDens")],
##                  function(x)
##               {
##                   data.frame("categDens" = as.character(unique(x$categDens)),
##                              "Station_ID" = tapply(as.character(x$stationID),
##                                                   as.character(x$stationID),
##                                                   unique),
##                              "mWeight" = tapply(x$Weight,
##                                                 as.character(x$stationID),
##                                                 mean, na.rm = TRUE))
##               }, simplify = FALSE)

## meanweight <- do.call(rbind, meanweight)

## head(meanweight)
## summary(meanweight)

## X11()
## par(mfcol = c(3, 2))

## tapply(meanweight$mWeight,
##        meanweight$categDens,
##        hist, breaks = 20)

biomDens <- statmp
biomDens$biomDensity <- biomDens$WeightDensityCorr_Ensis_arcuatus

head(biomDens)

## X11()
## par(mfcol = c(2, 2))

## tapply(biomDens$biomDensity,
##        biomDens$categDens,
##        hist, breaks = 10)


listBdens <- tapply(biomDens$biomDensity,
                    biomDens$categDens,
                    function(x) x[! is.na(x)], simplify = FALSE)

listBdens[1]

maxStat <- max(sapply(listBdens, length))

matBdens <- do.call(cbind,
                    sapply(listBdens,
                           function(x, maxL)
                    {
                        if (maxL > length(x)) x <- c(x, rep(NA, maxL - length(x)))
                        return(x)
                    },
                           maxL = maxStat,
                           simplify = FALSE))

Nbc <- apply(matBdens, 2, function(x) sum(! is.na(x)))
head(matBdens, 29)

areascl <- areas$area

dataList <- list(matBdens = matBdens,## + 0.001,
                 Nbc = Nbc,
                 areascl = areascl)

library(rjags)

model1 <- jags.model(file = "../5_tmp.jag",
                     data = dataList, n.chains = 3, n.adapt = 10000)

samples <- coda.samples(model1,
                        variable.names = c("mu", "sigma", "Bcl",
                                           "Btot"),
                        n.iter = 15000)

dic1 <- dic.samples(model1, n.iter = 15000)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples[ , c("sigma[1]", "sigma[2]",
                      "sigma[3]", ## "sigma[4]",
                      "mu[1]", "mu[2]",  "mu[3]"## ,
                      ## "mu[4]"
                      )],
         show.obs=FALSE)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples, trace = FALSE, density = TRUE, show.obs=FALSE)
traceplot(samples[ , c("sigma[1]", "sigma[2]",
                       "sigma[3]", #"sigma[4]",
                       "mu[1]", "mu[2]",  "mu[3]"## ,
                       ## "mu[4]"
                       )])

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples[ , c("Bcl[1]", "Bcl[2]",
                      "Bcl[3]", ## "Bcl[4]",
                      "Btot")],
         show.obs=FALSE)


model2 <- jags.model(file = "../5_tmp2.jag",
                     data = dataList, n.chains = 3, n.adapt = 500)

samples2 <- coda.samples(model2,
                         variable.names = c("mu", "sigma", "Bcl",
                                            "Btot"),
                         n.iter = 5000)

dic2 <- dic.samples(model2, n.iter = 5000)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples2[ , c("sigma[1]", "sigma[2]",
                       "sigma[3]", ## "sigma[4]",
                       "mu[1]", "mu[2]",  "mu[3]"## ,
                       ## "mu[4]"
                       )],
         show.obs=FALSE)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples2[ , c("Bcl[1]", "Bcl[2]",
                       "Bcl[3]", ## "Bcl[4]",
                       "Btot")],
         show.obs=FALSE)

diffdic(dic1 = dic1, dic2 = dic2)

summary(samples2)$quan["Btot", c(3,1,5)] * 1000
diff(summary(samples2)$quan["Btot", c(1,3)]) * 1000
diff(summary(samples2)$quan["Btot", c(3,5)]) * 1000

summary(samples)$quan["Btot", c(3,1,5)]
diff(summary(samples)$quan["Btot", c(1,3)])
diff(summary(samples)$quan["Btot", c(3,5)])


## ##########################################################################################
##

matBdens2 <- cbind(biomDens$biomDensity, NULL)
colnames(matBdens2) <- "Tot"

Nbc2 <- c(Tot = length(biomDens$biomDensity))
head(matBdens2, 5)

areascl2 <- sum(areas$area)

dataList2 <- list(matBdens = matBdens2 + 0.001,
                  Nbc = Nbc2,
                  areascl = areascl2)


model1 <- jags.model(file = "../5_tmp.jag",
                     data = dataList2, n.chains = 3, n.adapt = 10000)

samples <- coda.samples(model1,
                        variable.names = c("mu", "sigma", "Bcl",
                                           "Btot"),
                        n.iter = 15000)

dic1.2 <- dic.samples(model1, n.iter = 15000)

X11()
par(mfrow=c(2, 2))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples[ , c("sigma", "mu",
                       "Bcl", "Btot")],
         show.obs=FALSE)


model2 <- jags.model(file = "../5_tmp2.jag",
                     data = dataList2, n.chains = 3, n.adapt = 500)

samples2 <- coda.samples(model2,
                         variable.names = c("mu", "sigma", "Bcl",
                                            "Btot"),
                         n.iter = 5000)

dic2.2 <- dic.samples(model2, n.iter = 5000)

X11()
par(mfrow=c(2, 2))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples2[ , c("sigma", "mu",
                       "Bcl", "Btot")],
         show.obs=FALSE)

dic1.2 - dic2.2

summary(samples)$quan["Btot", c(3,1,5)] * 1000
diff(summary(samples)$quan["Btot", c(1,3)]) * 1000
diff(summary(samples)$quan["Btot", c(3,5)]) * 1000

summary(samples2)$quan["Btot", c(3,1,5)] * 1000
diff(summary(samples2)$quan["Btot", c(1,3)]) * 1000
diff(summary(samples2)$quan["Btot", c(3,5)]) * 1000

diffdic(dic2.2, dic2)

str(dic2)
str(dic2.2)

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
