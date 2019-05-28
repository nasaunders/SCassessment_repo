#-*- coding: latin-1 -*-

### File: 5_tmp.R
### Time-stamp: <2018-11-01 11:02:06 yreecht>
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

statmp <- stationsCat[[1]]@data[ , c("AbundanceDensity_Cerastoderma_edule", "Station_ID", "categDens")]
head(statmp)

bioltmp <- subset(biolCat[[1]],
                  LD_SpeciesID == "Cerastoderma edule")[ , c("Weight", "stationID", "categDens")]
head(bioltmp)


areas <- densInterpPolyFusion[[1]]@data

sum(areas$area)

X11()
par(mfcol = c(3, 2))

tapply(statmp$AbundanceDensity_Cerastoderma_edule,
       statmp$categDens,
       hist)

X11()
par(mfcol = c(3, 2))

tapply(bioltmp$Weight,
       bioltmp$categDens,
       hist)

meanweight <- by(bioltmp,
                 bioltmp[ , c("categDens")],
                 function(x)
              {
                  data.frame("categDens" = as.character(unique(x$categDens)),
                             "Station_ID" = tapply(as.character(x$stationID),
                                                  as.character(x$stationID),
                                                  unique),
                             "mWeight" = tapply(x$Weight,
                                                as.character(x$stationID),
                                                mean, na.rm = TRUE))
              }, simplify = FALSE)

meanweight <- do.call(rbind, meanweight)

head(meanweight)
summary(meanweight)

X11()
par(mfcol = c(3, 2))

tapply(meanweight$mWeight,
       meanweight$categDens,
       hist, breaks = 20)

biomDens <- merge(statmp, meanweight, all = TRUE)
biomDens$biomDensity <- biomDens$AbundanceDensity_Cerastoderma_edule * biomDens$mWeight

head(biomDens)

X11()
par(mfcol = c(3, 2))

tapply(biomDens$biomDensity,
       biomDens$categDens,
       hist, breaks = 10)

head(biomDens)

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

dataList <- list(matBdens = matBdens,
                 Nbc = Nbc,
                 areascl = areascl)

library(rjags)

model1 <- jags.model(file = "../5_tmp.jag",
                     data = dataList, n.chains = 3, n.adapt = 500)

samples <- coda.samples(model1,
                        variable.names = c("mu", "sigma", "Bcl",
                                           "Btot"),
                        n.iter = 5000)

dic1 <- dic.samples(model1, n.iter = 5000)

X11()
par(mfrow=c(4, 3))
## densplot(...)
## plot(samples, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples[ , c("sigma[1]", "sigma[2]",
                      "sigma[3]", "sigma[5]",
                      "sigma[5]", "mu[1]",
                       "mu[2]",  "mu[3]",
                       "mu[4]",  "mu[5]")],
         show.obs=FALSE)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples[ , c("Bcl[1]", "Bcl[2]",
                      "Bcl[3]", "Bcl[4]",
                      "Bcl[5]", "Btot")],
         show.obs=FALSE)


model2 <- jags.model(file = "../5_tmp2.jag",
                     data = dataList, n.chains = 3, n.adapt = 500)

samples2 <- coda.samples(model2,
                         variable.names = c("mu", "sigma", "Bcl",
                                            "Btot"),
                         n.iter = 5000)

dic2 <- dic.samples(model2, n.iter = 5000)

X11()
par(mfrow=c(4, 3))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples2[ , c("sigma[1]", "sigma[2]",
                       "sigma[3]", "sigma[5]",
                       "sigma[5]", "mu[1]",
                       "mu[2]",  "mu[3]",
                       "mu[4]",  "mu[5]")],
         show.obs=FALSE)

X11()
par(mfrow=c(3, 2))
## densplot(...)
## plot(samples2, trace = FALSE, density = TRUE, show.obs=FALSE)
densplot(samples2[ , c("Bcl[1]", "Bcl[2]",
                       "Bcl[3]", "Bcl[4]",
                       "Bcl[5]", "Btot")],
         show.obs=FALSE)

diffdic(dic1 = dic1, dic2 = dic2)

summary(samples2)$quan["Btot", c(3,1,5)]
summary(samples)$quan["Btot", c(3,1,5)]


fac <- c(1, 2, 2, 1, 3, 4, 1)
x <- c(1, 1, 2, NA, 4, NA, 3)

tapply(x, fac, sum, na.rm = TRUE)

sum(NA, na.rm = TRUE)

biomassCat$AbundanceDensity_Ensis_arcuatus$categDens

i <- "AbundanceDensity_Ensis_arcuatus"
funcInit()



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
