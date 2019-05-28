#-*- coding: latin-9 -*-

### File: 2_Data_preparation_QC.R
### Time-stamp: <2019-02-04 17:04:21 yreecht>
###
### Created: 12/09/2018	12:25:53
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

message("\n## Quality check:\n")
library(coda)

options("surveyPlot.dev" = "png")
ResultsPath <- outputDir


if (exists("weighLengthFile") &&
    ! is.null(weighLengthFile) &&
    file.exists(weighLengthFile))
{
    LWdata <- read.csv(file = weighLengthFile)

    colnames(LWdata)[1:3] <- c("species", "size", "weight") # may have further columns (if subset wanted upstream).
}else{
    nbMeas <- sapply(by(data = measurements,
                        INDICES = measurements$species,
                        FUN = function(x)
                        {
                            sum( (! is.na(x$length)) & (! is.na(x$weight)))
                        }, simplify=TRUE),
                     I)

    if (any(nbMeas >= 30))
    {
        LWdata <- subset(x = measurements,
                     subset = species %in% names(nbMeas)[nbMeas >= 30] &
                         ! is.na(length) & ! is.na(weight),
                     select = c("species", "length", "weight"))

        colnames(LWdata)[1:3] <- c("species", "size", "weight")
    }else{
        LWdata <- NULL
    }
}

if ( ! is.null(LWdata))
{
    ## Modelling the L-W relationships:
    .glmRes <- list()
    source(file.path(scriptPath, "2_Weight-length_relationships.R"))

    colBy <- c("species", "grade", "stationid", "gear", "vessel.team")

    colBy <- colBy[! sapply(measurements, function(x)all(is.na(x)))[colBy]]

    ## if (getOption("hasWeight"))
    ## {
    ## Simulation of sample weights from the L-W relationships:
    Nit <- nrow(unique(measurements[ , colBy]))

    pb <- init.progressBar(width = 50, max = Nit)

    resWLcatch <- by(data = measurements,
                     INDICES = as.list(measurements[ , colBy]),
                     FUN = function(x, glmRes, LWparam, pb)
                     {
                         Nsimu <- 3000
                         sp <- as.character(x$species[1])

                         stepProgressBar(pb)## ,
                         ## extra = paste0(" (",
                         ##                paste(apply(unique(x[c("species",
                         ##                                       "grade", "stationid")]),
                         ##                            1, as.character),
                         ##                      collapse = "+"),
                         ##                ")"))

                         if (sp %in% names(glmRes))
                         {
                             ## field "size" required for the GLM predictions:
                             x$size <- x$length

                             calcWgMean <- sum(LWparam[sp , "a"] *
                                               x$length^LWparam[sp , "b"] * 1e-3,
                                               na.rm = TRUE)

                             calcWgSimu <- sapply(1:Nsimu,
                                                  function(i, x, glmRes)
                                                  {
                                                      exp(predict(glmRes$glmObj, newdata = x) +
                                                          rnorm(n = nrow(x),
                                                                mean = 0,
                                                                sd = glmRes$glmPred$residual.scale))
                                                  },
                                                  x = x, glmRes = glmRes[[sp]])

                             calcWgSimuTot <- if(is.matrix(calcWgSimu))
                                              {
                                                  apply(calcWgSimu, 2, sum, na.rm = TRUE)
                                              }else{         # If only one indiv, no need for a sum.
                                                  calcWgSimu
                                              }

                             meanSimu <- mean(calcWgSimuTot * 1e-3, na.rm = TRUE)

                             CI95 <- as.data.frame(HPDinterval(as.mcmc(calcWgSimuTot), prob = 0.95) * 1e-3)
                             colnames(CI95) <- paste0("CI95%", colnames(CI95))

                             Nmeas <- sum(! is.na(x$length))

                             return(cbind(unique(x[ , colBy]),
                                          Nmeas = Nmeas,
                                          meanGLM = calcWgMean,
                                          meanSimu = meanSimu,
                                          CI95))
                         }else{
                             return(NULL)
                         }
                     }, glmRes = .glmRes, LWparam = LWparamCalc, pb = pb, simplify=FALSE)

    resWLcatch <- do.call(rbind, resWLcatch)

    samplesQC <- merge(samples, resWLcatch, all.x = TRUE)

    ## ####################################################################################################
    ## mean indiv weight by species:
    meanIndWg <- by(samplesQC, as.list(samplesQC[ , c("species", "stationid", "grade")]),
                    function(x) x[ , "meanGLM"] / x[ , "Nmeas"])


    Nind <- tapply(biological$Size,
                   as.list(biological[ , c("LD_SpeciesID", "LD_HaulNo", "grade")]),
                   function(x)sum(! is.na(x)))


    meanIndWgRed <- meanIndWg[dimnames(Nind)[[1]], dimnames(Nind)[[2]], dimnames(Nind)[[3]], drop = FALSE]


    meanIndWeightSp <- sapply(dimnames(meanIndWgRed)[[1]],
                              function(sp) weighted.mean(meanIndWgRed[sp, , ],
                                                         Nind[sp, , ],
                                                         na.rm = TRUE))

    ## ####################################################################################################
    ## Correction of simulated sample weights for broken razors:
    ##   if only broken (no measurable individual),
    ##   the global specific mean individual weight is used,
    ##   with a 0.1-0.9 interval for individual weight loss.
    samplesQC$sWg.GLM.corr <- ifelse(test = is.na(samplesQC$meanGLM) & is.na(samplesQC$Nmeas),
                                     yes = meanIndWeightSp[samplesQC$species] * 0.5 * samplesQC$broken,
                                     no = samplesQC$meanGLM * (1 + 0.5 * samplesQC$broken / samplesQC$Nmeas))
    samplesQC$sWg.simu.corr <- ifelse(test = is.na(samplesQC$meanSimu) & is.na(samplesQC$Nmeas),
                                      yes = meanIndWeightSp[samplesQC$species] * 0.5 * samplesQC$broken,
                                      no = samplesQC$meanSimu * (1 + 0.5 * samplesQC$broken / samplesQC$Nmeas))
    ## For CI, add between 10 and 90% of mean individual weight for broken:
    samplesQC$"WL_CI95%lower" <- ifelse(test = is.na(samplesQC$"CI95%lower") & is.na(samplesQC$Nmeas),
                                        yes = meanIndWeightSp[samplesQC$species] * 0.1 * samplesQC$broken,
                                        no = samplesQC$"CI95%lower" * (1 + 0.1 * samplesQC$broken / samplesQC$Nmeas))
    samplesQC$"WL_CI95%upper" <- ifelse(test = is.na(samplesQC$"CI95%upper") & is.na(samplesQC$Nmeas),
                                        yes = meanIndWeightSp[samplesQC$species] * 0.9 * samplesQC$broken,
                                        no = samplesQC$"CI95%upper" * (1 + 0.9 * samplesQC$broken / samplesQC$Nmeas))


    ## head(biological)
    ## head(samplesQC)
    ## head(subset(samplesQC, broken == 0 & is.na(Nmeas) & species == "Ensis siliqua"))
    ## samplesQC[which(is.na(samplesQC$"CI95%lower") & samplesQC$broken > 0), ]

    samplesQC$consistentWg <- samplesQC$sampleweight >= samplesQC$"WL_CI95%lower" &
        samplesQC$sampleweight <= samplesQC$"WL_CI95%upper"


    ## ################################################################################
    ## Formatting the results:
    extraCols <- c(if (length(unique(samplesQC$gear)) > 1)
                   {"gear"},
                   if (length(unique(samplesQC$vessel.team)) > 1)
                   {"vessel.team"})

    samplesQC.LW <- subset(samplesQC,
                           ! consistentWg,
                           select = c("species", "grade", "stationid",
                                      extraCols,
                                      "sampleweight", "WL_CI95%lower", "WL_CI95%upper"))

    samplesQC.LW$diff <- ifelse(test = samplesQC.LW$sampleweight < samplesQC.LW$"WL_CI95%lower",
                                yes = paste0(round((samplesQC.LW$sampleweight -
                                                    samplesQC.LW$"WL_CI95%lower") * 100 /
                                                   samplesQC.LW$"WL_CI95%lower"),
                                             "%"),
                                no = ifelse(test = samplesQC.LW$sampleweight > samplesQC.LW$"WL_CI95%upper",
                                            yes = paste0("+",
                                                         round((samplesQC.LW$sampleweight -
                                                                samplesQC.LW$"WL_CI95%upper") * 100 /
                                                               samplesQC.LW$"WL_CI95%upper"),
                                                         "%"),
                                            no = ""))


    samplesQC.LW <- samplesQC.LW[order(as.numeric(gsub("[-+]([[:digit:]]+)[%]",
                                                       "\\1",
                                                       samplesQC.LW$diff)),
                                       decreasing = TRUE),
                                 ]

    colnames(samplesQC.LW)[match(c("stationid", "sampleweight",
                                   "WL_CI95%lower", "WL_CI95%upper"),
                                 colnames(samplesQC.LW))] <-
        c("station", "sampleWg",
          "WL_CI95%low", "WL_CI95%up")

    samplesQC.LW$"WL_CI95%low" <- round(samplesQC.LW$"WL_CI95%low", 3)
    samplesQC.LW$"WL_CI95%up" <- round(samplesQC.LW$"WL_CI95%up", 3)
    ## }else{
    ##     samplesQC <- samples
    ## }

    samplesQC$consitentCount <- (samplesQC$Nmeas + samplesQC$broken) == samplesQC$samplecount

    samplesQC.Count <- subset(samplesQC,
                              ! consitentCount,
                              select = c("species", "grade", "stationid",
                                         extraCols,
                                         "samplecount", "broken", "Nmeas"))


    samplesQC.Count$diff <- samplesQC.Count$samplecount -
        (samplesQC.Count$broken + samplesQC.Count$Nmeas)

    samplesQC.Count <- samplesQC.Count[order(abs(samplesQC.Count$diff),
                                             decreasing = TRUE) , ]

    Warn <- FALSE

    if (getOption("hasWeight")){
        ## ####################################################################################################
        ## Reports:

        ## head(samplesQC)

        ## summary(samplesQC$sWg.GLM.corr / samplesQC$sWg.simu.corr)

        ## summary(100 * (samplesQC$sampleweight - samplesQC$sWg.GLM.corr) / samplesQC$sWg.GLM.corr)

        ## idx <- which(abs((samplesQC$sampleweight - samplesQC$sWg.GLM.corr) / samplesQC$sWg.GLM.corr) > 0.3)

        ## test <- samplesQC[idx, ]
        ## test$consistentWg <- test$sampleweight >= test$"WL_CI95%lower" &
        ##     test$sampleweight <= test$"WL_CI95%upper"

        meanBias <- mean(100 * (samplesQC$sampleweight - samplesQC$sWg.GLM.corr) / samplesQC$sWg.GLM.corr,
                         na.rm = TRUE)

        medianBias <- median(100 * (samplesQC$sampleweight - samplesQC$sWg.GLM.corr) / samplesQC$sWg.GLM.corr,
                             na.rm = TRUE)

        w.meanBias <- weighted.mean(x = 100 * (samplesQC$sampleweight -
                                               samplesQC$sWg.GLM.corr) / samplesQC$sWg.GLM.corr,
                                    w = samplesQC$sWg.GLM.corr,
                                    na.rm = TRUE)

        message("\n## Comparison of sample weight vs. mean L-W relationship based estimate",
                "\n## [bias = (sample - LW) / LW]:",
                "\n##\tmedian of bias = ", ifelse(meanBias < 0, "-", "+"),
                prettyNum(abs(medianBias), format = "fg", digits = 2), "%",
                "\n##\tweighted mean of bias = ", ifelse(meanBias < 0, "-", "+"),
                prettyNum(abs(w.meanBias), format = "fg", digits = 2), "%\n")



        Warn <- FALSE
        if (nrow(samplesQC.LW))
        {
            message("\nWarning message:",
                    "\n Some sample weight do not match the weight expected from the length-weight relationships:\n")

            print(samplesQC.LW)
            Warn <- TRUE
        }
    }

    if (nrow(samplesQC.Count))
    {
        message("\nWarning message:",
                "\n Some sample counts do not match the number measured + broken:\n")

        print(samplesQC.Count)
        Warn <- TRUE
    }

    if (Warn) warning("You should check the data!")

}else{
    message("\n## Not length-weight data to run the catch weight QC!")
}


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
