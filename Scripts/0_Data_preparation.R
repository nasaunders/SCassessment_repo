#-*- coding: latin-1 -*-

### File: 0_Data_preparation.R
### Time-stamp: <2019-02-04 15:38:19 yreecht>
###
### Created: 06/09/2017	08:43:03
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## ##################################################
## Packages:

mandatoryPack <- c("utils", "coda")

source(file = "./0_Load_packages.R")

source("./lib/Graphical_functions/1_Functions_graphics.R")
source("./lib/Base_functions/0_progressBar.R")

## ##################################################
## Load data:

## #############
## Measurements:
measurements <- read.csv(file = measurementsFile, header = TRUE)

colnames(measurements) <- tolower(gsub("_", ".", colnames(measurements)))
colnames(measurements)[colnames(measurements) == "track.id"] <- "stationid"

measurements$date <- as.Date(x = as.character(measurements$date),
                             format = "%d/%m/%Y")

## Cleaning empty rows (csv export issue):
idxEmptyMeas <- apply(measurements, 1, function(x)all(is.na(x) | x == ""))

## Needs to drop unused levels (empty string "" in particular):
measurements <- droplevels(subset(measurements,
                                  ! idxEmptyMeas))

## head(measurements)
## tail(measurements)

## ########
## Samples:
samples <- read.csv(file = samplesFile, header = TRUE)

colnames(samples) <- tolower(gsub("_", ".", colnames(samples)))
colnames(samples)[colnames(samples) == "track.id"] <- "stationid"

samples$date <- as.Date(x = as.character(samples$date),
                        format = "%d/%m/%Y")

## Cleaning empty rows (csv export issue):
idxEmptySamp <- apply(samples, 1, function(x)all(is.na(x) | x == ""))

## Needs to drop unused levels (empty string "" in particular):
samples <- droplevels(subset(samples,
                             ! idxEmptySamp))

## Seems sensible to assume no bulk splitting by default (exceptional enough):
samples$fraction.raising[is.na(samples$fraction.raising)] <- 1

## head(samples)
## tail(samples)

## Checking data types:
if(any(idxFact <-
           sapply(samples, is.factor)[c("broken",
                                        "totalweight", "totalcount",
                                        "sampleweight", "samplecount", "fraction.raising")]))
{
    warning("Some unexpected factor(s): ",
            paste(names(idxFact)[idxFact], collapse = ", "),
            "\n converted to numeric but you should investigate the cause.")

    for(fac in names(idxFact)[idxFact])
    {
        samples[ , fac] <- as.numeric(as.character(samples[ , fac]))
    }
}

## head(samples)


## ##################################################
## Matching stations:

## Checking consistency of vessels/teams in both data sets:
stopifnot((all(is.na(samples$vessel.team)) && all(is.na(measurements$vessel.team))) ||
          all(levels(measurements$vessel.team) %in% levels(samples$vessel.team)))


## De-ambiguation of stations if needed:
##   Attribute a different letter to each vessel/team and append it to the station number
if (is.factor(measurements$vessel.team) && nlevels(measurements$vessel.team) > 1)
{
    vessInitial <- toupper(substr(levels(samples$vessel.team), 1, 1))
    if (length(unique(vessInitial)) == nlevels(measurements$vessel.team))
    {
        measurements$stationid <- paste(toupper(substr(as.character(measurements$vessel.team), 1, 1)),
                                        toupper(measurements$stationid), sep = "")

        samples$stationid <- paste(toupper(substr(as.character(samples$vessel.team), 1, 1)),
                                   toupper(samples$stationid), sep = "")

        message("\nAttributed letters by vessel/team (station prefix):\n  * ",
                paste(paste(levels(samples$vessel.team),
                            vessInitial, sep = ": "),
                      collapse = "\n  * "))
    }else{
        measurements$stationid <- paste(LETTERS[as.integer(measurements$vessel.team)],
                                        toupper(measurements$stationid), sep = "")

        samples$stationid <- paste(LETTERS[as.integer(samples$vessel.team)],
                                   toupper(samples$stationid), sep = "")

        message("\nAttributed letters by vessel/team (station prefix):\n  * ",
                paste(paste(levels(samples$vessel.team),
                            LETTERS[seq_len(nlevels(samples$vessel.team))], sep = ": "),
                      collapse = "\n  * "))
    }
}

## Filling grades if ommitted:
if (all(is.na(samples$grade)) && all(is.na(measurements$grade)))
{
    message("\n## Assumed to be all ungraded!")
    samples$grade <- "Ungraded"
    measurements$grade <- "Ungraded"
}

## Filling broken if ommitted:
if (all(is.na(samples$broken)))
{
    message("\n## All broken assumed to be zero!!!")
    samples$broken <- 0
}

## ##################################################
## Check data:
source(file = "./1_Data_preparation_checks.R")


## ##################################################
## Calculation of catch data:

## Test if weight recorded at all:
options(hasWeight = any(! is.na(samples[ , "totalweight"])))

catches <- by(data = samples,
              INDICES = as.list(samples[ , c("date", "species", "stationid", "grade")]),
              FUN = function(catchSub, measurements)
           {
               ## print(unique(catchSub$stationid))
               ## if (catchSub$totalweight > catchSub$sampleweight) browser()
               ## if (catchSub$stationid == 112 && catchSub$species == "Ostrea edulis") browser()
               ##
               if((nrow(catchSub) > 1))
               {
                   message("Error:",
                           "\n##\tduplicate (station/species/grade) sample:")

                   print(catchSub)
               }
               stopifnot(nrow(catchSub) == 1)

               ## Subset of the corresponding measurements:
               measSub <- subset(x = measurements,
                                 subset = (date == catchSub$date &
                                           species == as.character(catchSub$species) &
                                           stationid == as.character(catchSub$stationid) &
                                           grade == as.character(catchSub$grade)))

               ## ...and number of individuals measured:
               nMeasured <- nrow(measSub)

               ## In case weight has not been recorded but all individuals (minus broken) have been measured:
               if (all(sapply(catchSub[ , c("totalweight", "totalcount", "sampleweight", "samplecount")],
                              is.na)))
               {
                   catchSub[ , "totalcount"] <-
                       catchSub[ , "samplecount"] <- nMeasured + catchSub[ , "broken"]
               }

               if (is.na(catchSub[ , "samplecount"]) &
                   ! is.na(nMeasured + catchSub[ , "broken"]))
               {
                   catchSub[ , "samplecount"] <- nMeasured + catchSub[ , "broken"]
               }

               ## Total catch in numbers (different ways):
               catchN <- with(catchSub,
                              expr = ifelse(test = ! is.na(totalcount),
                                            yes = totalcount, # total count in priority:
                                            no = ifelse(test = ! is.na(totalweight) & totalweight == 0,
                                                        yes = 0,
                                                        no = ifelse(test = (! is.na(samplecount) &
                                                                            ! is.na(sampleweight) &
                                                                            ! is.na(totalweight)),
                                                                    ## sample count raised by weights if relevant:
                                                                    ##   Include broken already!
                                                                    yes = samplecount * totalweight / sampleweight,
                                                                    ## ... measured+broken used instead otherwise:
                                                                    no = ifelse(test = (! is.na(sampleweight) &
                                                                                        ! is.na(totalweight)),
                                                                                yes = (nMeasured + broken) *
                                                                                    totalweight / sampleweight,
                                                                                ## if everything is NA, use count only:
                                                                                ## (should not be relevant anymore
                                                                                ##  since tested and filled earlier)
                                                                                no = nMeasured + broken)))) *
                                  fraction.raising) # additional raising if bulk splitted.

               ## Total catch in weight (different ways):
               catchW <- with(data = catchSub,
                              expr = ifelse(test = ! is.na(totalweight),
                                            yes = totalweight, # measured total catch in priority.
                                            no = ifelse(test = (! is.na(totalcount) & totalcount == 0 &
                                                                getOption("hasWeight")), # Only if weight recorded.
                                                        yes = 0,
                                                        no = ifelse(test = (! is.na(sampleweight) &
                                                                            ! is.na(samplecount) &
                                                                            ! is.na(totalcount)),
                                                                    ## sample weight raised by counts if relevant:
                                                                    yes = sampleweight * totalcount / samplecount,
                                                                    ## Not calculable otherwise:
                                                                    no = NA))) *
                                  fraction.raising) # additional raising if bulk splitted.

               ## Calculate the overall raising factor (priority to the count, omitted if NA):
               raisingFact <- as.vector(na.omit(c(catchN, catchW) /
                                                unlist(catchSub[ , c("samplecount", "sampleweight")])))[1]

               ## if (any(is.na(raisingFact))) browser()
               data.frame("LD_EventStartDate" = format(unique(catchSub$date),
                                                       format = "%d/%m/%Y"),
                          "LD_HaulNo" = unique(catchSub$stationid),
                          "SpeciesID" = unique(catchSub$species),
                          "CatchUnitID" = c("Number", "KGs"),
                          "Caught" = c(catchN, catchW),
                          "IsTargetY.N" = unique(catchSub$targetspecies),
                          "LD_GearID" = catchSub$gear,
                          "grade" = unique(catchSub$grade),
                          "sampleCountConsistent" = as.logical(catchSub$samplecount == nMeasured + catchSub$broken),
                          "sampleCountDiff" = (catchSub$samplecount -  (nMeasured + catchSub$broken)),
                          "raisingFactor" = raisingFact)
           },
           measurements = measurements, simplify=FALSE)

catches <- do.call(what = rbind, args = catches)


## tail(head(catches, 40), 10)
## summary(catches$Caught)
## table(catches$SpeciesID, catches$sampleCountConsistent, useNA = "always")
## catches[!is.na(catches$sampleCountConsistent), ]

catchesF <- aggregate(data = catches,
                      Caught ~ LD_EventStartDate + LD_HaulNo + SpeciesID + CatchUnitID + IsTargetY.N + LD_GearID,
                      FUN = sum)


## ## head(catchesF)
## catches[is.na(catches$Caught), ]
## catchesF[is.na(catchesF$Caught), ]
## catchesF[catchesF$LD_HaulNo == 613, ]


biological <- by(data = measurements,
                 INDICES = as.list(measurements[ , c("date", "species", "stationid", "grade")]),
                 FUN = function(x, catchesGrades)
              {
                  ## cat(".", unique(x$stationid), ".")

                  catchesGradesSub <- subset(catchesGrades,
                                             LD_EventStartDate == format(x$date, format = "%d/%m/%Y")[1] &
                                             SpeciesID == as.character(x$species[1]) &
                                             LD_HaulNo == as.character(x$stationid[1]) &
                                             grade == as.character(x$grade[1]))

                  ## At least original data are kept:
                  resMes <- x[ , c("length", "weight", "observation")]

                  if (is.na(abs(catchesGradesSub$raisingFactor[1] - 1) > 1e-2))
                  {
                      message("## The raising factor cannot be calculated for ",
                              "\n## \t* Sp:", catchesGradesSub$SpeciesID[1],
                              ", St:", catchesGradesSub$LD_HaulNo[1],
                              ",  Gd:", catchesGradesSub$grade[1])
                      stop("Raising factor issue")
                  }
                  ## Need for adding "data" by size for the fraction over 1 if raising factor > 1
                  if (abs(catchesGradesSub$raisingFactor[1] - 1) > 1e-2)
                  {
                      ## Addition of "dummy" individuals can be carried out if raising factor > 1:
                      if (catchesGradesSub$raisingFactor[1] > 1)
                      {
                          ## browser()
                          if ("Number" %in% catchesGradesSub$CatchUnitID)
                          {
                              ## Underestimation of numbers due to broken individuals tends
                              ##   to propagate proportionally to the raising factor...
                              ##   Best to correct it if possible:
                              raisingFactor <- unlist(subset(catchesGradesSub,
                                                             CatchUnitID == "Number",
                                                             select = "Caught")) / nrow(x)
                          }else{        # Should almost never happen.
                              raisingFactor <- catchesGradesSub$raisingFactor[1]
                          }

                          ## Number of repetitions by size class:
                          reps <- round(tapply(X = x$length,
                                               INDEX = round(x$length),
                                               FUN = length) * (raisingFactor - 1))

                          ## Mean values to repeat:
                          vals <- round(tapply(X = x$length,
                                               INDEX = round(x$length),
                                               FUN = mean), 1)

                          ## Repeated values:
                          resrep <- sapply(seq_along(reps),
                                           function(i) rep(vals[i],
                                                           reps[i]),
                                           simplify = FALSE)

                          resrep <- do.call(c, resrep)

                          ## Appending "dummy" individuals:
                          if (length(resrep)) # rounding may lead to actually no "dummy" indiv.
                          {
                              resMes <- rbind(resMes,
                                              data.frame(length = resrep,
                                                         weight = NA_real_,
                                                         observation = "Dummy indiv."))
                          }
                      }else{
                          stop("negative raising factor for station ", x$stationid[1],
                               ", species \"", x$species[1],
                               "\" and grade \"", x$grade[1], "\"")
                      }
                  }

                  ##
                  return(data.frame("LD_EventStartDate" = format(unique(x$date),
                                                                 format = "%d/%m/%Y"),
                                    "LD_HaulNo" = unique(x$stationid),
                                    "LD_SpeciesID" = unique(x$species),
                                    "Size" = resMes$length,
                                    "Weight" = resMes$weight,
                                    "IsTargetY.N" = unique(catchesGradesSub$IsTargetY.N),
                                    "LD_GearID" = unique(x$gear),
                                    "grade" = unique(x$grade),
                                    "observation" = resMes$observation))
              }, catchesGrades = catches, simplify=FALSE)

biological <- do.call(rbind, biological)

## head(biological, 5)
## dim(biological)
## dim(measurements)
## tapply(catchesF$Caught, catchesF$CatchUnitID, sum, na.rm = TRUE)

## nrow(biological) + sum(samples$broken)

## ##################################################
## Save results:

## Output directory check/create:
if ( ! dir.exists(outputDir))
{
    dir.create(outputDir, recursive = TRUE)
}else{}

if (file.access(outputDir, mode = 2) < 0)
{
    stop("You don't have the rights for writting in\n\t",
         normalizePath(outputDir, winslash = "/"))
}

write.csv(catchesF, file = file.path(outputDir,
                                     paste0(getOption("FilePrefix"),
                                            "Catches",
                                            getOption("FileSuffix"),
                                            ".csv")),
          row.names = FALSE, na = "")

write.csv(biological, file = file.path(outputDir,
                                       paste0(getOption("FilePrefix"),
                                              "Biological",
                                              getOption("FileSuffix"),
                                              ".csv")),
          row.names = FALSE, na = "")


## Quality Checks:
source(file = "./2_Data_preparation_QC.R")

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
