#-*- coding: latin-1 -*-

### File: 2_Check_table_data.R
### Time-stamp: <2018-10-08 16:23:06 yreecht>
###
### Created: 27/10/2016	16:56:48
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## Working on subset of data to avoid modification of the original ones:
catchDataNum <- dropLevels.f(subset(catchData,
                                    CatchUnitID == "Number" &
                                    ! is.na(CatchUnitID),
                                    select = c("LD_HaulNo", "SpeciesID", "LD_GearID", "Caught")))

biolCols <- c("LD_HaulNo", "LD_SpeciesID",
              "LD_GearID", "observation")[c("LD_HaulNo", "LD_SpeciesID",
                                            "LD_GearID", "observation") %in%
                                          colnames(biologicalData)]

biologicalDataSel <- dropLevels.f(subset(biologicalData,
                                         select = biolCols))


catchDataNum$SpeciesID <- factor(as.character(catchDataNum$SpeciesID),
                                 levels = union(levels(catchDataNum$"SpeciesID"),
                                                levels(biologicalDataSel$"LD_SpeciesID")))

biologicalDataSel$LD_SpeciesID <- factor(as.character(biologicalDataSel$LD_SpeciesID),
                                         levels = union(levels(catchDataNum$"SpeciesID"),
                                                        levels(biologicalDataSel$"LD_SpeciesID")))


catchDataNum$LD_GearID <- factor(catchDataNum$LD_GearID,
                                 levels = union(levels(catchDataNum$"LD_GearID"),
                                                levels(biologicalDataSel$"LD_GearID")))

biologicalDataSel$LD_GearID <- factor(biologicalDataSel$LD_GearID,
                                      levels = union(levels(catchDataNum$"LD_GearID"),
                                                     levels(biologicalDataSel$"LD_GearID")))

if (all(is.numeric(as.numeric(c(catchDataNum$LD_HaulNo,
                                biologicalDataSel$LD_HaulNo)))))
{
    haulLevels <- as.character(sort(as.numeric(unique(c(catchDataNum$LD_HaulNo,
                                                        biologicalDataSel$LD_HaulNo)))))
}else{
    haulLevels <- as.character(sort(unique(c(catchDataNum$LD_HaulNo,
                                             biologicalDataSel$LD_HaulNo))))
}

catchDataNum$LD_HaulNo <- factor(catchDataNum$LD_HaulNo,
                                 levels = haulLevels)

biologicalDataSel$LD_HaulNo <- factor(biologicalDataSel$LD_HaulNo,
                                      levels = haulLevels)


## Individuals count, NAs <- 0:
catchCount <- tapply(X = catchDataNum$Caught,
                     INDEX = as.list(catchDataNum[ , c("LD_HaulNo", "SpeciesID", "LD_GearID")]),
                     FUN = sum, na.rm = TRUE)

catchCount[is.na(catchCount)] <- 0

mensurationCount <- tapply(X = biologicalDataSel$LD_HaulNo,
                           INDEX = as.list(biologicalDataSel[ , c("LD_HaulNo", "LD_SpeciesID", "LD_GearID")]),
                           FUN = length)

mensurationCount[is.na(mensurationCount)] <- 0

## Detecting species-station-gear combination with probable raising:
if ( ! is.null(biologicalDataSel$observation))
{
    raised <- tapply(X = biologicalDataSel$observation,
                     INDEX = as.list(biologicalDataSel[ , c("LD_HaulNo", "LD_SpeciesID", "LD_GearID")]),
                     FUN = function(x)
                     {
                         any(grepl(pattern = "^(Fake|Dummy) indiv\\.",
                                   x = x))
                     }) |                   ## [!!!]
        abs(catchCount - round(catchCount)) > sqrt(.Machine$double.eps) # decimal count!

    raised[is.na(raised)] <- FALSE
}else{   ## If no observation field:
    raised <- abs(catchCount - round(catchCount)) > sqrt(.Machine$double.eps) # decimal count!

    raised[is.na(raised)] <- FALSE
}

## catchCount - mensurationCount

resCheck <- sapply(X = species,
                   FUN = function(sp, Catch, Mens, raised)
            {
                message("\n## #########", rep("#", nchar(sp)),
                        "\n## Species ", sp, ":")

                ## Less individual caught than measure AND without raising:
                lessCatch <- (Catch[ , sp, ] < Mens[ , sp, ]) & ! raised[ , sp, ]
                if (any(lessCatch) && isTRUE(speciesCatchUnit[sp, "Number"])) # Only if number measured.
                {
                    stationNames <- if (is.null(dim(lessCatch)))
                                    {
                                        names(lessCatch)[lessCatch]
                                    }else{
                                        paste(rep(row.names(lessCatch), ncol(lessCatch)),
                                              rep(colnames(lessCatch), each = nrow(lessCatch)),
                                              sep = "_")[lessCatch]
                                    }
                    warning("More individuals measured than caught in station (caught, measured):\n",
                            "      ",
                            paste(paste(stationNames, " (",
                                        Catch[ , sp, ][lessCatch], ", ",
                                        Mens[ , sp, ][lessCatch], ")", sep = ""), collapse = ", "))
                }

                noMeasure <- Catch[ , sp, ] > 0 & Mens[ , sp, ] == 0

                if (any(noMeasure) && any(Mens[ , sp, ] > 0))
                {
                    stationNames2 <- if (is.null(dim(noMeasure)))
                                     {
                                         names(noMeasure)[noMeasure]
                                     }else{
                                         paste(rep(row.names(noMeasure), ncol(noMeasure)),
                                               rep(colnames(noMeasure), each = nrow(noMeasure)),
                                               sep = "_")[noMeasure]
                                     }
                    warning("No measurement when individuals caught in station(caught, measured):\n",
                            "      ",
                            paste(paste(stationNames2, "(",
                                        Catch[ , sp, ][noMeasure], ", ",
                                        Mens[ , sp, ][noMeasure], ")", sep = ""), collapse = ", "))
                }

                if ( ! ((any(lessCatch) && isTRUE(speciesCatchUnit[sp, "Number"])) ||
                        (any(noMeasure) && any(Mens[ , sp, ] > 0))))
                {
                    message("No problem found!\n")
                }

                return(c(lessCatch = any(lessCatch) && isTRUE(speciesCatchUnit[sp, "Number"]),
                         noMeasure = any(noMeasure) && any(Mens[ , sp, ] > 0)))
            }, Catch = catchCount, Mens = mensurationCount, raised = raised, simplify = FALSE)

resCheck <- do.call(rbind, resCheck)

if (any(resCheck))
{
    def <- ifelse(any(resCheck[ , 1]), "N", "Y")
    ask <- readline(paste(ifelse(def == "Y",
                                 "Possible inconsistencies", "Inconsistencies"),
                          " found, do you want to keep going? ",
                          c(Y = "(Y)/N", N = "Y/(N)")[def],
                          ": ", sep = ""))
    answ <- ifelse(def == "N",
                   ifelse(is.element(tolower(substr(ask, 1, 1)), c("y", "o")), "Y", "N"),
                   ifelse(is.element(tolower(substr(ask, 1, 1)), c("n")), "N", "Y"))

    if (answ == "N") stop("Data need corrections!")
}


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
