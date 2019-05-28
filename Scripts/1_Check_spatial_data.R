#-*- coding: latin-1 -*-

### File: 1_Check_spatial_data.R
### Time-stamp: <2016-10-25 13:05:26 yreecht>
###
### Created: 25/10/2016	12:51:56
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

if(! (class(stationData) %in% c("SpatialLinesDataFrame", "SpatialPointsDataFrame")))
{
    stop("Station data must be of class \"SpatialLinesDataFrame\"",
         " or  \"SpatialPointsDataFrame\"!")
}

if ( ! all(c("TRACKID", "TIMESTART") %in% colnames(stationData@data)))
{
    stop("Following field missing in the station data:\n  * \"",
         paste(c("TRACKID", "TIMESTART")[ ! c("TRACKID", "TIMESTART") %in%
                                          colnames(stationData@data)],
               collapse = "\"\n  * \""),
         "\"")
}






### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
