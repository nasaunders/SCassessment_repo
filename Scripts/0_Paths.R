#-*- coding: latin-1 -*-

### File: 0_Pathes.R
### Time-stamp: <2018-11-19 11:51:39 yreecht>
###
### Created: 23/03/2016	10:14:05
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Check and create (if needed) the result tree.
####################################################################################################


if ( ! dir.exists(ResultsPath))
{
    dir.create(ResultsPath, recursive = TRUE)
}else{}

if (file.access(ResultsPath, mode = 2) < 0)
{
    stop("You don't have the rights for writting in\n\t",
         normalizePath(ResultsPath, winslash = "/"))
}

if ( ! dir.exists(file.path(ResultsPath, "shapeFiles")))
{
    dir.create(file.path(ResultsPath, "shapeFiles"), recursive = TRUE)
}else{}

if (file.access(file.path(ResultsPath, "shapeFiles"), mode = 2) < 0)
{
    stop("You don't have the rights for writting in the directory\n\t",
         normalizePath(file.path(ResultsPath, "shapeFiles"), winslash = "/"))
}

if ( ! dir.exists(file.path(ResultsPath, "Bayesian/Diagnostics")))
{
    dir.create(file.path(ResultsPath, "Bayesian/Diagnostics"), recursive = TRUE)
}else{}

if (file.access(file.path(ResultsPath, "Bayesian/Diagnostics"), mode = 2) < 0)
{
    stop("You don't have the rights for writting in the directory\n\t",
         normalizePath(file.path(ResultsPath, "Bayesian/Diagnostics"), winslash = "/"))
}

if ( ! dir.exists(file.path(ResultsPath, "Geostatistics/Diagnostics")))
{
    dir.create(file.path(ResultsPath, "Geostatistics/Diagnostics"), recursive = TRUE)
}else{}

if (file.access(file.path(ResultsPath, "Geostatistics/Diagnostics"), mode = 2) < 0)
{
    stop("You don't have the rights for writting in the directory\n\t",
         normalizePath(file.path(ResultsPath, "Geostatistics/Diagnostics"), winslash = "/"))
}

if ( ! dir.exists(file.path(ResultsPath, "Geostatistics/Tiff")))
{
    dir.create(file.path(ResultsPath, "Geostatistics/Tiff"), recursive = TRUE)
}else{}

if (file.access(file.path(ResultsPath, "Geostatistics/Tiff"), mode = 2) < 0)
{
    stop("You don't have the rights for writting in the directory\n\t",
         normalizePath(file.path(ResultsPath, "Geostatistics/Tiff"), winslash = "/"))
}

if ( ! dir.exists(file.path(ResultsPath, "Stratified_(deprecated)")))
{
    dir.create(file.path(ResultsPath, "Stratified_(deprecated)"), recursive = TRUE)
}else{}

if (file.access(file.path(ResultsPath, "Stratified_(deprecated)"), mode = 2) < 0)
{
    stop("You don't have the rights for writting in the directory\n\t",
         normalizePath(file.path(ResultsPath, "Stratified_(deprecated)"), winslash = "/"))
}


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
