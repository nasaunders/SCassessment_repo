#-*- coding: latin-1 -*-

### File: 0_Load_packages.R
### Time-stamp: <2019-02-04 15:29:42 yreecht>
###
### Created: 14/03/2016	14:22:48
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
### Mandatory packages for analysing survey data
####################################################################################################


library(tcltk)

if ( ! exists("mandatoryPack"))
{
    mandatoryPack <- c("alphahull", "rgeos", "gstat", "geoR",
                       "mapdata", "maptools", "maps", "sp", "rgdal",
                       "utils", "raster", "OpenStreetMap", "rjags", "R2jags",
                       "RColorBrewer", "MASS")
}

installPack.f <- function(pack)
{
    ## Purpose: Installation of missing packages.
    ## ----------------------------------------------------------------------
    ## Arguments: pack : list of missing packages (character string).
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  3 sept. 2010, 13:43

    ## Packages available online:
    ## (if no repository yet defined, the user will have to choose one).
    availablePack <- available.packages()[ , "Package"]

    ## Warning if any missing package is not available:
    if (any(!is.element(pack, availablePack)))
    {
        tkmessageBox(message=paste("The package(s) : \n\n    * '",
                                   paste(pack[!is.element(pack, availablePack)], collapse="'\n    * '"),
                                   "'\n\nis (are) not available!", sep=""),
                     icon="warning")

        res <- "error"
    }else{
        res <- "ok"
    }

    ## Installation des packages disponibles :
    if (any(is.element(pack, availablePack)))
    {
        install.packages(pack[is.element(pack, availablePack)],
                         dependencies=TRUE, lib = .libPaths()[1])
    }else{}

    return(res)
}

loadPackages.f <- function(requiredPack)
{
    ## Purpose: Load required packages and install missing ones
    ## ----------------------------------------------------------------------
    ## Arguments: requiredPack: required packages (character strings)
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  3 sept. 2010, 12:58


    require(tcltk)

    ## Some versions do not switch automatically to user libraries when
    ## admin rights are not given.
    ## ...Windows only:
    if (.Platform$OS.type == "windows")
    {
        env <- environment(.libPaths)
        assign(".lib.loc", shortPathName(get(".lib.loc", envir=env)), envir=env)
        ##  -> All paths in short versions
    }else{}

    ## Installed packages
    installedPack <- installed.packages()[ , "Package"]

    if (any(!is.element(requiredPack, installedPack)))
    {
        on.exit(tkdestroy(WinInstall))

        Done <- tclVar("0")             # User action status.

        ## Missing packages:
        missingPack <- requiredPack[!is.element(requiredPack, installedPack)]

        ## Interface:
        WinInstall <- tktoplevel()
        tkwm.title(WinInstall, "Missing packages")

        B.Install <- tkbutton(WinInstall, text="   Install   ", command=function(){tclvalue(Done) <- "1"})
        B.Cancel <- tkbutton(WinInstall, text="   Cancel    ", command=function(){tclvalue(Done) <- "2"})

        tkgrid(tklabel(WinInstall, text=" "))
        tkgrid(tklabel(WinInstall,
                       text=paste("Missing package(s) :\n\    * '",
                                  paste(missingPack, collapse="'\n    * '"),
                                  "'\n\n Do you want to install them? \n",
                                  "(active internet connection & administration rights required).", sep=""),
                       justify="left"),
               column=1, columnspan=3, sticky="w")

        tkgrid(tklabel(WinInstall, text=" "))

        tkgrid(B.Install, column=1, row=3, sticky="e")
        tkgrid(tklabel(WinInstall, text="       "), column=2, row=3)
        tkgrid(B.Cancel, column=3, row=3, sticky="w")

        tkgrid(tklabel(WinInstall, text=" "), column=4)

        tkbind(WinInstall, "<Destroy>", function(){tclvalue(Done) <- "2"})

        ## Waiting for user action:
        tkwait.variable(Done)

        if (tclvalue(Done) == "1")
        {
            tkdestroy(WinInstall)

            ## Packages installation:
            res <- installPack.f(pack=missingPack)
        }else{
            res <- "abord"
        }

    }else{
        ## All packages installed... nothing else to do:
        res <- "ok"
    }


    ## Testing for the exit status:
    switch(res,
           ok = invisible(sapply(requiredPack, library, character.only=TRUE)),
           stop(paste("You must install manually the following package(s) :\n\n    * '",
                      paste(requiredPack[!is.element(requiredPack, installed.packages()[ , "Package"])],
                            collapse="\n    * '"),
                      "'",
                      "\n\ninstall.packages(c('",
                      paste(requiredPack[!is.element(requiredPack, installed.packages()[ , "Package"])],
                            collapse="', '"),
                      "'))",
                      sep="")))
}


loadPackages.f(requiredPack = mandatoryPack)


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
