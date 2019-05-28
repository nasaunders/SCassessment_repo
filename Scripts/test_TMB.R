#-*- coding: latin-1 -*-

### File: test_TMB.R
### Time-stamp: <2016-12-06 09:56:23 yreecht>
###
### Created: 06/12/2016	09:47:10
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

library(TMB)


compile("../testTMBnormal.cpp")
dyn.load(dynlib("testTMBnormal"))







### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
