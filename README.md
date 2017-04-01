
<!-- README.md is generated from README.Rmd. Please edit that file -->
CRSSIO
======

[![Travis-CI Build Status](https://travis-ci.org/rabutler/CRSSIO.svg?branch=master)](https://travis-ci.org/rabutler/CRSSIO)

R Package to manage code for manipulating the input and output data for CRSS.

Usage
-----

### Creating or modifying input files

-   `createCRSSDNFInputFiles` creates CRSS natural flow input files from posted natural flow data (<http://www.usbr.gov/lc/region/g4000/NaturalFlow/current.html>) or from the `CoRiverNF` package ([CoRiverNF](https://www.github.com/BoulderCodeHub/CoRiverNF))
-   These files can also be created with a GUI through an R Studio Addin (see `?createCRSSInputAddIn`)
-   `changeStartDate` changes the start date of the natural flow input files
-   `trimCCNFFiles` trims the climate change hydrology input files to a specified time period
-   Vectors of the natural flow gage names (`nfGageNames`), along with corresponding CRSS natural inflow input slot names (`CRSSNFInputNames`), corresponding CRSS natural salt input slot names (`CRSSNatSaltInputNames`), and corresponding short, i.e., variable, names (`nfShortNames`)

### Processing CRSS output

-   Functions (`sysCondSALMatrix` and `createSysCondTable`) to create the standard System Conditions Table from CRSS output. Commonly refered to as the "5-year table" but it can go through as many years as simulation data exists. Ex:

``` r
library(CRSSIO)
library(RWDataPlyr) # install_github("BoulderCodeHub/RWDataPlyr")
# create the slot aggregation list
slotAggList <- RWDataPlyr::createSlotAggList(CRSSIO::sysCondSALMatrix())
# use example data in RWDataPlyr to create system condition table
# first get all of the data
scenFolder <- 'DNF,CT,IG'
scenName <- 'DNF Hydrology'
scenPath <- system.file('extdata','Scenario/',package = 'RWDataPlyr')
sysData <- RWDataPlyr::getDataForAllScens(scenFolder, scenName, slotAggList,
                                          scenPath, 'tmp.feather', TRUE)
# then create the system condition table
sysCondTable <- createSysCondTable(sysData, 2017:2021)
# sysCondTable[['limitedTable']] to access results
```

Installation
------------

Only available from GitHub. Use the following to install:

``` r
install.packages('devtools')
library(devtools)
devtools::install_github('BoulderCodeHub/CRSSIO')
```

Log:
----

For details, see the [News](NEWS.md)

-   2017-03-31: version 0.4.0 available
-   2016-10-04: version 0.3 available
-   2016-05-30: version 0.2.1 available
-   2016-05-05: version 0.2 available
-   2015-02-10: version 0.1 available
