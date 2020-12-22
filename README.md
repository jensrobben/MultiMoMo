# Multi-population Mortality Modelling

# The `MultiMoMo` package <img src="man/figures/flags2.png" alt="" align="right" height="200">

This is the source code for the  `MultiMoMo` package, which is still under development.

## Installation
To install `MultiMoMo` from GitHub, first install the `devtools` package:

``` r
install.packages('devtools')
devtools::install_github('RobbenJ/MultiMoMo')
```

## Overview
The goal of the package `MultiMoMo` is to build a multipopulation mortality model using the Li-Lee approach.
In short, the following concepts are tackled: 
1. Downloading annual mortality data from HMD and Eurostat: period death counts and period exposures-to-risk per gender, per age, per country and per year. 
Period exposures from Eurostat are not available, but are calculated using the HMD-protocol.
2. Calibrating the multi-population mortality model (Li-Lee model) which consists out of two Lee-Carter models: one for the common trend (multi-population trend) and one for the
country-specific deviation from this trend.
3. Fitting a time series process (RWD or AR(k)-process) to the time-dependent parameters in the Li-Lee model.
4. Projecting the time dependent parameters to future years using the fitted time series processes and using multiple simulations.
5. Closing (the simulations of) the mortality rates.
6. Calculating period and cohort life expectancies. 
