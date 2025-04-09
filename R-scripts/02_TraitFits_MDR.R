## Lilian Chan, University of Guelph
## Arctic Vector-borne disease transmission suitability model
##
## Purpose: use Bayesian inference (JAGS) to fit TPCs for mosquito development 
## rate (MDR) for Aedes nigripes
##     1) with uniform priors; and 
##     2) with informative priors using TPCs from Aedes sierrensis data
## 
## Table of content:
##    0. Set-up workspace
##    1. Load and process data
##    2. Fitting TPC
##        A. Briere function (Truncated normally-distributed)
##        B. Quadratic function (Truncated normally-distributed)
##    3. Plotting


##########
###### 0. Set-up workspace ----
##########

library(tidyverse)
library(readxl)
library(janitor)
library(R2jags)
library(mcmcplots) # Diagnostic plots for fits
library(rTPC)
library(nls.multstart)
library(broom)

setwd("~/Documents/UofG/Arctic-VBD")

##########
###### 1. Load and process data ----
##########

nigripes <- 


