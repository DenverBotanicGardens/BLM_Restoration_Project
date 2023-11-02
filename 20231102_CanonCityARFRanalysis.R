## April Goebl
## Script started 2023-11-02
## BLM Restoration project at Denver Botanic Gardens
## Analyze data from Canyon City genetic diversity experiment   




rm(list=ls())
dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
library(lme4)
library(plotrix)
library(EnvStats)
library(car)
library(effects)
library(tidyverse)
library(AICcmodavg)

calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------



## SET WORKING DIRECTORY --------------------------------------------------------------------------
#setwd("C:/Users/april.goebl/OneDrive - Denver Botanic Gardens/AGoebl_NonProject_WorkRelated/BLM_Restoration_Project")
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/CanyonCity")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv(file="20231102_CanonCityData_ARFR.csv", sep=",", header=TRUE, dec=".")
## ------------------------------------------------------------------------------------------------



## LOOK AT AND MODIFY STRUCTURE OF DATA -----------------------------------------------------------
str(dats)
dats$Plot.Number <- as.factor(dats$Plot.Number)
dats$PlotTagNum <- as.factor(dats$PlotTagNum)
dats$Site <- as.factor(dats$Site)
dats$Treatment <- as.factor(dats$Treatment)
dats$Contents <- as.factor(dats$Contents)
## ------------------------------------------------------------------------------------------------



## ----
unique(dats$Contents)
#Remove empty, control, and pilot categories 

