## April Goebl
## Script started 2023-08-02
## BLM Restoration project at Denver Botanic Gardens
## Perform data checks on Chatfield Common Garden raw data




rm(list=ls())



## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------


## SET WORKING DIRECTORY --------------------------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
#PEVI <- read.csv(file="Chatfield/20230303_ChatfieldData_PEVI.csv", sep=",", header=TRUE, dec=".")



## DATA CHECKS
#If was dead previously is it alive now? This may mean late greening. But since wasn't collecting data on this per se, could just change to always alive?
#For some species (like PEVI and maybe BOGR), may be dead, then alive, then dead, etc. For other species (like ERNA and prob ARFR) should not fluctuate like this.
#If flowering previously, should always be 'flowering' after initial flowering
#Look at notes


