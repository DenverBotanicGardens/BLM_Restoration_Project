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
#.. in CGdataAnalysis script?? **

#Excel formula for surv cols: =IF(C2="",1,0) **Note from CGdataAnalysis script **

## For ERNA, if flowering previously, shoudl always be 'flowering'. If comment says 'seed' than phenology should be 1.
## For BOGR height, "NA" or 0 from excel should be NA in R. 
## If dead, phenology should be NA
