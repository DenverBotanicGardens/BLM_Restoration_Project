## April Goebl
## Script started 2024-02-04 (modified from 20231215_ChatfieldCGdataAnalysis_ERNA)
## BLM Restoration project at Denver Botanic Gardens
## Analyze PEVI data from Chatfield Common Garden  


rm(list=ls())
dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(Hmisc)
library(dplyr)
library(stringr)
library(tidyr)
library(lme4)
library(emmeans)
library(plotrix)
library(EnvStats)
library(car)
library(effects)
library(reshape2)
library(gplots)
library(tidyverse)
library(AICcmodavg)
library(corrplot)
library(PerformanceAnalytics)
calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------





## SET WORKING DIRECTORY --------------------------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## ------------------------------------------------------------------------------------------------




## LOAD DATA --------------------------------------------------------------------------------------
#PEVI.SdZn <- read.csv(file="AGoebl/Seeds/20231215_PEVI_LatLongSdZn_hexcodes.csv", sep=",", header=TRUE, dec=".")
PEVI23 <- read.csv(file="Chatfield/20240129_ChatfieldData2023_PEVI.csv", sep=",", header=TRUE, dec=".")
## ----------------------------------------------------------------------------------------------




## PEVI - DATA CLEAN UP ---------------------------------------------
str(PEVI23)
PEVI23$Source <- as.factor(PEVI23$Source)

## ** Look at clean up steps in files for other species to see what applies for PEVI ** 

#If plt alive and 1 was not entered in flowering col, enter 0 (not NA)
PEVI23$Flowering_20230522[PEVI23$Survival_20230522==1 & is.na(PEVI23$Flowering_20230522)] <- 0
PEVI23$Flowering_20230606[PEVI23$Survival_20230606==1 & is.na(PEVI23$Flowering_20230606)] <- 0
PEVI23$Flowering_20230613[PEVI23$Survival_20230613==1 & is.na(PEVI23$Flowering_20230613)] <- 0
PEVI23$Flowering_20230620[PEVI23$Survival_20230620==1 & is.na(PEVI23$Flowering_20230620)] <- 0
PEVI23$Flowering_20230627[PEVI23$Survival_20230627==1 & is.na(PEVI23$Flowering_20230627)] <- 0
PEVI23$Flowering_20230718[PEVI23$Survival_20230718==1 & is.na(PEVI23$Flowering_20230718)] <- 0
PEVI23$Flowering_20230808[PEVI23$Survival_20230808==1 & is.na(PEVI23$Flowering_20230808)] <- 0 
## ---------------------------------------------------------------



## Checks 
## Check that surv is only 1, 0 and maybe NA
PEVI23[(PEVI23$Survival_20230808 < 0 | PEVI23$Survival_20230808 > 1) & !is.na(PEVI23$Survival_20230808),]

#Check that surv is only integers
PEVI23[PEVI23$Survival_20230718 - floor(PEVI23$Survival_20230718) != 0,]

#Flowering cols should only be 1, 0, NA 
min(PEVI23$Flowering_20230606, na.rm=TRUE)
max(PEVI23$Flowering_20230606, na.rm=TRUE)
min(PEVI23$Flowering_20230613, na.rm=TRUE)
max(PEVI23$Flowering_20230613, na.rm=TRUE)
min(PEVI23$Flowering_20230620, na.rm=TRUE)
max(PEVI23$Flowering_20230620, na.rm=TRUE)
min(PEVI23$Flowering_20230627, na.rm=TRUE)
max(PEVI23$Flowering_20230627, na.rm=TRUE)
min(PEVI23$Flowering_20230718, na.rm=TRUE)
max(PEVI23$Flowering_20230718, na.rm=TRUE)
min(PEVI23$Flowering_20230808, na.rm=TRUE)
max(PEVI23$Flowering_20230808, na.rm=TRUE)

# ** Check that phenology only increases or stays the same; once flowering=1, it shouldn't go back to 0 

# ** Check that if surv=0 for a given date, there are no phenology or height values for that date
## ----------------------------------------------------------------------------------------------



## PEVI - DATA MODS -----------------------------------------------------------------------------
# ** obtain seed zones and biovars for PEVI **
## ----------------------------------------------------------------------------------------------




## PEVI - COMBINE DATA TYPES --------------------------------------------------------------------
## ** once have the seed zone and biovar files ** 
#PEVI.biovar <- left_join(PEVI.biovar, PEVI.SdZn, by="Source")
#PEVI23 <- left_join(PEVI23, PEVI.biovar, by="Source")
## ----------------------------------------------------------------------------------------------




## PEVI 2023 - CREATE FLOWERING/ PHENOLOGY VARIABLES ---------------------------------------------
## For 2023, estimate days to 1st flwr & if flowered at all based on '1' in any pheno survey 

## Calculate days to first flower 
PEVI.StartDate <- as.Date("2023-01-01")  #Use first day of the year as arbitrary date 
PEVI.PhenoCol.List <- colnames(PEVI23)[grepl("Flowering*", colnames(PEVI23))]   #Obtain phenology column names
PEVI.Pheno.List <- str_replace(PEVI.PhenoCol.List, "Flowering_", "")            #Obtain just date from phenology columns
PEVI.Pheno.List <- as.Date(PEVI.Pheno.List, "%Y%m%d")
PEVI.DaysToFlwr <- PEVI.Pheno.List - PEVI.StartDate                             #Calculate days from Jan 1 to each phenology survey 

## Loop over each phenology column & enter the num days since Jan 1 when a 1 (bud or later repro stage) first appears
PEVI23$DaysToFlwr <- NA
for (pp in 1:length(PEVI.PhenoCol.List)) {
  PEVI23$DaysToFlwr[PEVI23[,PEVI.PhenoCol.List[pp]]==1 & is.na(PEVI23$DaysToFlwr)] <- as.integer(PEVI.DaysToFlwr)[pp]
}
PEVI23 %>% group_by(Source) %>% dplyr::summarise(Pheno_Avg=mean(DaysToFlwr,na.rm=TRUE), Pheno_SD=sd(DaysToFlwr,na.rm=TRUE))

## Make a flowered Yes or No column
PEVI23$FlwrYesNo <- NA 
PEVI23$FlwrYesNo[!is.na(PEVI23$DaysToFlwr)] <- 1                               #If there's a value in Days to Flwr, then enter Yes
## **If plt alive and didn't flwr, enter No **Find better way, as not sure what date(s) to use for 
PEVI23$FlwrYesNo[is.na(PEVI23$DaysToFlwr) & (PEVI23$Survival_20230808==1 | PEVI23$Survival_20230718==1 |
                                             PEVI23$Survival_20230627==1)] <- 0  
nrow(PEVI23[PEVI23$FlwrYesNo==0,])                                             #Num that survived but didn't flower 
PEVI23 %>% group_by(Source) %>% dplyr::summarise(FlwrYesNo_Avg=mean(FlwrYesNo,na.rm=TRUE)) #Or could try sum and/or NUM=n()
## ---------------------------------------------------------------