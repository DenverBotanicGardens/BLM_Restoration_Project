## April Goebl
## Script started 2023-12-14 (modified from 20231122_ChatfieldCGdataAnalysis_BOGR)
## BLM Restoration project at Denver Botanic Gardens
## Analyze ARFR data from Chatfield Common Garden  


rm(list=ls())
#dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
#library(Hmisc)
library(dplyr)
library(stringr)
#library(tidyr)
library(lme4)
library(plotrix)
library(EnvStats)
#library(car)
library(effects)
library(reshape2)
library(gplots)
library(tidyverse)
library(AICcmodavg)
#library(corrplot)
library(PerformanceAnalytics)
library(vcfR)
#library(SNPRelate)
calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------





## SET WORKING DIRECTORY --------------------------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
ARFR22 <- read.csv(file="Chatfield/2022_data/20230616_ChatfieldData2022_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
ARFR23 <- read.csv(file="Chatfield/2023_data/20240301_ChatfieldData2023_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
ARFR24 <- read.csv(file="Chatfield/2024_data/20241219_ChatfieldData2024_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
ARFR24.surv <- read.csv(file="Chatfield/2024_data/CommonGarden/ARFR/20250321_ChatfieldCGsurveyData2024_ARFR_amgUpdates.csv", sep=",", header=TRUE, dec=".", na.strings="")

ARFR.SdZn <- read.csv(file="AGoebl/Seeds/20231212_ARFR_LatLong_hexcodes.csv", sep=",", header=TRUE, dec=".")
ARFR.biovar <- readRDS("AGoebl/Seeds/20230814_ARFR_BiovarsAvg1980_2021")

## ** load ARFR 2024 field data, from survey123 **

## Load ARFR filtered vcf table 
vcf_filt <- read.vcfR(file="DNA_Seq/AlysonEmery/AGartemisia_v5.recode.vcf")
## ----------------------------------------------------------------------------------------------




## ARFR - DATA CLEAN UP ---------------------------------------------
str(ARFR22)
str(ARFR23)
str(ARFR24)

ARFR22$Source <- as.factor(ARFR22$Source)
ARFR23$Source <- as.factor(ARFR23$Source)
ARFR24$Source <- as.factor(ARFR24$Source)

##If OrigPltSurv_20220527 = 0 & plt not replaced (N), ignore data for this plt, i.e. future surv should be NA, not 0
ARFR22.cl <- ARFR22[ARFR22$OrigPltSurvival_20220527==1 | (ARFR22$OrigPltSurvival_20220527==0 & ARFR22$Replaced_YorN_20220531=="Y"),]

#Note: Don't use OrigPltSurv_20220527 data in days to mort or other field analyses,
#this surv may not correspond to plt names in Source/Pop col (may correspond to orig planted or assigned)

#Some plts that weren't dead & had sz measure on 5/27 were replaced. In this case, don't use length_0527.
#If length_0527 is not NA and Replaced = Y, then ignore length_527 as it doesn't correspond to plt names in source col
ARFR22.cl$Length_cm_20220527[!is.na(ARFR22.cl$Length_cm_20220527) & ARFR22.cl$Replaced_YorN_20220531=="Y"] <- NA

#If not replaced, but died before planting, don't use subsequent surv data
#ARFR22.cl[!is.na(ARFR22.cl$DateMortalityObservedPreTransplant),] #All plts that died before planting were replaced

## Add other checks/ edits to address replacements if needed 


## Checks 
#Were some plts dead that were selected to be harvested? Current datasheet only has "Harvest" marked for plts there were alive
ARFR22.coll <- ARFR22.cl[!is.na(ARFR22.cl$Harvest_20221014) | !is.na(ARFR22.cl$Harvest_20221110),]
ARFR22.coll[is.na(ARFR22.coll$AGB_MinusBag),] #Two harvested plts do not have final AGB: 1107 (alive?) and 1274 (dead?)
ARFR22.MissinfBM <- ARFR22.coll[is.na(ARFR22.coll$InfBM_Wobag_g),]
ARFR22.MissinfBM$ID[ARFR22.MissinfBM$Phenology_20220922==3] #All except 1107 included in 2024 datasheet

#Checks to 2024 re-weighing of 2022 repro BM data
colnames(ARFR24)
ARFR.coll <- ARFR24[!is.na(ARFR24$Harvested_AGB),]
ARFR.infBM <- ARFR24[!is.na(ARFR24$InfBM2022_Wobag_g),]
ARFR.infBMnewH <- ARFR24[!is.na(ARFR24$InfBM2022smpls_HEADS_2024weigh),]
#Check samples w info in 2022 inf weight column that lack data in 2024 chaff/head columns
ARFR.infCheck <- ARFR24[!is.na(ARFR24$InfBM2022_Wobag_g) & is.na(ARFR24$InfBM2022smpls_HEADS_2024weigh),]
#Most have small weights in 2022 column, so maybe tissue not saved. However, ID 524 and 908 are larger. Should these be re-weighed?
# Check all these IDs for flowering: 524, 885, 908, 1853, 1918, 1950
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==524] #Yes, keep this inf BM data
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==885] #Yes, keep this inf BM data
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==908] #Yes, keep this inf BM data
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==1853] #No, don't include inf BM data
ARFR22.cl$Phenology_20220909[ARFR22.cl$ID==1853] 
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==1918] #No, don't include inf BM data
ARFR22.cl$Phenology_20220909[ARFR22.cl$ID==1918] 
ARFR22.cl$Phenology_20220922[ARFR22.cl$ID==1950] #No, don't include inf BM data
ARFR22.cl$Phenology_20220909[ARFR22.cl$ID==1950] 

length(ARFR24$InfBM2022_Wobag_g[!is.na(ARFR24$InfBM2022_Wobag_g)])
ARFR24$InfBM2022smpls_HEADS_2024weigh[!is.na(ARFR24$InfBM2022smpls_HEADS_2024weigh)]
ARFR24$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24$InfBM2022smpls_CHAFF_2024weigh)]
length(ARFR24$InfBM2022smpls_HEADS_2024weigh[!is.na(ARFR24$InfBM2022smpls_HEADS_2024weigh)])
length(ARFR24$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24$InfBM2022smpls_CHAFF_2024weigh)]) #~10 samples had chaff combined w flwr heads for weighing
hist(ARFR24$InfBM2022smpls_HEADS_2024weigh)
hist(ARFR24$InfBM2022smpls_CHAFF_2024weigh)
str(ARFR24$InfBM2022smpls_HEADS_2024weigh)
str(ARFR24$InfBM2022smpls_CHAFF_2024weigh)


## Check that surv is only 1, 0 and maybe NA
ARFR22.cl[ARFR22.cl$Survival_20220622 < 0 | ARFR22.cl$Survival_20220622 > 1,]
ARFR23[(ARFR23$Survival_20230927 < 0 | ARFR23$Survival_20230927 > 1) & !is.na(ARFR23$Survival_20230927),] #Don't use 2023 surv data w/o 2024 data to corroborate 
ARFR23[(ARFR23$Survival_20230615 < 0 | ARFR23$Survival_20230615 > 1) & !is.na(ARFR23$Survival_20230615),]
## ** Add checks to 2024 surv **


## Check that pheno, surv are only integers
ARFR22.cl[ARFR22.cl$Survival_20220622 - floor(ARFR22.cl$Survival_20220622) != 0,]
ARFR22.cl[ARFR22.cl$Survival_20220922 - floor(ARFR22.cl$Survival_20220922) != 0,]
ARFR22.cl[ARFR22.cl$Phenology_20220715 - floor(ARFR22.cl$Phenology_20220715) != 0 & !is.na(ARFR22.cl$Phenology_20220715),]
ARFR23[ARFR23$Survival_20230615 - floor(ARFR23$Survival_20230615) != 0,]
ARFR23[ARFR23$Survival_20230927 - floor(ARFR23$Survival_20230927) != 0 & !is.na(ARFR23$Survival_20230927),]


## Check that length and other relevant traits are only numeric
ARFR22.cl[!is.numeric(ARFR22.cl$Length_cm_20220527),]
ARFR22.cl[!is.numeric(ARFR22.cl$Length_cm_20220622),]
ARFR22.cl[!is.numeric(ARFR22.cl$Length_cm_20220726),]
ARFR23[!is.numeric(ARFR23$Height_20230927),]
ARFR24[!is.numeric(ARFR24$SLA_mm2permg),]
ARFR24[!is.numeric(ARFR24$InfBM2022smples_HEADS_2024weigh) & !is.na(ARFR24$InfBM2022smples_HEADS_2024weigh),]
ARFR24[!is.numeric(ARFR24$InfBM2022smpls_CHAFF_2024weigh),]


## Check that if surv=0 for a given date, there are no height values for that date
ARFR22.cl$Length_cm_20220527[ARFR22.cl$OrigPltSurvival_20220527==0 & !is.na(ARFR22.cl$Length_cm_20220527)]
ARFR22.cl$Length_cm_20220622[ARFR22.cl$Survival_20220622==0 & !is.na(ARFR22.cl$Length_cm_20220622)]
ARFR22.cl$Length_cm_20220726[ARFR22.cl$Survival_20220726==0 & !is.na(ARFR22.cl$Length_cm_20220726)]
ARFR22.cl$Length_cm_20220726[ARFR22.cl$Survival_20220726==0 & !is.na(ARFR22.cl$Length_cm_20220726)]
ARFR23$Height_20230927[ARFR23$Survival_20230927==0 & !is.na(ARFR23$Height_20230927)]


#If plt alive and >0 was not entered in height col, enter NA (not 0)
ARFR22.cl$Length_cm_20220527[ARFR22.cl$OrigPltSurvival_20220527==1 & !is.na(ARFR22.cl$Length_cm_20220527) & ARFR22.cl$Length_cm_20220527==0]
ARFR22.cl$Length_cm_20220622[ARFR22.cl$Survival_20220622==1 & !is.na(ARFR22.cl$Length_cm_20220622) & ARFR22.cl$Length_cm_20220622==0]
ARFR22.cl$Length_cm_20220622[ARFR22.cl$Survival_20220622==1 & !is.na(ARFR22.cl$Length_cm_20220622) & ARFR22.cl$Length_cm_20220622==0] <- NA
ARFR22.cl[ARFR22.cl$Survival_20220622==1 & !is.na(ARFR22.cl$Length_cm_20220622) & ARFR22.cl$Length_cm_20220622==0,]
ARFR22.cl$Length_cm_20220726[ARFR22.cl$Survival_20220726==1 & !is.na(ARFR22.cl$Length_cm_20220726) & ARFR22.cl$Length_cm_20220726==0]
ARFR23$Height_20230927[ARFR23$Survival_20230927==1 & !is.na(ARFR23$Height_20230927) & ARFR23$Height_20230927==0]


#Check that once surv=0, future surv stays zero (if becomes 1 later, could be data entry error)
#To start, could look at all rows with inconsistent survival data. Then either mark as Remove or correct error if obvious
## CONSOLIDATE SURVIVAL DATA
ARFR22.surv <- ARFR22.cl %>% dplyr::select(c(starts_with("Survival_")))
ARFR23.surv <- ARFR23 %>% dplyr::select(c(starts_with("Survival_")))

#Check for survival data inconsistencies
ARFR22.surv <- ARFR22.surv %>% mutate(Check1 = ifelse(Survival_20220715 - Survival_20220622 <= 0, "", "Remove?"),
                                      Check2 = ifelse(Survival_20220721 - Survival_20220715 <= 0, "", "Remove?"),
                                      Check3 = ifelse(Survival_20220726 - Survival_20220721 <= 0, "", "Remove?"),
                                      Check4 = ifelse(Survival_20220804 - Survival_20220726 <= 0, "", "Remove?"),
                                      Check5 = ifelse(Survival_20220817 - Survival_20220804 <= 0, "", "Remove?"),
                                      Check6 = ifelse(Survival_20220830 - Survival_20220817 <= 0, "", "Remove?"),
                                      Check7 = ifelse(Survival_20220909 - Survival_20220830 <= 0, "", "Remove?"),
                                      Check8 = ifelse(Survival_20220922 - Survival_20220909 <= 0, "", "Remove?"))
ARFR22.surv[ARFR22.surv$Check1=="Remove?",]
## ** Combine check columns into 1 and add that to main datasheet **

ARFR23.surv <- ARFR23.surv %>% mutate(Check1 = ifelse(Survival_20230620 - Survival_20230615 <= 0, "", "Remove?"),
                                      Check2 = ifelse(Survival_20230629 - Survival_20230620 <= 0, "", "Remove?"),
                                      Check3 = ifelse(Survival_20230705 - Survival_20230629 <= 0, "", "Remove?"),
                                      Check4 = ifelse(Survival_20230728 - Survival_20230705 <= 0, "", "Remove?"),
                                      Check5 = ifelse(Survival_20230801 - Survival_20230728 <= 0, "", "Remove?"),
                                      Check6 = ifelse(Survival_20230927 - Survival_20230801 <= 0, "", "Remove?"))
ARFR23.surv[ARFR23.surv$Check1=="Remove?",]
ARFR23.surv[ARFR23.surv$Check2=="Remove?",]
ARFR23.surv[ARFR23.surv$Check3=="Remove?",]
ARFR23.surv[ARFR23.surv$Check4=="Remove?",]
ARFR23.surv[ARFR23.surv$Check5=="Remove?",]
ARFR23.surv[ARFR23.surv$Check6=="Remove?",]
## Address NAs in surv data from 9/27, could remove this survey, or just treat NAs appropriately **
## ----------------------------------------------------------------------------------------------






## ARFR - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ARFR.SdZn$Source <- str_replace(ARFR.SdZn$SOURCE_CODE, "4-SOS", "")
ARFR.biovar$Source <- str_replace(ARFR.biovar$Pop, "4-SOS", "")

## Add Population name abbreviation column. Within state should be ordered by increasing lat
ARFR.SdZn$PopAbbrev[grepl("ARFR-AZ930-423-NAVAJO-18", ARFR.SdZn$Source)] = "AZ.1" #"A.AZ.1"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-AZ930-422-NAVAJO-18", ARFR.SdZn$Source)] = "AZ.2" #"A.AZ.2"    
ARFR.SdZn$PopAbbrev[grepl("ARFR-NM930N-66-11", ARFR.SdZn$Source)] = "NM.1" #"A.NM.1"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY930-44-LASANIMAS-13", ARFR.SdZn$Source)] = "CO.1" #"A.CO.1"       
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-314-JEFFERSON-12", ARFR.SdZn$Source)] = "CO.2" #"A.CO.2"  
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-316-JEFFERSON-12", ARFR.SdZn$Source)] = "CO.3" #"A.CO.3"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-UT080-109-UINTAH-12", ARFR.SdZn$Source)] = "UT.1" #"A.UT.1"    
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-294-11", ARFR.SdZn$Source)] = "CO.4" #"A.CO.4"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY040-71-10", ARFR.SdZn$Source)] = "WY.1" #"A.WY.1"       
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY050-151-FREMONT-16", ARFR.SdZn$Source)] = "WY.2" #"A.WY.2" 
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY050-49-FREMONT-12", ARFR.SdZn$Source)] = "WY.3" #"A.WY.3"  

## Edit column names for biovariables
colnames(ARFR.biovar) <- c("Pop","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                           "bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                           "bio17","bio18","bio19","Source")

## Add colour columns that corresponds to pop or seed zone
#SdZn.list <- unique(ARFR.SdZn$SdZone) #If coloring by seed zone, make distinct from lat/ unique pop colors**
ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "powderblue"        #semi-humid, cold
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen"     #semi-humid, cool
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen4"    #semi-humid, warm
ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "darkgoldenrod1"   #semi-arid, cold
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "orange3"         #semi-arid, cool
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "tomato2"         #semi-arid, warm

## Add seed zone name abbreviation column
#ARFR.SdZn$SdZnAbbrev[grepl("5 - 10 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "semi-humid, cold"
#ARFR.SdZn$SdZnAbbrev[grepl("10 - 15 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "semi-humid, cool"
#ARFR.SdZn$SdZnAbbrev[grepl("15 - 20 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "semi-humid, warm"
#ARFR.SdZn$SdZnAbbrev[grepl("5 - 10 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "semi-arid, cold"
#ARFR.SdZn$SdZnAbbrev[grepl("10 - 15 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "semi-arid, cool"
#ARFR.SdZn$SdZnAbbrev[grepl("15 - 20 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "semi-arid, warm"



## Color by latitude (based on Alyson's map and SNP PCA colors)
pop.list <- unique(ARFR.SdZn$Source)
ARFR.SdZn$PopCol[grepl("ARFR-AZ930-423-NAVAJO-18", ARFR.SdZn$Source)] = "#9E0142"        
ARFR.SdZn$PopCol[grepl("ARFR-AZ930-422-NAVAJO-18", ARFR.SdZn$Source)] = "#D53E4F"   
ARFR.SdZn$PopCol[grepl("ARFR-NM930N-66-11", ARFR.SdZn$Source)] = "#F46D43"    
ARFR.SdZn$PopCol[grepl("ARFR-WY930-44-LASANIMAS-13", ARFR.SdZn$Source)] = "#FDAE61"   
ARFR.SdZn$PopCol[grepl("ARFR-CO932-314-JEFFERSON-12", ARFR.SdZn$Source)] = "#FEE08B"        
ARFR.SdZn$PopCol[grepl("ARFR-CO932-316-JEFFERSON-12", ARFR.SdZn$Source)] = "#FFFF85"        
ARFR.SdZn$PopCol[grepl("ARFR-UT080-109-UINTAH-12", ARFR.SdZn$Source)] = "#aabe7f"        
ARFR.SdZn$PopCol[grepl("ARFR-CO932-294-11", ARFR.SdZn$Source)] = "#ABDDA4"     
ARFR.SdZn$PopCol[grepl("ARFR-WY040-71-10", ARFR.SdZn$Source)] = "#66C2A5"    
ARFR.SdZn$PopCol[grepl("ARFR-WY050-151-FREMONT-16", ARFR.SdZn$Source)] = "#3288BD"   
ARFR.SdZn$PopCol[grepl("ARFR-WY050-49-FREMONT-12", ARFR.SdZn$Source)] = "#5E4FA2"

## Add seed zone name abbreviation column (look at BOGR script)..?

## Add column with seed zone or pop 'order' (look at BOGR script)..?
ARFR.SdZn$PopOrder[grepl("ARFR-AZ930-423-NAVAJO-18", ARFR.SdZn$Source)] = "A"        
ARFR.SdZn$PopOrder[grepl("ARFR-AZ930-422-NAVAJO-18", ARFR.SdZn$Source)] = "B"   
ARFR.SdZn$PopOrder[grepl("ARFR-NM930N-66-11", ARFR.SdZn$Source)] = "C"    
ARFR.SdZn$PopOrder[grepl("ARFR-WY930-44-LASANIMAS-13", ARFR.SdZn$Source)] = "D"   
ARFR.SdZn$PopOrder[grepl("ARFR-CO932-314-JEFFERSON-12", ARFR.SdZn$Source)] = "E"        
ARFR.SdZn$PopOrder[grepl("ARFR-CO932-316-JEFFERSON-12", ARFR.SdZn$Source)] = "F"        
ARFR.SdZn$PopOrder[grepl("ARFR-UT080-109-UINTAH-12", ARFR.SdZn$Source)] = "G"        
ARFR.SdZn$PopOrder[grepl("ARFR-CO932-294-11", ARFR.SdZn$Source)] = "H"     
ARFR.SdZn$PopOrder[grepl("ARFR-WY040-71-10", ARFR.SdZn$Source)] = "I"    
ARFR.SdZn$PopOrder[grepl("ARFR-WY050-151-FREMONT-16", ARFR.SdZn$Source)] = "J"   
ARFR.SdZn$PopOrder[grepl("ARFR-WY050-49-FREMONT-12", ARFR.SdZn$Source)] = "K"
## ----------------------------------------------------------------------------------------------




## ARFR - COMBINE DATA TYPES --------------------------------------------
ARFR.biovar <- left_join(ARFR.biovar, ARFR.SdZn, by="Source")
ARFR22.cl <- left_join(ARFR22.cl, ARFR.biovar, by="Source")
ARFR23 <- left_join(ARFR23, ARFR.biovar, by="Source")
ARFR24 <- left_join(ARFR24, ARFR.biovar, by="Source")
## ----------------------------------------------------------------------





## ARFR - ADD GROWTH RATE VARIABLES ------------------------------
## 'Late' growth (June to July)
#ARFR22.cl$GrwthRate_Specific <- log(ARFR22.cl$Length_cm_20220726/ARFR22.cl$Length_cm_20220622)
#ARFR22.cl$GrwthRate_Absolute <- ARFR22.cl$Length_cm_20220726-ARFR22.cl$Length_cm_20220622
ARFR22.cl$GrwthRate_Relative <- (ARFR22.cl$Length_cm_20220726-ARFR22.cl$Length_cm_20220622)/ARFR22.cl$Length_cm_20220622

## Could look at early vs late growth if can use pre-replacement early height measurement 20220527
## ---------------------------------------------------------------



## COMBINE AND ADD REPRO BM DATA ---------------------------------
#ARFR22.cl$InfBM_Wobag_g
#ARFR22.cl$InfBM_Wobag_g[!is.na(ARFR22.cl$InfBM_Wobag_g)] #Number of plts with repro data
#length(ARFR22.cl$InfBM_Wobag_g[!is.na(ARFR22.cl$InfBM_Wobag_g)])
#length(ARFR22.cl$InfBM_Wbag[!is.na(ARFR22.cl$InfBM_Wbag)])
#hist(ARFR22.cl$InfBM_Wobag_g)
#str(ARFR22.cl$InfBM_Wobag_g)

## Combine flwr head and chaff/seed weights + any missed smpls from initial 2022 weights
ARFR24$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24$InfBM2022smpls_HEADS_2024weigh) & is.na(ARFR24$InfBM2022smpls_CHAFF_2024weigh)] <- 0
length(ARFR24$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24$InfBM2022smpls_CHAFF_2024weigh)]) 
ARFR24$InfBM2022_2024updated <- ARFR24$InfBM2022smpls_HEADS_2024weigh + ARFR24$InfBM2022smpls_CHAFF_2024weigh
length(ARFR24$InfBM2022_2024updated[!is.na(ARFR24$InfBM2022_2024updated)])
## ** add several indivs from 2022 col **
## ---------------------------------------------------------------



## COULD ADD AGB VARIABLES AS WELL? 
## For now, I think height is better since more data and correlated to AGB
length(ARFR22.cl$AGB_MinusBag[!is.na(ARFR22.cl$AGB_MinusBag)])
length(ARFR24$AGB2022_MinusBag[!is.na(ARFR24$AGB2022_MinusBag)])



## CLEAN 2023 PLT SZ FIELD MEASUREMENTS
nrow(ARFR23[!is.na(ARFR23$Height_20230927>0),])
ARFR23$Height_20230927[ARFR23$ExcludeSzDueToUncertainty=="Y" & !is.na(ARFR23$ExcludeSzDueToUncertainty)] <- NA
nrow(ARFR23[!is.na(ARFR23$Height_20230927>0),])


## ** When looking at Surv, used exclude Suv col in 2023 to clean entries **
## Use tablet survey data from 2024 **


## COMBINE RELEVANT 2022, 2023, 2024 DATA
# ** could add other lengths and early growth from 2022 later; make sure replacements corrected for **
ARFR22.sel <- ARFR22.cl %>% dplyr::select(c("ID", "Length_cm_20220726", "GrwthRate_Specific", "GrwthRate_Absolute", "GrwthRate_Relative"))               
ARFR23.sel <- ARFR23 %>% dplyr::select(c("ID","Height_20230927")) 
ARFR.cl <- left_join(ARFR24, ARFR23.sel, by="ID") 
ARFR.cl <- left_join(ARFR.cl, ARFR22.sel, by="ID") 
ARFR.cl$Source <- as.factor(ARFR.cl$Source)
## ---------------------------------------------------------------




## Do something similar for ARFR if needed ----------------------------------------------------------
## MAKE FILE WITH JUST COLUMNS/ PHENOTYPES OF INTEREST FOR ANALYSIS WITH GENOTYPE DATA --------------
#ERNA.clSel <- ERNA.cl %>% dplyr::select(c("Source","ID","SentForSequencing_20230803","Length_cm_20220608",
#                                          "Length_cm_20220719","Length_cm_20220915","DaysToFlwr","SLA_mm2permg",
#                                          "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
#                                          "bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19",
#                                          "Longitude","Latitude","HexCode","seed_zone","SdZnAbbrev","PopAbbrev")) 

#write.csv(ERNA.clSel, "Chatfield/20241018_ChatfieldPhenotypes_ERNA.csv", row.names=FALSE)
## --------------------------------------------------------------------------------------------------







## ARFR - TEST FOR TREATMENT EFFECT -------------------------------------------------------
## Review and update models as needed **
#hist(log(ARFR.cl$AGB_MinusBag))
#hist(ARFR$AGB_MinusBag)
#ARFR.tx.mod <- aov(ARFR$AGB_MinusBag ~ ARFR$Source + ARFR$Treatment)
#ARFR.tx.mod <- lmer(log(AGB_MinusBag) ~ Source + Treatment + (1|Block), data=ARFR)
#summary(ARFR.tx.mod)
#ARFR.pop.mod <- lmer(log(AGB_MinusBag) ~ Source + (1|Block), data=ARFR)

#models <- list(ARFR.tx.mod, ARFR.pop.mod)
#mod.names <- c('IncldTx', 'JustPop')
#aictab(cand.set = models, modnames = mod.names )
# No support for treatment 

#hist(log(ARFR$Length_cm_20220726))
#hist(ARFR$Length_cm_20220726)
#ARFR.tx.mod <- aov(log(ARFR$Length_cm_20220726) ~ ARFR$Source + ARFR$Treatment)
#ARFR.tx.modlog <- lmer(log(ARFR$Length_cm_20220726) ~ Source + Treatment + (1|Block), data=ARFR)
#ARFR.pop.modlog <- lmer(log(ARFR$Length_cm_20220726) ~ Source + (1|Block), data=ARFR)
#ARFR.tx.mod <- lmer(ARFR$Length_cm_20220726 ~ Source + Treatment + (1|Block), data=ARFR)
#ARFR.pop.mod <- lmer(ARFR$Length_cm_20220726 ~ Source + (1|Block), data=ARFR)

#models <- list(ARFR.tx.mod, ARFR.pop.mod)
#mod.names <- c('IncldTx', 'JustPop')
#aictab(cand.set = models, modnames = mod.names )
# No support for treatment? 
## -----------------------------------------------------------------------------------------




## ARFR - VISUALIZE RAW DATA ---------------------------------------------------------------

## Order populations for plotting 
## Order by average size or lat or other trait(s)
ARFR.htByMed <- with(ARFR.cl, reorder(Source, Length_cm_20220726, median, na.rm=TRUE))
ARFR.ht23ByMed <- with(ARFR.cl, reorder(Source, Height_20230927, median, na.rm=TRUE))
ARFR.infByMed <- with(ARFR.cl, reorder(Source, InfBM2022smpls_2024reweigh, median, na.rm=TRUE))
ARFR.latByMed <- with(ARFR.cl, reorder(Source, Lat, median, na.rm=TRUE))

ARFR.meds <- ARFR.cl %>% group_by(Source) %>% 
             dplyr::summarise(Height22_MD=median(Length_cm_20220726,na.rm=TRUE), AGB22_MD=median(AGB2022_MinusBag,na.rm=TRUE),
             ReproBMrw_MD=median(InfBM2022smpls_2024reweigh,na.rm=TRUE), Height23_MD=median(Height_20230927,na.rm=TRUE),
             GrowthRe_MD=median(GrwthRate_Relative,na.rm=TRUE), ReproBM22_MD=median(InfBM2022_Wobag_g,na.rm=TRUE),
             GrowthSp_MD=median(GrwthRate_Specific,na.rm=TRUE), GrowthAb_MD=median(GrwthRate_Absolute,na.rm=TRUE))
ARFR.meds <- left_join(ARFR.meds, ARFR.SdZn, by="Source")


## Boxplots of raw data 
ARFR.meds <- ARFR.meds[order(ARFR.meds$Latitude),] #Order by lat

par(mfrow=c(2,3))

## Size 2022
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Height22_MD),] #Order by median sz
#boxplot(Length_cm_20220726 ~ ARFR.htByMed, data=ARFR.cl,
#        xlab=NA, ylab="Height (cm)", cex.lab=1.25,
#        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2,
#        main="FINAL SIZE", cex.main=1.5, col=ARFR.meds$PopCol)
boxplot(Length_cm_20220726 ~ ARFR.latByMed, data=ARFR.cl,
        xlab="Height (cm)", ylab=NA, cex.lab=1.25, horizontal=TRUE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2,
        main="FINAL SIZE 2022", cex.main=1.5, col=ARFR.meds$PopCol)


## Growth rate(s) ** Add time interval to growth rate calcs? **
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Growth_MD),]
boxplot(GrwthRate_Relative ~ ARFR.latByMed, data=ARFR.cl, las=2, horizontal=TRUE,
        xlab="Plant relative growth", ylab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, ylim=c(-0.5,5.5),
        cex.main=1.5, col=ARFR.meds$PopCol, ylim=c(-0.4,7), main="GROWTH RATE 2022")

#boxplot(GrwthRate_Absolute ~ ARFR.htByMed, data=ARFR.cl, las=2,
#        xlab=NA, ylab="Plant absolute growth", cex.lab=1, cex.axis=0.9, names=ARFR.meds$PopAbbrev,
#        cex.main=1.5, col=ARFR.meds$PopCol)

#boxplot(GrwthRate_Specific ~ ARFR.htByMed, data=ARFR.cl, las=2,
#        xlab=NA, ylab="Plant specific growth", cex.lab=1, cex.axis=0.9, names=ARFR.meds$PopAbbrev,
#        cex.main=1.5, col=ARFR.meds$PopCol)
#ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ ARFR.grsByMed, data=ARFR.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=ARFR.SdZn$PopCol)

## Blank plot
#plot.new()
#legend("center", unique(ARFR.meds$Source[order(ARFR.meds$PopOrder, decreasing=TRUE)]), 
#       col=unique(ARFR.meds$PopCol[order(ARFR.meds$PopOrder, decreasing=TRUE)]), cex=1.1, pch=19)


## Repro
#boxplot(InfBM_Wobag_g ~ ARFR.htByMed, data=ARFR.cl, las=2,
#        xlab=NA, ylab="Reproductive biomass", cex.lab=1.25, cex.axis=0.99, names=ARFR.meds$PopAbbrev,
#        cex.main=1.5, col=ARFR.meds$PopCol, main="REPRODUCTIVE OUTPUT")

boxplot(InfBM2022smpls_2024reweigh ~ ARFR.latByMed, data=ARFR.cl, las=2,
        xlab="Reproductive biomass (g)", ylab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=TRUE, ylim=c(0,60),
        cex.main=1.5, col=ARFR.meds$PopCol, main="REPRODUCTIVE OUTPUT 2022")

#ARFR.meds <- ARFR.meds[order(ARFR.meds$ReproBMrw_MD),]
#boxplot(InfBM2022smpls_2024reweigh ~ ARFR.infByMed, data=ARFR.cl, las=2,
#        xlab=NA, ylab="Reproductive biomass", cex.lab=1.25, names=ARFR.meds$PopAbbrev,
#        cex.axis=0.79, main="Artemisia frigida", cex.main=1.5, col=ARFR.meds$PopCol)

#boxplot(InfBM_Wobag_g ~ ARFR.infByMed, data=ARFR.cl, las=2,
#        xlab=NA, ylab="Reproductive biomass", cex.lab=1.25, names=ARFR.meds$PopAbbrev,
#        cex.axis=0.79, main="Artemisia frigida", cex.main=1.5, col=ARFR.meds$PopCol)


## size 2023
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Height23_MD),] #Order by median 2023 sz
boxplot(Height_20230927 ~ ARFR.latByMed, data=ARFR.cl,
        xlab="Height (cm)", ylab=NA, cex.lab=1.25, horizontal=TRUE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2, ylim=c(15,90),
        main="FINAL SIZE 2023", cex.main=1.5, col=ARFR.meds$PopCol)


## SLA 2024
boxplot(SLA_mm2permg ~ ARFR.latByMed, data=ARFR.cl, las=2,
        xlab="Specific leaf area (mm2/mg)", ylab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=TRUE, ylim=c(0,40),
        cex.main=1.5, col=ARFR.meds$PopCol, main="SPECIFIC LEAF AREA 2024")

plot.new()
legend("center", unique(ARFR.meds$Source[order(ARFR.meds$PopOrder, decreasing=TRUE)]), 
       col=unique(ARFR.meds$PopCol[order(ARFR.meds$PopOrder, decreasing=TRUE)]), cex=1.5, pch=19)
## ---------------------------------------------------





## Trait PCA ----------------------------------------
#ARFR.traits <- ARFR.cl %>% dplyr::select(c("Length_cm_20220726","Survival_20220922","Phenology_20220922",
#                                           "AGB_MinusBag","InfBM_Wbag"))
#covMat.traits <- cov(ARFR.traits, use="pairwise.complete.obs")
#pca.results <- prcomp(covMat.traits)
#cols <- viridis(5)
#plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, col=cols, cex=1.2)
#legend("topleft", colnames(ARFR.traits), col=cols, cex=0.75, pch=19)
## ---------------------------------------------------







### VCR table and PCA  --------------------------------------------------------------------------
## Get list of sample names from filtered vcf table 
genotype_mxFilt <- vcfR::extract.gt(vcf_filt, as.numeric=TRUE)
indvNames <- as.data.frame(colnames(genotype_mxFilt))
colnames(indvNames) <- "Sample"

#make column with ID using string replace or something
indvNames$Temp <- str_replace(indvNames$Sample, "ARFR_", "")
indvNames$ID <- as.integer(str_replace(indvNames$Temp, "_sorted", ""))
#join by ID to get source (pop ID)
indvNames <- left_join(indvNames, ARFR24, by="ID")



## TRY PRCOMP PCA FUNCTION 
## Extract, then convert genotypes to allele counts (function suggested by chatGPT)
#gts <- extract.gt(vcf_filt)
#convert_gt_to_allele_counts <- function(gt_matrix) {
#  apply(gt_matrix, 2, function(col) {
#    sapply(col, function(gt) {
#      #if (is.na(gt) || gt == "." || gt == "./.") return(NA)
#      alleles <- unlist(strsplit(gt, "[|/]"))
#      sum(as.numeric(alleles))
#    })
#  })
#}

#allele_counts <- convert_gt_to_allele_counts(gts)

## Calculate covariance matrix (from old PCA GSD code)
#gts <- t(gts) #Note: gts is 012 mx
#covMat<- cov(gts, use="pairwise.complete.obs")
#covMat<- cov(allele_counts, use="pairwise.complete.obs")

## Run PCA
#pca.results = prcomp(covMat) 
#plot(x=pca.results$x[,1],y=pca.results$x[,2],pch=19) #, col=indvNames$ColPop)
#plot(x=pca.results$x[,2],y=pca.results$x[,3],pch=19, col=indvNames$ColPop)
#plot(x=pca.results$x[,3],y=pca.results$x[,4],pch=19, col=indvNames$ColPop)

## Doesn't look right .... 

## Scree plot
#pca.results$sdev[1]**2/sum(pca.results$sdev**2)
#barplot(pca.results$sdev[1:10]**2/sum(pca.results$sdev**2))





## TRY PCADAPT
library(pcadapt)
path_to_vcf <- "DNA_Seq/AlysonEmery/AGartemisia_v5.recode.vcf"

vcf_tbl <- read.pcadapt(path_to_vcf, type="vcf")                            #Read in data from vcf table
x <- pcadapt(vcf_tbl, K=5)                                                  #Run PCAdapt 

plot(x, option="screeplot")
plot(x, option="scores", pop=indvNames$Source)



## Extract PCA scores  
## All samples
xScores <- x$scores    
dfScores <- as.data.frame(xScores[,1:4])                    #Isolate scores for all or select PCs
colnames(dfScores) <- c("PC1score","PC2score","PC3score","PC4score")
dfScores$Sample <- indvNames$Sample               
dfScores$Source <- indvNames$Source                          #Add a column with pop ID
popNames <- unique(indvNames$PopID)


## Calculate mean PC1 values for each population
PC1.mean <- dfScores %>% group_by(Source) %>% summarise(PC1mean = mean(PC1score), n=n())

## ** Look into what 'NA' is can clean up ****
PC1.mean <- PC1.mean[1:11,]

## Create color gradient and assign colors based on numeric continuous PC1 mean values

colorAccording2(
  x,
  gradTy = "logGray",
  nStartOmit = NULL,
  nEndOmit = "sep",
  revCol = FALSE,
  alpha = 1,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)

plot(1:11,PC1.mean$PC1mean,pch=16,cex=2,col=colorAccording2(PC1.mean$PC1mean))
plot(PC1.mean$PC1mean, col=colorAccording2(PC1.mean$PC1mean), pch=19, cex=1.5)

PC1.mean$HexCode <- colorAccording2(PC1.mean$PC1mean)
col <- c("#00FF2EFF", "#00FFB9FF", "#FF008BFF", "#5D00FFFF", "#002EFFFF", "#FF0000FF", "#00B9FFFF")
# Create the color ramp function (from AI)
#color_function <- colorRampPalette(c("blue", "red"))
# Generate colors for the data
#PC1.mean$color <- color_function(PC1.mean$PC1mean)

## --------------------------------------------------------------------------------------------
## Plots source pops on a map and color by seed zone or other characteristics (e.g. PCA scores). 









