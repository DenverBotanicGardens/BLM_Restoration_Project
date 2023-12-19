## April Goebl
## Script started 2023-11-22 (modified from 20230129_ChatfieldCGdataAnalysis)
## BLM Restoration project at Denver Botanic Gardens
## Analyze BOGR data from Chatfield Common Garden  




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
#setwd("C:/Users/april/Denver Botanic Gardens/Conservation - BLM-Grassland")
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
BOGR <- read.csv(file="Chatfield/20230301_ChatfieldData_BOGR.csv", sep=",", header=TRUE, dec=".")

BOGR.SdZn <- read.csv(file="AGoebl/Seeds/20230824_BOGR_LatLongSdZn.csv", sep=",", header=TRUE, dec=".")
BOGR.biovar <- readRDS("AGoebl/Seeds/20230825_BOGR_BiovarsAvg1980_2021")
## ----------------------------------------------------------------------------------------------



## BOGR - DATA CLEAN UP ---------------------------------------------------------------------------
#str(BOGR)
BOGR$Source <- as.factor(BOGR$Source)

## If OrigPltSurv_20220531=0 & plt not replaced (N), ignore data for this plt, i.e. future surv should be NA, not 0
BOGR.cl <- BOGR[BOGR$OrigPltSurvival_20220531==1 | (BOGR$OrigPltSurvival_20220531==0 & BOGR$Replaced_YorN_20220606=="Y"),]

## Note: Don't use OrigPltSurv_20220531 data in days to mort or other field analyses,
#this surv may not correspond to plt names in Source/Pop col (may correspond to orig planted or assigned)

## Check that all numinf_coll_0927 are entered in numinf_0927 - Yes 20230824

## Check that surv is only 1, 0 and maybe NA
#BOGR.cl[BOGR.cl$Survival_20220713 < 0 | BOGR.cl$Survival_20220713 > 1,]

#Check that pheno, surv, numInf are only integers
#BOGR.cl[BOGR.cl$Survival_20220627 - floor(BOGR.cl$Survival_20220627) != 0,]
#BOGR.cl[BOGR.cl$Survival_20220927 - floor(BOGR.cl$Survival_20220927) != 0,]
BOGR.cl[BOGR.cl$Phenology_20220701 - floor(BOGR.cl$Phenology_20220701) != 0 & 
       !is.na(BOGR.cl$Phenology_20220701),]
BOGR.cl[BOGR.cl$NumInflorescence_20220927 - floor(BOGR.cl$NumInflorescence_20220927) != 0 & 
       !is.na(BOGR.cl$NumInflorescence_20220927),]

# ** Check that length is only numeric
# ** Check that once zero in surv on 6/27 or later, stays zero (cannot become 1 later) ** Look at attempt for ARFR **
# ** Check that if surv=0 for a given date, there are no phenology or height values for that date

# ** Trying to figure out for block 3 how measured heights from 5/31 correspond to plts names in Source col
#given replacements made on 6/6. There was an error in data entry (one row off; a non-BOGR measured at the end of a row/ block) around 5/31-6/2
#but correction to the datasheet was not done until around 6/28. 

## ** Is there an analog for BOGR?
#Some plts that weren't dead & had sz measure on 5/27 were replaced. In this case, don't use length_0527.
#If length_0527 is not NA and Replaced = Y, then ignore length_527 as it doesn't correspond to plt names in source col
#ARFR.cl$Length_cm_20220527[!is.na(ARFR.cl$Length_cm_20220527) & ARFR.cl$Replaced_YorN_20220531=="Y"] <- NA
## ----------------------------------------------------------------------------------------------




## BOGR - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
BOGR.SdZn$Source <- str_replace(BOGR.SdZn$Code, "2-SOS", "")
BOGR.biovar$Source <- str_replace(BOGR.biovar$Pop, "2-SOS", "")
BOGR.biovar$Source[16] <- str_replace(BOGR.biovar$Source[16], "BOGRNM-930", "BOGR-NM930")

## Add Population name abbreviation column. W/in state ordered by increasing latitude
BOGR.SdZn$PopAbbrev[grepl("BOGR-AZ040-302-SANTACRUZ-17", BOGR.SdZn$Source)] = "B.AZ.1"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-AZ040-284-SANTACRUZ-17", BOGR.SdZn$Source)] = "B.AZ.2" 
BOGR.SdZn$PopAbbrev[grepl("BOGR-AZ040-278-SANTACRUZ-17", BOGR.SdZn$Source)] = "B.AZ.3"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-AZ930-302-10", BOGR.SdZn$Source)] = "B.AZ.4" 
BOGR.SdZn$PopAbbrev[grepl("BOGR-AZ930-303-10", BOGR.SdZn$Source)] = "B.AZ.5"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-360-FREMONT-17", BOGR.SdZn$Source)] = "B.CO.1"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-362-CHAFFEE-17", BOGR.SdZn$Source)] = "B.CO.2"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-307-JEFFERSON-12", BOGR.SdZn$Source)] = "B.CO.3"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-274-11", BOGR.SdZn$Source)] = "B.CO.4"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-363-CLEARCREEK-17", BOGR.SdZn$Source)] = "B.CO.5"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-CO932-276-11", BOGR.SdZn$Source)] = "B.CO.6" 
BOGR.SdZn$PopAbbrev[grepl("BOGR-NM930-071-08", BOGR.SdZn$Source)] = "B.NM.1"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-NM930-145-10", BOGR.SdZn$Source)] = "B.NM.2"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-NM080-67-CHAVES-16", BOGR.SdZn$Source)] = "B.NM.3"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-NM930-140-10", BOGR.SdZn$Source)] = "B.NM.4"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-NM930N-95-SANDOVAL-12", BOGR.SdZn$Source)] = "B.NM.5"   
BOGR.SdZn$PopAbbrev[grepl("BOGR-UT030-215-GARFIELD-13", BOGR.SdZn$Source)] = "B.UT.1"    
BOGR.SdZn$PopAbbrev[grepl("BOGR-UT080-163-CARBON-14", BOGR.SdZn$Source)] = "B.UT.2"    
BOGR.SdZn$PopAbbrev[grepl("BOGR-WY050-133-FREMONT-16", BOGR.SdZn$Source)] = "B.WY.1"       
BOGR.SdZn$PopAbbrev[grepl("BOGR-WY050-113-FREMONT-15", BOGR.SdZn$Source)] = "B.WY.2" 
BOGR.SdZn$PopAbbrev[grepl("BOGR-WY070-71-JOHNSON-15", BOGR.SdZn$Source)] = "B.WY.3"  

## Edit column names for biovariables
colnames(BOGR.biovar) <- c("Pop","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                           "bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                           "bio17","bio18","bio19","Source")

## Add colour columns that corresponds to seed zone
unique(BOGR.SdZn$seed_zone)
BOGR.SdZn$SdZnColful[grepl("10 - 15 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen2"    #semi-humid, cool #hex#B4EEB4
BOGR.SdZn$SdZnColful[grepl("15 - 20 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen"     #semi-humid, warm
BOGR.SdZn$SdZnColful[grepl("20 - 25 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen4"    #semi-humid, v.warm
BOGR.SdZn$SdZnColful[grepl("10 - 15 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "lightgoldenrod1" #semi-arid, cool
BOGR.SdZn$SdZnColful[grepl("15 - 20 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "goldenrod2"      #semi-arid, warm
BOGR.SdZn$SdZnColful[grepl("20 - 25 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orange3"         #semi-arid, v.warm
BOGR.SdZn$SdZnColful[grepl("25 - 30 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orangered"       #semi-arid, hot
BOGR.SdZn$SdZnColful[grepl("30 - 35 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orangered3"      #semi-arid, v.hot
BOGR.SdZn$SdZnColful[grepl("5 - 10 Deg. F. / 12 - 30", BOGR.SdZn$seed_zone)] = "pink2"           #arid, cold


## 'Order'/color by latitude 
#unique(BOGR.SdZn$Source)
#BOGR.popCols <- colfunc(length(unique(BOGR.SdZn$Source)))
#BOGR.SdZn$PopCol[grepl("BOGR-AZ040-302-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[1]   
#BOGR.SdZn$PopCol[grepl("BOGR-AZ040-284-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[2] 
#BOGR.SdZn$PopCol[grepl("BOGR-AZ040-278-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[3]   
#BOGR.SdZn$PopCol[grepl("BOGR-AZ930-302-10", BOGR.SdZn$Source)] = BOGR.popCols[8] 
#BOGR.SdZn$PopCol[grepl("BOGR-AZ930-303-10", BOGR.SdZn$Source)] = BOGR.popCols[9]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-360-FREMONT-17", BOGR.SdZn$Source)] = BOGR.popCols[12]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-362-CHAFFEE-17", BOGR.SdZn$Source)] = BOGR.popCols[13]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-307-JEFFERSON-12", BOGR.SdZn$Source)] = BOGR.popCols[14]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-274-11", BOGR.SdZn$Source)] = BOGR.popCols[15]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-363-CLEARCREEK-17", BOGR.SdZn$Source)] = BOGR.popCols[16]   
#BOGR.SdZn$PopCol[grepl("BOGR-CO932-276-11", BOGR.SdZn$Source)] = BOGR.popCols[18] 
#BOGR.SdZn$PopCol[grepl("BOGR-NM-930-071-08", BOGR.SdZn$Source)] = BOGR.popCols[4]   
#BOGR.SdZn$PopCol[grepl("BOGR-NM930-145-10", BOGR.SdZn$Source)] = BOGR.popCols[5]   
#BOGR.SdZn$PopCol[grepl("BOGR-NM080-67-CHAVES-16", BOGR.SdZn$Source)] = BOGR.popCols[6]   
#BOGR.SdZn$PopCol[grepl("BOGR-NM930-140-10", BOGR.SdZn$Source)] = BOGR.popCols[7]   
#BOGR.SdZn$PopCol[grepl("BOGR-NM930N-95-SANDOVAL-12", BOGR.SdZn$Source)] = BOGR.popCols[10]   
#BOGR.SdZn$PopCol[grepl("BOGR-UT030-215-GARFIELD-13", BOGR.SdZn$Source)] = BOGR.popCols[11]    
#BOGR.SdZn$PopCol[grepl("BOGR-UT080-163-CARBON-14", BOGR.SdZn$Source)] = BOGR.popCols[17]    
#BOGR.SdZn$PopCol[grepl("BOGR-WY050-133-FREMONT-16", BOGR.SdZn$Source)] = BOGR.popCols[19]       
#BOGR.SdZn$PopCol[grepl("BOGR-WY050-113-FREMONT-15", BOGR.SdZn$Source)] = BOGR.popCols[20] 
#BOGR.SdZn$PopCol[grepl("BOGR-WY070-71-JOHNSON-15", BOGR.SdZn$Source)] = BOGR.popCols[21]  


## Add seed zone name abbreviation column
BOGR.SdZn$SdZnAbbrev[grepl("5 - 10 / arid", BOGR.SdZn$New_label)] = "arid, cold"   
BOGR.SdZn$SdZnAbbrev[grepl("30 - 35 / semi-arid", BOGR.SdZn$New_label)] = "semi-arid, v.hot"   
BOGR.SdZn$SdZnAbbrev[grepl("25 - 30 / semi-arid", BOGR.SdZn$New_label)] = "semi-arid, hot"  
BOGR.SdZn$SdZnAbbrev[grepl("20 - 25 / semi-arid", BOGR.SdZn$New_label)] = "semi-arid, v.warm"   
BOGR.SdZn$SdZnAbbrev[grepl("15 - 20 / semi-arid", BOGR.SdZn$New_label)] = "semi-arid, warm"   
BOGR.SdZn$SdZnAbbrev[grepl("10 - 15 / semi-arid", BOGR.SdZn$New_label)] = "semi-arid, cool"   
BOGR.SdZn$SdZnAbbrev[grepl("20 - 25 / semi-humid", BOGR.SdZn$New_label)] = "semi-humid, v.warm"   
BOGR.SdZn$SdZnAbbrev[grepl("15 - 20 / semi-humid", BOGR.SdZn$New_label)] = "semi-humid, warm"   
BOGR.SdZn$SdZnAbbrev[grepl("10 - 15 / semi-humid", BOGR.SdZn$New_label)] = "semi-humid, cool"   

## Add column with seed zone 'order'
BOGR.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "A"    #semi-humid, cool
BOGR.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "B"    #semi-humid, warm
BOGR.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "C"    #semi-humid, v.warm
BOGR.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "D"   #semi-arid, cool
BOGR.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "E"   #semi-arid, warm
BOGR.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "F"   #semi-arid, v.warm
BOGR.SdZn$SdZnOrder[grepl("25 - 30 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "G"   #semi-arid, hot
BOGR.SdZn$SdZnOrder[grepl("30 - 35 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "H"   #semi-arid, v.hot
BOGR.SdZn$SdZnOrder[grepl("5 - 10 Deg. F. / 12 - 30", BOGR.SdZn$seed_zone)] = "I"   #arid, cold
## ----------------------------------------------------------------------------------------------




## BOGR - COMBINE DATA TYPES --------------------------------------------
BOGR.cl <- left_join(BOGR.cl, BOGR.SdZn, by="Source")
BOGR.cl <- left_join(BOGR.cl, BOGR.biovar, by="Source")
## ---------------------------------------------------------------




## BOGR - CONSOLIDATE PHENOLOGY DATA AND DO CHECKS ------------------------------------
#Pheno cols should only be b/w 2-6 
BOGR.Pheno <- BOGR.cl %>% dplyr::select(starts_with("Phenology"))
min(BOGR.Pheno, na.rm=TRUE)
max(BOGR.Pheno, na.rm=TRUE)

#min(BOGR.cl$Phenology_20220627, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220701, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220705, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220713, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220718, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220725, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220801, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220809, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220818, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220829, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220915, na.rm=TRUE)
#min(BOGR.cl$Phenology_20220927, na.rm=TRUE)
## Mins all should be and are 2

#max(BOGR.cl$Phenology_20220627, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220701, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220705, na.rm=TRUE)
#BOGR.cl[which.max(BOGR.cl$Phenology_20220705),]
#max(BOGR.cl$Phenology_20220713, na.rm=TRUE) 
#BOGR.cl[BOGR.cl$Phenology_20220705>3 & !is.na(BOGR.cl$Phenology_20220705),]
#max(BOGR.cl$Phenology_20220718, na.rm=TRUE)
#BOGR.cl[which.max(BOGR.cl$Phenology_20220718),]
#max(BOGR.cl$Phenology_20220725, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220801, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220809, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220818, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220829, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220915, na.rm=TRUE)
#max(BOGR.cl$Phenology_20220927, na.rm=TRUE)
#BOGR.cl[which.max(BOGR.cl$Phenology_20220927),]
## Typos for max (e.g. 55 instead of 5) corrected in csv file 

# ** Check that phenology only increases or stays the same, or at least that once it's a 3 it doesn't decrease 


## Calculate days to first flower 
BOGR.PlantingDate <- as.Date("2022-05-12")
BOGR.PhenoCol.List <- colnames(BOGR.cl)[grepl("Phenology*", colnames(BOGR.cl))] #Obtain phenology column names
BOGR.Pheno.List <- str_replace(BOGR.PhenoCol.List, "Phenology_", "")            #Obtain just date from phenology columns
BOGR.Pheno.List <- as.Date(BOGR.Pheno.List, "%Y%m%d")
BOGR.DaysToFlwr <- BOGR.Pheno.List - BOGR.PlantingDate                          #Calculate days from planting to each phenology survey 

## Loop over each phenology column & enter the num days since planting when a 3 (flower or later phenophase) first appears
BOGR.cl$DaysToFlwr <- NA
for (pp in 1:length(BOGR.PhenoCol.List)) {
  BOGR.cl$DaysToFlwr[BOGR.cl[,BOGR.PhenoCol.List[pp]]>=3 & is.na(BOGR.cl$DaysToFlwr)] <- as.integer(BOGR.DaysToFlwr)[pp]
}
## ---------------------------------------------------------------




## BOGR - MAKE COMBINED NUM INF VARIABLE --------------------------------
#Make a combined numInf col that combines numinf_927 w numinfCol1014 when the latter are from blocks 7-10 only
BOGR.cl$NumInf <- c(BOGR.cl$NumInflorescence_20220927[BOGR.cl$Block<7], BOGR.cl$NumInf_Coll_20221014[BOGR.cl$Block>=7])
#BOGR.cl %/% unite(BOGR.cl, col='NumInfAll', c('NumInf','NumInc_Coll_20221108'), sep='-')

# ** If numInf col has value greater than 0, phenol for corresponding date should be >2 **
## ---------------------------------------------------------------



## BOGR - ADD GROWTH RATE VARIABLES -------------------------------------
#BOGR.cl$GrwthRate_Specific <- log(BOGR.cl$Length_cm_20220801/BOGR.cl$Length_cm_20220627)
#BOGR.cl$GrwthRate_Absolute <- BOGR.cl$Length_cm_20220801-BOGR.cl$Length_cm_20220627
BOGR.cl$GrwthRate_Relative <- (BOGR.cl$Length_cm_20220801-BOGR.cl$Length_cm_20220627)/BOGR.cl$Length_cm_20220627
## ---------------------------------------------------------------




## BOGR - TEST FOR TREATMENT EFFECT -------------------------------------
hist(log(BOGR.cl$Length_cm_20220801))
hist(BOGR.cl$Length_cm_20220801)
BOGR.tx.mod <- lmer(Length_cm_20220801 ~ Source + Treatment + (1|Block), data=BOGR.cl)
BOGR.pop.mod <- lmer(Length_cm_20220801 ~ Source + (1|Block), data=BOGR.cl)
models <- list(BOGR.tx.mod, BOGR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# Little support for treatment 

hist(log(BOGR.cl$GrwthRate_Specific))
hist(BOGR.cl$GrwthRate_Specific)
BOGR.tx.mod <- lmer(GrwthRate_Specific ~ Source + Treatment + (1|Block), data=BOGR.cl)
BOGR.pop.mod <- lmer(GrwthRate_Specific ~ Source + (1|Block), data=BOGR.cl)
models <- list(BOGR.tx.mod, BOGR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment

hist(log(BOGR.cl$NumInf))
hist(BOGR.cl$NumInf)
BOGR.tx.mod <- glmer(log(NumInf) ~ Source + Treatment + (1|Block), data=BOGR.cl)
BOGR.pop.mod <- glmer(log(NumInf) ~ Source + (1|Block), data=BOGR.cl)
BOGR.tx.modP <- glmer(NumInf ~ Source + Treatment + (1|Block), data=BOGR.cl, family=poisson (link="log"))
BOGR.pop.modP <- glmer(NumInf ~ Source + (1|Block), data=BOGR.cl, family=poisson (link="log"))

models <- list(BOGR.tx.mod, BOGR.pop.mod)
models <- list(BOGR.tx.modP, BOGR.pop.modP)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment
## -----------------------------------------------------------------------------------------






## BOGR - VISUALIZE RAW DATA -----------------------------------------------------------------------------
## Order populations for plotting 
## EITHER --- Order by seed zone -----
## Order based on aridity & temp (cool & wet to hot & dry), w/in each cat, order from N to S


## OR --- Order by average trait value ---
BOGR.htByMed <- with(BOGR.cl, reorder(Source, Length_cm_20220801, median, na.rm=TRUE))
#BOGR.infByMed <- with(BOGR.cl, reorder(Source, NumInf, median, na.rm=TRUE))
#BOGR.graByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Absolute, median, na.rm=TRUE))
#BOGR.grrByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Relative, median, na.rm=TRUE))
#BOGR.grsByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Specific, median, na.rm=TRUE))
#BOGR.dfByMed <- with(BOGR.cl, reorder(Source, DaysToFlwr, median, na.rm=TRUE))

BOGR.meds <- BOGR.cl %>% group_by(Source) %>% 
             dplyr::summarise(Height_MD=median(Length_cm_20220801,na.rm=TRUE),
                              Inf_MD=median(NumInf,na.rm=TRUE),Growth_MD=median(GrwthRate_Relative,na.rm=TRUE),
                              DaysToFlwr_MD=median(DaysToFlwr,na.rm=TRUE))
BOGR.meds <- left_join(BOGR.meds, BOGR.SdZn, by="Source") 

## Boxplots of raw data 
par(mfrow=c(1,1))
## Size
#BOGR.SdZn <- BOGR.SdZn[order(BOGR.SdZn$Lat),] #Order by lat
#BOGR.cl$Source <- factor(BOGR.cl$Source, levels=BOGR.SdZn$Source)
#boxplot(Length_cm_20220801 ~ Source, data=BOGR.cl,
#        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5, col=BOGR.SdZn$PopCol,
#        cex.axis=0.5, names=BOGR.SdZn$PopAbbrev, main="Bouteloua gracilis", cex.main=1.5)

BOGR.meds <- BOGR.meds[order(BOGR.meds$Height_MD),] #Order by median 
boxplot(Length_cm_20220801 ~ BOGR.htByMed, data=BOGR.cl,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.25,
        cex.axis=0.7, names=BOGR.meds$PopAbbrev, las=2,
        main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)
legend("topleft", unique(BOGR.meds$SdZnAbbrev[order(BOGR.meds$SdZnOrder)]), 
       col=unique(BOGR.meds$SdZnColful[order(BOGR.meds$SdZnOrder)]), cex=1.95, pch=19)

## Growth rate(s) ** Add time interval to growth rate calcs? **
#BOGR.meds <- BOGR.meds[order(BOGR.meds$Growth_MD),]
boxplot(GrwthRate_Relative ~ BOGR.htByMed, data=BOGR.cl, las=2,
        xlab="Population", ylab="Plant relative growth", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)

#BOGR.SdZn <- BOGR.SdZn[order(BOGR.SdZn$GR_AB),]
#boxplot(GrwthRate_Absolute ~ BOGR.graByMed, data=BOGR.cl, las=2,
#        xlab="Population", ylab="Plant Absolute growth", cex.lab=1.25, names=BOGR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.SdZn$PopCol)

#BOGR.SdZn <- BOGR.SdZn[order(BOGR.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ BOGR.grsByMed, data=BOGR.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=BOGR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.SdZn$PopCol)

## Repro
#BOGR.meds <- BOGR.meds[order(BOGR.meds$Inf_MD),]
boxplot(NumInf ~ BOGR.htByMed, data=BOGR.cl, las=2,
        xlab="Population", ylab="Number of inflorescences", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)

## Phenology
#BOGR.meds <- BOGR.meds[order(BOGR.meds$DaysToFlwr_MD),]
boxplot(DaysToFlwr ~ BOGR.htByMed, data=BOGR.cl, las=2,
        xlab="Population", ylab="Days to flower", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)
## ---------------------------------------------------





## BOGR - ESTIMATE VARIATION B/W POPULATIONS -------------------------------------------------------------
## Use Emmeans ----------------
## Final size (height)
#BOGR.htMod <- lmer(Length_cm_20220801 ~ Source + (1|Block), data=BOGR.cl)
#BOGR.pair.ht <- emmeans(BOGR.htMod, specs = pairwise ~ Source)
#summary(BOGR.pair.ht)
#plot(BOGR.pair.ht)

#BOGR.ht.pw <- as.data.frame(BOGR.pair.ht$contrasts)
#BOGR.ht.pwDiffs <- BOGR.ht.pw$estimate
#max(BOGR.ht.pwDiffs) 


## Plant growth
#BOGR.grMod <- lmer(GrwthRate_Relative ~ Source + (1|Block), data=BOGR.cl)
#BOGR.pair.gr <- emmeans(BOGR.grMod, specs = pairwise ~ Source)
#summary(BOGR.gr.em)
#plot(BOGR.pair.gr)

#BOGR.gr.means <- as.data.frame(BOGR.pair.gr$emmeans)
#BOGR.gr.meanEsts <- BOGR.gr.means$emmean
#max(BOGR.gr.meanEsts)
#min(BOGR.gr.meanEsts)

#BOGR.gr.pw <- as.data.frame(BOGR.pair.gr$contrasts)
#BOGR.gr.pwDiffs <- exp(BOGR.gr.pw$estimate)
#max(BOGR.gr.pwDiffs)
## ----------------------


## GLMs ------------
## Survival
BOGR.survMod <- glmer(Survival_20220927 ~ Source + (1|Block), data=BOGR.cl, family=binomial (link="logit"))
Anova(BOGR.survMod) #No support for Source in Surv model

## Days until flower

## Number of inflorescences
BOGR.infMod <- glmer(NumInf ~ Source + (1|Block), data= BOGR.cl, family=poisson (link="log"))
Anova(BOGR.infMod) #Support for Source in Inf model
## ----------------------
## ---------------------------------------------------------------------------------------------------





## BOGR - ESTIMATE VARIATION WITHIN POPULATIONS ------------------------------------------------------
## Calculate the coefficient of variation 
BOGR.cv <- BOGR.cl %>% group_by(Source) %>% summarise(Height_CV=cv(Length_cm_20220801, na.rm=TRUE),
                                                         Growth_CV=cv(GrwthRate_Relative, na.rm=TRUE),
                                                         DaysToFlwr_CV=cv(DaysToFlwr, na.rm=TRUE),
                                                         Inf_CV=cv(NumInf, na.rm=TRUE))

BOGR.cv <- left_join(BOGR.cv, BOGR.SdZn, by="Source")
## ---------


## Barplots
## Order populations for plotting BY VALUE
BOGR.cv <- BOGR.cv[order(BOGR.cv$Height_CV),]
barplot(BOGR.cv$Height_CV, xlab="Population", ylab="CV in plant height", cex.lab=1.3, las=2, 
        names=BOGR.cv$PopAbbrev, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.cv$SdZnCol, cex.names=0.6)

BOGR.cv <- BOGR.cv[order(BOGR.cv$Growth_CV),]
barplot(BOGR.cv$Growth_CV, xlab="Population", ylab="CV in plant relative growth", cex.lab=1.3, cex.axis=1, las=2, 
        names=BOGR.cv$PopAbbrev, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.cv$SdZnCol, cex.names=0.6)

BOGR.cv <- BOGR.cv[order(BOGR.cv$DaysToFlwr_CV),]
barplot(BOGR.cv$DaysToFlwr_CV, xlab="Population", ylab="CV in days to first flower", cex.lab=1.3, cex.axis=1, las=2,
        names=BOGR.cv$PopAbbrev, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.cv$SdZnCol, cex.names=0.6)

BOGR.cv <- BOGR.cv[order(BOGR.cv$Inf_CV),]
barplot(BOGR.cv$Inf_CV, xlab="Population", ylab="CV in inflorescence number", cex.lab=1.3, cex.axis=1, las=2,
        names=BOGR.cv$PopAbbrev, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.cv$SdZnCol, cex.names=0.6)


## Make a vertical stripchart showing the CVs for all traits for each population
#For ordering the plot, calculate the difference b/w the min and max CV value for each pop
#BOGR.cvMin <- BOGR.cv %>% dplyr::select(c(ends_with("_CV"))) %>% apply(1, FUN=min)
#BOGR.cvMax <- BOGR.cv %>% dplyr::select(c(ends_with("_CV"))) %>% apply(1, FUN=max)
#BOGR.cv$CVdif <- BOGR.cvMax-BOGR.cvMin
#BOGR.cv <- BOGR.cv[order(BOGR.cv$CVdif),]

BOGR.cvH <- as.data.frame(cbind(BOGR.cv$Source, BOGR.cv$Height_CV))
BOGR.cvH$Trait <- "Height_CV"
BOGR.cvG <- as.data.frame(cbind(BOGR.cv$Source, BOGR.cv$Growth_CV))
BOGR.cvG$Trait <- "Growth_CV"
BOGR.cvR <- as.data.frame(cbind(BOGR.cv$Source, BOGR.cv$Inf_CV))
BOGR.cvR$Trait <- "Inf_CV"
BOGR.cvP <- as.data.frame(cbind(BOGR.cv$Source, BOGR.cv$DaysToFlwr_CV))
BOGR.cvP$Trait <- "DaysToFlwr_CV"
BOGR.cvAll <- rbind(BOGR.cvH, BOGR.cvG, BOGR.cvR, BOGR.cvP)
colnames(BOGR.cvAll) <- c("Source","CV","Trait")

stripchart(as.numeric(CV) ~ Source, data=BOGR.cvAll, vertical=TRUE, pch=19, group.names=BOGR.cv$PopAbbrev, 
           las=2, ylab="Coefficient of variation", cex.lab=1.4, cex=1.25)
BOGR.Gcol <- BOGR.cvAll$Trait == "Growth_CV"
stripchart(as.numeric(CV) ~ Source, data=BOGR.cvAll[BOGR.Gcol,], col="dodgerblue", vertical=TRUE, pch=19,cex=1.25, add=TRUE)
BOGR.Rcol <- BOGR.cvAll$Trait == "Inf_CV"
stripchart(as.numeric(CV) ~ Source, data=BOGR.cvAll[BOGR.Rcol,], col="blueviolet", vertical=TRUE, pch=19,cex=1.25, add=TRUE)
BOGR.Pcol <- BOGR.cvAll$Trait == "DaysToFlwr_CV"
stripchart(as.numeric(CV) ~ Source, data=BOGR.cvAll[BOGR.Pcol,], col="bisque3", vertical=TRUE, pch=19,cex=1.25, add=TRUE)


#test <- BOGR.cv %>% dplyr::select(c(ends_with("_CV")))
BOGR.cvT <- as.data.frame(rbind(as.vector(BOGR.cv$Height_CV), as.vector(BOGR.cv$Growth_CV),
                  as.vector(BOGR.cv$DaysToFlwr_CV), as.vector(BOGR.cv$Inf_CV)))
#BOGR.cvT[5,] <- colSums(BOGR.cvT)
colnames(BOGR.cvT) <- BOGR.cv$Source
barplot(BOGR.cvT, names=BOGR.cv$PopAbbrev, las=2, col=c("black","dodgerblue","bisque3","mediumpurple3"))
## ------------------------------------------------------------------------------




## BOGR - LOOK AT RELATIONSHIP BETWEEN TRAIT VAR AND AVERAGE ---------------------------
BOGR.ht.mn <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Height_MN=mean(Length_cm_20220801, na.rm=TRUE))
BOGR.gr.mn <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Growth_MN=mean(GrwthRate_Relative, na.rm=TRUE))
BOGR.df.mn <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(DaysToFlwr_MN=mean(DaysToFlwr, na.rm=TRUE))
BOGR.inf.mn <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Inf_MN=mean(NumInf, na.rm=TRUE))

BOGR.ht.sd <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Height_SD=sd(Length_cm_20220801, na.rm=TRUE))
BOGR.gr.sd <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Growth_SD=sd(GrwthRate_Relative, na.rm=TRUE))
BOGR.df.sd <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(DaysToFlwr_SD=sd(DaysToFlwr, na.rm=TRUE))
BOGR.inf.sd <- BOGR.cl %>% group_by(Source) %>% dplyr::summarise(Inf_SD=sd(NumInf, na.rm=TRUE))

plot(BOGR.ht.sd$Height_SD, BOGR.ht.mn$Height_MN, pch=19, col=BOGR.SdZn$SdZnColful, cex=1.3)
plot(BOGR.gr.sd$Growth_SD, BOGR.gr.mn$Growth_MN, pch=19, col=BOGR.SdZn$SdZnColful, cex=1.2)
plot(BOGR.df.sd$DaysToFlwr_SD, BOGR.df.mn$DaysToFlwr_MN, pch=19, col=BOGR.SdZn$SdZnColful, cex=1.2)
plot(BOGR.inf.sd$Inf_SD, BOGR.inf.mn$Inf_MN, pch=19,col=BOGR.SdZn$SdZnColful, cex=1.1)
## ----------------------------------------------------------


## BOGR - LOOK AT CORRELATION BETWEEN TRAITS ----------------------------------------------
par(mfrow=c(3,2))
plot(BOGR.ht.mn$Height_MN, BOGR.gr.mn$Growth_MN, pch=19, col="black", cex=1.3,
     xlab="Mean plant height", ylab="Mean relative growth rate", cex.lab=1.2)
plot(BOGR.ht.mn$Height_MN, BOGR.df.mn$DaysToFlwr_MN, pch=19, col="black", cex=1.3,
     xlab="Mean plant height", ylab="Mean days to first flower", cex.lab=1.2)
plot(BOGR.ht.mn$Height_MN, BOGR.inf.mn$Inf_MN, pch=19, col="black", cex=1.3,
     xlab="Mean plant height", ylab="Mean number of inflorescences", cex.lab=1.2)
plot(BOGR.df.mn$DaysToFlwr_MN, BOGR.gr.mn$Growth_MN, pch=19, col="black", cex=1.3,
     xlab="Mean days to first flower", ylab="Mean relative growth rate", cex.lab=1.2)
plot(BOGR.inf.mn$Inf_MN, BOGR.gr.mn$Growth_MN, pch=19, col="black", cex=1.3,
     xlab="Mean number of inflorescences", ylab="Mean relative growth rate", cex.lab=1.2)
plot(BOGR.inf.mn$Inf_MN, BOGR.df.mn$DaysToFlwr_MN, pch=19, col="black", cex=1.3,
     xlab="Mean number of inflorescences", ylab="Mean days to first flower", cex.lab=1.2)
## ----------------------------------------------------------------------------------------




## BOGR - EVALUATE RELATIONSHIPS B/W TRAITS AND SOURCE CLIMATE -------------------------------------------
## 19 Bio-climate variables 
## Models
#BOGR.ht.mod <- lmer(Length_cm_20220801 ~ Latitude.x + Longitude.x + (1|Block), data=BOGR.cl)
#Anova(BOGR.ht.mod)
#plot(allEffects(BOGR.ht.mod))
#plot(predictorEffects(BOGR.ht.mod))

#BOGR.ht.bv.mod <- lmer(Length_cm_20220801 ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
#                       + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
#                       + scale(bio17) + (1|Block), data=BOGR.cl)
#Anova(BOGR.ht.bv.mod)
#plot(allEffects(BOGR.ht.bv.mod))
#plot(predictorEffects(BOGR.ht.bv.mod))

#BOGR.inf.bv.mod <- glmer(NumInf ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
#                       + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
#                       + scale(bio17) + (1|Block), data=BOGR.cl, family=poisson(link="log"))
#Anova(BOGR.inf.bv.mod)
#plot(allEffects(BOGR.inf.bv.mod))
#plot(predictorEffects(BOGR.inf.bv.mod))

## ** Pick different model form? **
#BOGR.df.bv.mod <- lmer(log(DaysToFlwr) ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
#                        + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
#                        + scale(bio17) + (1|Block), data=BOGR.cl)
#Anova(BOGR.df.bv.mod)
#plot(allEffects(BOGR.df.bv.mod))
#plot(predictorEffects(BOGR.df.bv.mod))
## --------------------------------------------------------------


## Look at PCA of 19 bioclim variables to reduce number of predictors -------------------
BOGR.biovar <- left_join(BOGR.biovar, BOGR.SdZn, by="Source")

#BOGR.biovarCen <- scale(BOGR.biovar[,2:20], center=TRUE, scale=TRUE)
#BOGR.biovarComb <- as.data.frame(cbind(BOGR.biovar[,21],BOGR.biovarCen))
BOGR.pcaBiovar <- prcomp(BOGR.biovar[,2:20], scale=TRUE)
par(pty="s")
par(mfrow=c(1,1))
plot(x=BOGR.pcaBiovar$x[,1], y=BOGR.pcaBiovar$x[,2], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,2], y=BOGR.pcaBiovar$x[,3], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,3], y=BOGR.pcaBiovar$x[,4], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)

## Look at scree plot
BOGR.pcaBiovar$sdev[1]**2/sum(BOGR.pcaBiovar$sdev**2)
BOGR.pcaBV.varExpl <- BOGR.pcaBiovar$sdev^2/sum(BOGR.pcaBiovar$sdev^2)
barplot(BOGR.pcaBiovar$sdev[1:19]**2/sum(BOGR.pcaBiovar$sdev**2))
sum(BOGR.pcaBV.varExpl[1:3]) #top 6 PC axes explain over 98% of variation

## Add arrows on PCA plot and look at loadings 
biplot(BOGR.pcaBiovar)
BOGR.pcaBiovar$rotation #loadings
#BOGR.pc1 <- BOGR.pcaBiovar$rotation[,1][order(BOGR.pcaBiovar$rotation[,1])]
#BOGR.pc2 <- BOGR.pcaBiovar$rotation[,2][order(BOGR.pcaBiovar$rotation[,2])]
## ----------------------------------------


## Models
## Use 'top' biovars as determined by PCA
## The following bioclim vars selected due to high loadings in orthogonal directions along PC1 and PC2 
##BIO4, BIO5, BIO11, BIO12, BIO17

BOGR.ht.bv.mod <- lmer(Length_cm_20220801 ~ scale(bio4) + scale(bio5) + scale(bio11) + scale(bio12) + scale(bio17)
                         + (1|Block), data=BOGR.cl)
Anova(BOGR.ht.bv.mod)
plot(allEffects(BOGR.ht.bv.mod))
plot(predictorEffects(BOGR.ht.bv.mod))
 

BOGR.inf.bv.mod <- glmer(NumInf ~ scale(bio4) + scale(bio5) + scale(bio11) + scale(bio12) + scale(bio17)
                         + (1|Block), data=BOGR.cl, family=poisson(link="log"))
Anova(BOGR.inf.bv.mod)
plot(allEffects(BOGR.inf.bv.mod))
plot(predictorEffects(BOGR.inf.bv.mod))

## Is negative binomial an appropriate model form?
hist(log(BOGR.cl$DaysToFlwr))
BOGR.df.bv.mod <- glmer.nb(DaysToFlwr ~ scale(bio4) + scale(bio5) + scale(bio11) + scale(bio12) + scale(bio17)
                       + (1|Block), data=BOGR.cl)
Anova(BOGR.df.bv.mod)
plot(allEffects(BOGR.df.bv.mod))
plot(predictorEffects(BOGR.df.bv.mod))
## ---------------------------------------


## Visualize raw data relationships with bioclimate vars
plot(BOGR.cl$Length_cm_20220801 ~ BOGR.cl$bio17)
plot(BOGR.cl$Length_cm_20220801 ~ BOGR.cl$bio5)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio5)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio11)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio4)
plot(BOGR.cl$DaysToFlwr ~ BOGR.cl$bio11)


## Calculate trait means
BOGR.summ <- BOGR.cl %>% group_by(Source) %>% 
             dplyr::summarise(DaysToFlwr_MN=mean(DaysToFlwr, na.rm=TRUE), DaysToFlwr_SE=calcSE(DaysToFlwr),
                   Height_MN=mean(Length_cm_20220801, na.rm=TRUE), Height_SE=calcSE(Length_cm_20220801),
                   Inf_MN=mean(NumInf, na.rm=TRUE), Inf_SE=calcSE(NumInf), 
                   Growth_MN=mean(GrwthRate_Relative, na.rm=TRUE), Growth_SE=calcSE(GrwthRate_Relative), 
                   BIO4=mean(bio4, na.rm=TRUE), BIO5=mean(bio5, na.rm=TRUE),BIO11=mean(bio11, na.rm=TRUE),
                   BIO12=mean(bio12, na.rm=TRUE), BIO17=mean(bio17, na.rm=TRUE))
BOGR.summ <- left_join(BOGR.summ, BOGR.SdZn, by="Source")

#Loop over columns with biovars
colnames(BOGR.summ)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
for (bb in 10:14) { 
  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Height_MN, col=BOGR.summ$SdZnColful, pch=19, cex=1.2, main="B. gracilis", 
     cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Height (cm)", cex.lab=1.5, cex.axis=1.1)
  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Height_MN-BOGR.summ$Height_SE, as.vector(t(BOGR.summ[,bb])), 
         BOGR.summ$Height_MN+BOGR.summ$Height_SE, 
         angle=90, col=BOGR.summ$SdZnColful, code=3, length=0, lwd=1.6) }
plot.new()
legend("center", unique(BOGR.meds$SdZnAbbrev[order(BOGR.meds$SdZnOrder)]), 
       col=unique(BOGR.meds$SdZnColful[order(BOGR.meds$SdZnOrder)]), cex=0.95, pch=19)


#for (bb in 10:14) { 
#  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Growth_MN, col=BOGR.summ$SdZnColful, pch=19, cex=1.2, main="B. gracilis", 
#       cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Relative plant growth", cex.lab=1.5, cex.axis=1.1)
#  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Growth_MN-BOGR.summ$Growth_SE, as.vector(t(BOGR.summ[,bb])), 
#         BOGR.summ$Growth_MN+BOGR.summ$Growth_SE, 
#         angle=90, col=BOGR.summ$SdZnColful, code=3, length=0, lwd=1.6) }


for (bb in 10:14) { 
  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Inf_MN, col=BOGR.summ$PopCol, pch=19, cex=1.2, main="B. gracilis", 
       cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Number of inflorescences", cex.lab=1.5, cex.axis=1.1)
  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Inf_MN-BOGR.summ$Inf_SE, as.vector(t(BOGR.summ[,bb])), 
         BOGR.summ$Inf_MN+BOGR.summ$Inf_SE, 
         angle=90, col=BOGR.summ$PopCol, code=3, length=0, lwd=1.6) }

for (bb in 10:28) { 
  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$DaysToFlwr_MN, col=BOGR.summ$PopCol, pch=19, cex=1.2, main="B. gracilis", 
       cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Days until first flower", cex.lab=1.5, cex.axis=1.1)
  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$DaysToFlwr_MN-BOGR.summ$DaysToFlwr_SE, as.vector(t(BOGR.summ[,bb])), 
         BOGR.summ$DaysToFlwr_MN+BOGR.summ$DaysToFlwr_SE, 
         angle=90, col=BOGR.summ$PopCol, code=3, length=0, lwd=1.6) }

#plot(BOGR.summ$Latitude, BOGR.summ$Height_MN, col=BOGR.summ$SdZnCol, pch=19, cex=1.25, main="",
#     cex.main=1.5, xlab="Latitude", ylab="Plant Height (cm)", cex.lab=1.5, cex.axis=1.1)
#arrows(BOGR.summ$Latitude, BOGR.summ$Height_MN-BOGR.summ$Height_SE, BOGR.summ$Latitude, BOGRsumm$Height_MN+BOGR.summ$Height_SE,
#       angle=90, col=BOGR.summ$SdZnCol, code=3, length=0, lwd=1.5)
## ---------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------








