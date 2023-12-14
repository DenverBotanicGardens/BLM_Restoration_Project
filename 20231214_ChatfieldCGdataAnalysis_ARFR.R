## April Goebl
## Script started 2023-12-14 (modified from 20231122_ChatfieldCGdataAnalysis_BOGR)
## BLM Restoration project at Denver Botanic Gardens
## Analyze ARFR data from Chatfield Common Garden  


rm(list=ls())
#dev.off()


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
ARFR <- read.csv(file="Chatfield/20230616_ChatfieldData2022_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
ARFR.SdZn <- read.csv(file="AGoebl/Seeds/20230823_ARFR_LatLong.csv", sep=",", header=TRUE, dec=".")
ARFR.biovar <- readRDS("AGoebl/Seeds/20230814_ARFR_BiovarsAvg1980_2021")
## ----------------------------------------------------------------------------------------------




## ARFR - DATA CLEAN UP ---------------------------------------------
str(ARFR)
ARFR$Source <- as.factor(ARFR$Source)

##If OrigPltSurv_20220527 = 0 & plt not replaced (N), ignore data for this plt, i.e. future surv should be NA, not 0
ARFR.cl <- ARFR[ARFR$OrigPltSurvival_20220527==1 | (ARFR$OrigPltSurvival_20220527==0 & ARFR$Replaced_YorN_20220531=="Y"),]

#Note: Don't use OrigPltSurv_20220527 data in days to mort or other field analyses,
#this surv may not correspond to plt names in Source/Pop col (may correspond to orig planted or assigned)

#Some plts that weren't dead & had sz measure on 5/27 were replaced. In this case, don't use length_0527.
#If length_0527 is not NA and Replaced = Y, then ignore length_527 as it doesn't correspond to plt names in source col
ARFR.cl$Length_cm_20220527[!is.na(ARFR.cl$Length_cm_20220527) & ARFR.cl$Replaced_YorN_20220531=="Y"] <- NA

#If not replaced, but died before planting, don't use subsequent surv data
#ARFR.cl[!is.na(ARFR.cl$DateMortalityObservedPreTransplant),] #All plts that died before planting were replaced



## Checks 
#Some plts were dead that were selected to be harvested, check? 
ARFR.Coll <- ARFR.cl[!is.na(ARFR.cl$Harvest_20221014) | !is.na(ARFR.cl$Harvest_20221110),]
#Current data sheet only have "Harvest" marked for plts there were alive
ARFR.Coll[is.na(ARFR.Coll$AGB_MinusBag),] #Two harvested plts do not have final AGB: 1107 and 1274. **Look in lab
ARFR.MissinfBM <- ARFR.Coll[is.na(ARFR.Coll$InfBM_Wobag_g),]
ARFR.MissinfBM$ID[ARFR.MissinfBM$Phenology_20220922==3] #Looks for these plts and get inf weights

## ** Add checks listed in BOGR ** 
## Check that surv is only 1, 0 and maybe NA
ARFR.cl[ARFR.cl$Survival_20220622 < 0 | ARFR.cl$Survival_20220622 > 1,]

#Check that pheno, surv, numInf are only integers
ARFR.cl[ARFR.cl$Survival_20220622 - floor(ARFR.cl$Survival_20220622) != 0,]
ARFR.cl[ARFR.cl$Survival_20220922 - floor(ARFR.cl$Survival_20220922) != 0,]
ARFR.cl[ARFR.cl$Phenology_20220715 - floor(ARFR.cl$Phenology_20220715) != 0 & 
          !is.na(ARFR.cl$Phenology_20220715),]

# ** Check that length is only numeric
# ** Check that if surv=0 for a given date, there are no phenology or height values for that date

## *** Need to work on this ****
#Check that once zero in surv on 6/22 or later, stays zero (if becomes 1 later, could be data entry error)
## CONSOLIDATE SURVIVAL DATA
ARFR.Surv <- ARFR.cl %>% dplyr::select(c(starts_with("Survival_")))
#for (rr in 1:nrow(ARFR.Surv)) {
#  for (cc in 1:(ncol(ARFR.Surv)-1)) {
#    if(ARFR.Surv[rr,cc]==0 & ARFR.Surv[rr,cc+1]==1) {
#     ARFR.Surv$Check[rr] <- "Y"
#    }
#  }
#}
ARFR.Surv %>% filter_all(any_vars(.==0))
#ARFR.cl %>% filter_all(starts_with("Survival_")==0)
## ----------------------------------------------------------------------------------------------




## ARFR - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ARFR.SdZn$Source <- str_replace(ARFR.SdZn$SOURCE_CODE, "4-SOS", "")
ARFR.biovar$Source <- str_replace(ARFR.biovar$Pop, "4-SOS", "")

## Add Population name abbreviation column
ARFR.SdZn$PopAbbrev[grepl("ARFR-AZ930-423-NAVAJO-18", ARFR.SdZn$Source)] = "A.AZ.1"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-AZ930-422-NAVAJO-18", ARFR.SdZn$Source)] = "A.AZ.2"    
ARFR.SdZn$PopAbbrev[grepl("ARFR-NM930N-66-11", ARFR.SdZn$Source)] = "A.NM.1"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY930-44-LASANIMAS-13", ARFR.SdZn$Source)] = "A.CO.1"       
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-314-JEFFERSON-12", ARFR.SdZn$Source)] = "A.CO.2"  
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-316-JEFFERSON-12", ARFR.SdZn$Source)] = "A.CO.3"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-UT080-109-UINTAH-12", ARFR.SdZn$Source)] = "A.UT.1"    
ARFR.SdZn$PopAbbrev[grepl("ARFR-CO932-294-11", ARFR.SdZn$Source)] = "A.CO.4"   
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY040-71-10", ARFR.SdZn$Source)] = "A.WY.1"       
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY050-151-FREMONT-16", ARFR.SdZn$Source)] = "A.WY.2" 
ARFR.SdZn$PopAbbrev[grepl("ARFR-WY050-49-FREMONT-12", ARFR.SdZn$Source)] = "A.WY.3"  

## Edit column names for biovariables
colnames(ARFR.biovar) <- c("Pop","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                           "bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                           "bio17","bio18","bio19","Source")

## Add colour columns that corresponds to pop or seed zone
#SdZn.list <- unique(ARFR.SdZn$SdZone) #If coloring by seed zone, change to match BOGR
#ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "powderblue"        #semi-humid, cold
#ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen"     #semi-humid, cool
#ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen4"    #semi-humid, warm
#ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "darkgoldenrod1"   #semi-arid, cold
#ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "orange3"         #semi-arid, cool
#ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "tomato2"         #semi-arid, warm

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
#BOGR.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "A"    #semi-humid, cool
#BOGR.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "B"    #semi-humid, warm
#BOGR.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "C"    #semi-humid, v.warm
#BOGR.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "D"   #semi-arid, cool
#BOGR.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "E"   #semi-arid, warm
#BOGR.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "F"   #semi-arid, v.warm
#BOGR.SdZn$SdZnOrder[grepl("25 - 30 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "G"   #semi-arid, hot
#BOGR.SdZn$SdZnOrder[grepl("30 - 35 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "H"   #semi-arid, v.hot
#BOGR.SdZn$SdZnOrder[grepl("5 - 10 Deg. F. / 12 - 30", BOGR.SdZn$seed_zone)] = "I"   #arid, cold
## ----------------------------------------------------------------------------------------------




## ARFR - COMBINE DATA TYPES --------------------------------------------
ARFR.cl <- left_join(ARFR.cl, ARFR.SdZn, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.biovar, by="Source")
## ----------------------------------------------------------------------



## ARFR - CONSOLIDATE PHENOLOGY DATA AND DO CHECKS ---------------------------------------------
## Add from 20230129_ChatfieldCGdataAnalysis and BOGR script if want to include.
## Exclude for now due to lack of variation and questionable data (inability to distinguish buds etc.)
## ---------------------------------------------------------------



## ARFR - ADD GROWTH RATE VARIABLES ------------------------------
ARFR.cl$GrwthRate_Specific <- log(ARFR.cl$Length_cm_20220726/ARFR.cl$Length_cm_20220622)
ARFR.cl$GrwthRate_Absolute <- ARFR.cl$Length_cm_20220726-ARFR.cl$Length_cm_20220622
ARFR.cl$GrwthRate_Relative <- (ARFR.cl$Length_cm_20220726-ARFR.cl$Length_cm_20220622)/ARFR.cl$Length_cm_20220622

## ** Look at early vs late growth if can use pre-replacement early height measurement 20220527 ** 
## ---------------------------------------------------------------



## LOOK AT REPRO BM DATA 
## Add any mods? 
ARFR.cl$InfBM_Wobag_g
ARFR.cl$InfBM_Wobag_g[!is.na(ARFR.cl$InfBM_Wobag_g)]
length(ARFR.cl$InfBM_Wobag_g[!is.na(ARFR.cl$InfBM_Wobag_g)])
length(ARFR.cl$InfBM_Wbag[!is.na(ARFR.cl$InfBM_Wbag)])
hist(ARFR.cl$InfBM_Wobag_g)
str(ARFR.cl$InfBM_Wobag_g)
## ---------------------------------------------------------------



## COULD ADD AGB VARIABLES AS WELL? *
## For now, I think height is better since more data and correlated to AGB




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
## Order by average size 
ARFR.htByMed <- with(ARFR.cl, reorder(Source, Length_cm_20220726, median, na.rm=TRUE))
ARFR.meds <- ARFR.cl %>% group_by(Source) %>% 
             dplyr::summarise(Height_MD=median(Length_cm_20220726,na.rm=TRUE), AGB_MD=median(AGB_MinusBag,na.rm=TRUE),
             ReproBM_MD=median(InfBM_Wobag_g,na.rm=TRUE), GrowthRe_MD=median(GrwthRate_Relative,na.rm=TRUE),
             GrowthSp_MD=median(GrwthRate_Specific,na.rm=TRUE), GrowthAb_MD=median(GrwthRate_Absolute,na.rm=TRUE))
ARFR.meds <- left_join(ARFR.meds, ARFR.SdZn, by="Source")


## Boxplots of raw data 
par(mfrow=c(1,1))

## Size
ARFR.meds <- ARFR.meds[order(ARFR.meds$Height_MD),] #Order by median 
boxplot(Length_cm_20220726 ~ ARFR.htByMed, data=ARFR.cl,
        xlab=NA, ylab="Plant height (cm)", cex.lab=1.25,
        cex.axis=0.95, names=ARFR.meds$PopAbbrev, las=2,
        main="Artemisia frigida", cex.main=1.5, col=ARFR.meds$PopCol.x)

## UPDATE ONCE ADD POP ORDER **
#legend("topleft", unique(ARFR.meds$SdZnAbbrev[order(ARFR.meds$SdZnOrder)]), 
#       col=unique(ARFR.meds$PopCol[order(ARFR.meds$SdZnOrder)]), cex=1.95, pch=19)


## UPDATE TO ORDER BY AVG SIZE ********************************
## Growth rate(s) ** Add time interval to growth rate calcs? **
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Growth_MD),]
#boxplot(GrwthRate_Relative ~ ARFR.grrByMed, data=ARFR.cl, las=2,
#        xlab="Population", ylab="Plant relative growth", cex.lab=1.25, names=ARFR.meds$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=ARFR.meds$SdZnColful)

#ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$GR_AB),]
#boxplot(GrwthRate_Absolute ~ ARFR.graByMed, data=ARFR.cl, las=2,
#        xlab="Population", ylab="Plant Absolute growth", cex.lab=1.25, names=ARFR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=ARFR.SdZn$PopCol)

#ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ ARFR.grsByMed, data=ARFR.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=ARFR.SdZn$PopCol)

## Repro
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Inf_MD),]
#boxplot(NumInf ~ ARFR.infByMed, data=ARFR.cl, las=2,
#        xlab="Population", ylab="Number of inflorescences", cex.lab=1.25, names=ARFR.meds$PopAbbrev,
#        cex.axis=0.7, main="Artemisia frigida", cex.main=1.5, col=ARFR.meds$SdZnColful)
## ---------------------------------------------------


