## April Goebl
## Script started 2023-01-29
## BLM Restoration project at Denver Botanic Gardens
## Analyze data from Chatfield Common Garden  




rm(list=ls())



## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
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
calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------



## SET WORKING DIRECTORY --------------------------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - BLM-Grassland")
setwd("C:/Users/april/Denver Botanic Gardens/Conservation - BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
PEVI <- read.csv(file="Chatfield/20230303_ChatfieldData_PEVI.csv", sep=",", header=TRUE, dec=".")
ARFR <- read.csv(file="Chatfield/20230616_ChatfieldData2022_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
ERNA <- read.csv(file="Chatfield/20230302_ChatfieldData_ERNA.csv", sep=",", header=TRUE, dec=".")
BOGR <- read.csv(file="Chatfield/20230301_ChatfieldData_BOGR.csv", sep=",", header=TRUE, dec=".")

#ERNA.SdZn <- read.csv(file="AGoebl/Seeds/Join_Features_to_SeedZones/ERNA_latLong_SdZone.csv", sep=",", header=TRUE, dec=".")
ARFR.SdZn <- read.csv(file="AGoebl/Seeds/20220929_ARFR_LatLong.csv", sep=",", header=TRUE, dec=".")
ARFR.ppt <- read.csv(file="AGoebl/Seeds/20230311_ARFR_pptAnnual.csv", sep=",", header=TRUE, dec=".")
ARFR.tmin <- read.csv(file="AGoebl/Seeds/20230311_ARFR_tminWinter.csv", sep=",", header=TRUE, dec=".")
## ------------------------------------------------------------------------------------------------



## DATA CLEAN-UP ---------------------------------------

## BOGR ------------------------------------------------
#BOGR <- BOGR[1:1135,]
#If OrigPltSurv_20220531 = 0 and plt not replaced (N),
#ignore subsequent data for this plant, i.e. future survival should be NA and not 0
#Don't use OrigPltSurv_20220531 data in days to mort or other field analyses, 
#this surv may not correspond to plt names in Code/Pop column (may correspond to orig planted or assigned)
#Check that all numinf_coll_0927 are entered in numinf_0927
#Check that once zero in surv on 6/27 or later, stays zero (cannot become 1 later)
#Check that phenology only increases or stays the same
#Check that if surv=0 for a given date, there are no phenology or height values for that date
#Check that phenology and surv are only integers and length is only numeric
#If numInf col has value greater than 0, phenol for corresponding date should be >2
#Pheno cols should only be b/w 2-6 (check upper value)
#NumInf cols should only be integers
#Surv should only be 1, 0 and maybe NA
#Make a combine num inf col that combines numinf_927 with numinfCol1014 and 1108 
#if the latter 2 are from blocks 7-10

#Trying to figure out for block 3 how measured heights from 5/31 correspond to plts names in Source col
#given replacements made on 6/6. There was an error in data entry (one row off; a non-BOGR measured at the end of a row/ block) around 5/31-6/2
#but correction to the datasheet was not done until around 6/28. 



## ARFR DATA CLEAN UP ---------------------------------------------
ARFR$Source <- as.factor(ARFR$Source)

##If OrigPltSurv_20220527 = 0 & plt not replaced (N), ignore data for this plt, i.e. future surv should be NA, not 0
ARFR.cl <- ARFR[ARFR$OrigPltSurvival_20220527==1 | (ARFR$OrigPltSurvival_20220527==0 & ARFR$Replaced_YorN_20220531=="Y"),]

#Note: Don't use OrigPltSurv_20220527 data in days to mort or other field analyses,
#this surv may not correspond to plt names in Source/Pop col (may correspond to orig planted or assigned)

#Some plts that weren't dead & had sz measure on 5/27 were replaced. In this case, don't use length_0527.
#If length_0527 is not NA and Replaced = Y, then ignore length_527 as it doesn't correspond to plt names in source col
ARFR.cl$Length_cm_20220527[!is.na(ARFR.cl$Length_cm_20220527) & ARFR.cl$Replaced_YorN_20220531=="Y"] <- NA

#If not replaced, but died before planting, don't use subsequent surv data
ARFR.cl[!is.na(ARFR.cl$DateMortalityObservedPreTransplant),] #All plts that died before planting were replaced

## *** Need to work on this more ****
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


#Some plts were dead that were selected to be harvested, check? 

#Excel formula for surv cols: =IF(C2="",1,0)

#Find bag size for ARFR ID 864. And look at entire AGB col for any missing



## ARFR DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ARFR.SdZn$Source <- str_replace(ARFR.SdZn$SOURCE_CODE, "4-SOS", "")
ARFR.ppt$Source <- str_replace(ARFR.ppt$Source, "4-SOS", "")
ARFR.tmin$Source <- str_replace(ARFR.tmin$Code, "4-SOS", "")

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


## COMBINE DATA TYPES --------------------------------------------
ARFR.cl <- left_join(ARFR.cl, ARFR.SdZn, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.ppt, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.tmin, by="Source")


## Add a colour column that corresponds to seed zone
#SdZn.list <- unique(ARFR.SdZn$SdZone_approx)
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 3 - 6", ARFR.SdZn$SdZone_approx)] = "darkseagreen"     #semi-humid, cool
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 3 - 6", ARFR.SdZn$SdZone_approx)] = "darkseagreen4"    #semi-humid, warm
ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 6 - 12", ARFR.SdZn$SdZone_approx)] = "darkgoldenrod1"   #semi-arid, cold
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 6 - 12", ARFR.SdZn$SdZone_approx)] = "orange2"         #semi-arid, warm
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 12 - 30", ARFR.SdZn$SdZone_approx)] = "tomato2"        #arid, cool
## ---------------------------------------------------------------


## CONSOLIDATE PHENOLOGY DATA
#ARFR.Pheno <- ARFR %>% dplyr::select(c(Source, Treatment, Block,SdZone_approx,PopAbbrev, Lat, starts_with("Phenology")))
ARFR.PlantingDate <- as.Date("2022-05-03")
ARFR.PhenoCol.List <- colnames(ARFR)[grepl("Phenology*", colnames(ARFR))] #Obtain phenology column names
ARFR.Pheno.List <- str_replace(ARFR.PhenoCol.List, "Phenology_", "")      #Obtain just date from phenology columns
ARFR.Pheno.List <- as.Date(ARFR.Pheno.List, "%Y%m%d")
ARFR.DaysToFlwr <- ARFR.Pheno.List - ARFR.PlantingDate                         #Calculate days from planting to each phenology survey 

## Loop over each phenology column & enter the num days since planting when a 3 (flower or later phenophase) first appears
ARFR$DaysToFlwr <- NA
for (pp in 1:length(ARFR.PhenoCol.List)) {
  ARFR$DaysToFlwr[ARFR[,ARFR.PhenoCol.List[pp]]>=3 & is.na(ARFR$DaysToFlwr)] <- as.integer(ARFR.DaysToFlwr)[pp]
}
## ---------------------------------------------------------------


## TEST FOR TREATMENT EFFECT -------------------------------------
hist(log(ARFR$AGB_MinusBag))
#ARFR.tx.mod <- aov(ARFR$AGB_MinusBag ~ ARFR$Source + ARFR$Treatment)
ARFR.tx.mod <- lmer(log(AGB_MinusBag) ~ Source + Treatment + (1|Block), data=ARFR)
summary(ARFR.tx.mod)
ARFR.pop.mod <- lmer(log(AGB_MinusBag) ~ Source + (1|Block), data=ARFR)

models <- list(ARFR.tx.mod, ARFR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment 


hist(log(ARFR$Length_cm_20220726))
#ARFR.tx.mod <- aov(log(ARFR$Length_cm_20220726) ~ ARFR$Source + ARFR$Treatment)
ARFR.tx.mod <- lmer(log(ARFR$Length_cm_20220726) ~ Source + Treatment + (1|Block), data=ARFR)
summary(ARFR.tx.mod)
ARFR.pop.mod <- lmer(log(ARFR$Length_cm_20220726) ~ Source + + (1|Block), data=ARFR)

models <- list(ARFR.tx.mod, ARFR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment 
## -----------------------------------------------------------------



## VISUALIZE RAW DATA -----------------------------------------------------------------------------
## ARFR --------------
## Order populations for plotting 

## EITHER --- Order by seed zone -----
## Order based on aridity & temp (cool & wet to hot & dry), w/in each cat, order from S to N [*CONFIRM*])
ARFR.SdZn$Source <- factor(ARFR.SdZn$Source, levels=c("ARFR-UT080-109-UINTAH-12","ARFR-CO932-316-JEFFERSON-12",
                                                      "ARFR-CO932-314-JEFFERSON-12","ARFR-NM930N-66-11",
                                                      "ARFR-WY930-44-LASANIMAS-13","ARFR-CO932-294-11","ARFR-WY040-71-10",
                                                      "ARFR-WY050-49-FREMONT-12","ARFR-AZ930-422-NAVAJO-18",
                                                      "ARFR-AZ930-423-NAVAJO-18", "ARFR-WY050-151-FREMONT-16"))
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Source),]
ARFR$Source <- factor(ARFR$Source, levels=ARFR.SdZn$Source)


## OR --- Order by latitude ----
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Lat),]
ARFR$Source <- factor(ARFR$Source, levels=ARFR.SdZn$Source)


## Boxplots of raw data 
boxplot(Length_cm_20220726 ~ Source, data=ARFR,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        cex.axis=1.1, names=ARFR.SdZn$PopAbbrev,
        main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

## Add time interval to growth rate calcs?
## Specific Growth Rate (takes into account exponential growth of some plants)
boxplot(log(Length_cm_20220726/Length_cm_20220622) ~ Source, data=ARFR.cl, names=ARFR.SdZn$PopAbbrev,
        xlab="Population", ylab="Plant growth", cex.lab=1.5, ylim=c(0,3),
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## Absolute Growth Rate
boxplot((Length_cm_20220726-Length_cm_20220622) ~ Source, data=ARFR.cl, names=ARFR.SdZn$PopAbbrev,
        xlab="Population", ylab="Plant growth", cex.lab=1.5, ylim=c(0,60),
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## Relative Growth Rate
boxplot(((Length_cm_20220726-Length_cm_20220622)/Length_cm_20220622) ~ Source, data=ARFR.cl, names=ARFR.SdZn$PopAbbrev,
        xlab="Population", ylab="Plant growth", cex.lab=1.5, ylim=c(0,15),
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)


boxplot(AGB_MinusBag ~ Source, data=ARFR.cl,
        xlab="Population", ylab="Above-ground biomass (g)", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

#boxplot(DaysToFlwr ~ Source, data=ARFR,
#        xlab="Population", ylab="Days to first flower", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
#        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

boxplot(InfBM_Wbag ~ Source, data=ARFR.cl,
        xlab="Population", ylab="Reproductive biomass (g)", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## ---------------------------------------   


## ESTIMATE VARIATION B/W POPULATIONS -------------------------------------------------------------
## ARFR
## Final size (height and biomass)
ARFR.htMod <- lmer(log(Length_cm_20220726) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.ht <- emmeans(ARFR.htMod, specs = pairwise ~ Source)
summary(ARFR.pair.ht)
plot(ARFR.pair.ht)

ARFR.bmMod <- lmer(log(AGB_MinusBag) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.bm <- emmeans(ARFR.bmMod, specs = pairwise ~ Source)
summary(ARFR.pair.bm)
plot(ARFR.pair.bm)

## Plant growth
ARFR.grMod <- lmer(log(Length_cm_20220726/Length_cm_20220622) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.gr <- emmeans(ARFR.grMod, specs = pairwise ~ Source)
summary(ARFR.pair.gr)
plot(ARFR.pair.gr)

## Repro BM
ARFR.rbmMod <- lmer(log(InfBM_Wbag) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.rbm <- emmeans(ARFR.rbmMod, specs = pairwise ~ Source)
summary(ARFR.pair.bm)
plot(ARFR.pair.rbm)


## ESTIMATE VARIATION WITHIN POPULATIONS -------------------------------------------------------------
## ARFR -------------
## Order populations for plotting 
## EITHER --- Order by seed zone -----
## Order based on aridity & temp (cool & wet to hot & dry), w/in each cat, order from S to N [*CONFIRM*])
ARFR.SdZn$Source <- factor(ARFR.SdZn$Source, levels=c("ARFR-UT080-109-UINTAH-12","ARFR-CO932-316-JEFFERSON-12",
                                                      "ARFR-CO932-314-JEFFERSON-12","ARFR-NM930N-66-11",
                                                      "ARFR-WY930-44-LASANIMAS-13","ARFR-CO932-294-11","ARFR-WY040-71-10",
                                                      "ARFR-WY050-49-FREMONT-12","ARFR-AZ930-422-NAVAJO-18",
                                                      "ARFR-AZ930-423-NAVAJO-18", "ARFR-WY050-151-FREMONT-16"))
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Source),]
ARFR$Source <- factor(ARFR$Source, levels=ARFR.SdZn$Source)

## OR --- Order by latitude ----
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Lat),]
ARFR$Source <- factor(ARFR$Source, levels=ARFR.SdZn$Source)



## Calculate the coefficient of variation 
ARFR.ht.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Height_CV=cv(Length_cm_20220726, na.rm=TRUE))
ARFR.gr.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Growth_CV=cv(log(Length_cm_20220726/Length_cm_20220622), na.rm=TRUE))
ARFR.bm.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Biomass_CV=cv(AGB_MinusBag, na.rm=TRUE))
#ARFR.dtf.cv <- ARFR %>% group_by(Source) %>% summarise(DaysToFlwr_CV=cv(DaysToFlwr, na.rm=TRUE))
ARFR.rbm.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Rbiomass_CV=cv(InfBM_Wbag, na.rm=TRUE))
## ---------



## Barplots
barplot(ARFR.ht.cv$Height_CV, xlab="Population", ylab="CV in plant height", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.gr.cv$Growth_CV, xlab="Population", ylab="CV in plant growth", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.bm.cv$Biomass_CV, xlab="Population", ylab="CV in above-ground biomass", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

#barplot(ARFR.dtf.cv$DaysToFlwr_CV, xlab="Population", ylab="CV in days to first flower", cex.lab=1.5, cex.axis=1.1, 
#        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.rbm.cv$Rbiomass_CV, xlab="Population", ylab="CV in reproductive biomass", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## ---------



## LOOK AT RELATIONSHIP BETWEEN TRAIT CV AND AVERAGE ---------------------------
ARFR.ht.mn <- ARFR.cl %>% group_by(Source) %>% summarise(Height_MN=mean(Length_cm_20220726, na.rm=TRUE))
ARFR.gr.mn <- ARFR.cl %>% group_by(Source) %>% summarise(Growth_MN=mean(log(Length_cm_20220726/Length_cm_20220622), na.rm=TRUE))
ARFR.bm.mn <- ARFR.cl %>% group_by(Source) %>% summarise(Biomass_MN=mean(AGB_MinusBag, na.rm=TRUE))
ARFR.rbm.mn <- ARFR.cl %>% group_by(Source) %>% summarise(Rbiomass_MN=mean(InfBM_Wbag, na.rm=TRUE))

plot(ARFR.ht.cv$Height_CV, ARFR.ht.mn$Height_MN)
plot(ARFR.gr.cv$Growth_CV, ARFR.gr.mn$Growth_MN)
plot(ARFR.bm.cv$Biomass_CV, ARFR.bm.mn$Biomass_MN)
plot(ARFR.rbm.cv$Rbiomass_CV, ARFR.rbm.mn$Rbiomass_MN)




## EVALUATE RELATIONSHIPS B/W TRAITS AND SOURCE CLIMATE -------------------------------------------
## ARFR -------------
## Visualize raw data
ARFR.summ <- ARFR %>% group_by(Source) %>% 
  summarise(AGB_MEAN=mean(AGB_MinusBag, na.rm=TRUE), AGB_SE=calcSE(AGB_MinusBag), 
            PPT=mean(Mean_Annual_Ppt, na.rm=TRUE), TMIN=mean(Mean_MinWinter_Temp, na.rm=TRUE))

ARFR.summ <- left_join(ARFR.summ, ARFR.SdZn, by="Source")

## Look at correlations in climate, latitude, etc.
plot(ARFR.summ$TMIN, ARFR.summ$Lat, pch=19)
summary(lm(ARFR.summ$TMIN ~ ARFR.summ$Lat))


plot(ARFR.summ$PPT, ARFR.summ$AGB_MEAN, col=ARFR.SdZn$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida", 
     cex.main=1.5, xlab="Mean annual precipitation", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$PPT, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$PPT, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE, 
       angle=90, col=ARFR.SdZn$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$TMIN, ARFR.summ$AGB_MEAN, col=ARFR.SdZn$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida", 
     cex.main=1.5, xlab="Mean winter temperature", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$TMIN, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$TMIN, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE, 
       angle=90, col=ARFR.SdZn$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$Lat, ARFR.summ$AGB_MEAN, col=ARFR.SdZn$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida",
     cex.main=1.5, xlab="Latitude", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$Lat, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$Lat, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE,
       angle=90, col=ARFR.SdZn$SdZnCol, code=3, length=0, lwd=1.5)


## Models
ARFR$Ppt.sq <- (ARFR$Mean_Annual_Ppt)^2
ARFR.bm.mod <- lmer(log(AGB_MinusBag) ~ Mean_Annual_Ppt + Ppt.sq + Mean_MinWinter_Temp + (1|Block), data=ARFR)
summary(ARFR.bm.mod)
Anova(ARFR.bm.mod)

plot(allEffects(ARFR.bm.mod))
plot(predictorEffects(ARFR.bm.mod))
## -------------------------------------------------------------------------------



## EMMEANS -----------------------------------------------------------------------
ARFR.htMod <- lmer(log(Length_cm_20220726) ~ Source + (1|Block), data = ARFR)
ARFR.ht.em <- emmeans(ARFR.htMod, specs = pairwise ~ Source)
#summary(ARFR.ht.em)
#str(ARFR.ht.em)
ARFR.ht.pw <- as.data.frame(ARFR.ht.em$contrasts)
#ARFR.ht.pw$Pop2 <- ARFR.ht.pw$contrast
#colnames(ARFR.ht.pw) <- c("Pop1", colnames(ARFR.ht.pw[2:7]))
#ARFR.ht.pw$Pop1 <- str_replace(ARFR.ht.pw$contrast, " - *", "")
ARFR.ht.pw$Pop1 <- gsub("\\ - .+?\\)", "", ARFR.ht.pw$contrast)
ARFR.ht.pw$Pop2 <- gsub("\\(.+?\\ - ", "", ARFR.ht.pw$contrast)
ARFR.ht.pw.pvals <- ARFR.ht.pw %>% dplyr::select(Pop1,Pop2, p.value)
ARFR.ht.pw.pvals <- as.matrix(acast(ARFR.ht.pw.pvals, Pop1~Pop2))
heatmap.2(ARFR.ht.pw.pvals, dendrogram="none", cexRow=0.7, cexCol=0.7)











## ERNA ----------------------------------------------
ERNA <- ERNA[ERNA$Replaced_YorN !="Y",]

#OrigPltLength_cm_Greenhouse_20220323 may not correspond to plt names in Source col 
#(may correspond to orig planted or assigned. This will be the case if the plant was replaced)'
#so probably don't use OrigPltLen in growth rates or other field analyses (especially if Replaced = Y)

#Seedlings planted 4/20-4/28
#First surv 5/18
#First replacements 5/18-5/19. Replaced if dead and missing and more of source available;
#not replaced if appeared dead (no green tissue) but not missing.
#Second surv/sz 6/8-6/10
#Second replacements 6/13. Replaced if dead on both 5/18 and 6/8. Not replaced if only dead 6/8.
#If 1 in both OrigPltSurv cols (Replaced=N or NA), then these and orig sizes correspond to plant named in Source.
#If 1 in first OrigSurv and 0 in second OrigSurv (Replaced=N or NA), then no data other than first OrigLen and first OrigSurv,
#and subsequent survs. The days to mort for these plants could be counted as the date of the second OrigSurv survey.
#If 0 in first OrigSurv and 0 in second OrigSurv (Replaced=Y if Source available), 
#then first OrigLen does not correspond to plant named and no data for second OrigLen. Subsequent data corresponds to plt named.
#If 0 in first OrigSurv and 1 in second OrigSurv (Replaced=Y or N), if Replaced=Y then first OrigLen does not correspond to plt named;
#and second OrigLen and subsequent data does correspond to plt named. If Replaced=N then first OrigSurv is wrong,
#but all other data (including first and second OrigLen and second OrigSurv) corresponds to plt named.
#Source should match Replacement (when value entered) minus unique plt id 

#Due to how confusing this is, it is probably easiest to ignore all plants with Replaced=Y (even if subsequent data corresponds
#to plt named, the data for these plts may vary do to being planted 1-2 months later).
#Also don't use first OrigSurv in days to mort or other field analyses, these deaths are due to transplant shock.
#Can use second OrigSurv in days to mort if first OrigSurv=1, ignore plt if first OrigSurv=0. If this plt wasn't replaced then it died from transplant.

#ERNA <- left_join(ERNA, ERNA.SdZn, by="Source")



## PEVI ------------------
#If dead on OrigSurv, ignore plt if not replaced. If replaced, use subsequent surv and sz data
#Don't use OrgiSurv for days to mort 




## VISUALIZE RAW DATA -----------------------------------------------------------------------------
par(mfrow=c(4,1))
par(mfrow=c(1,1))

## SIZE, GROWTH, DAYS TO FLOWER -------------------------------------------------------------------

## ERNA ------------------
#ERNA <- ERNA[ERNA$Code!="ERNA10-DBG-Chatfield-21",]
boxplot(Length_cm_20220915 ~ Source, data=ERNA,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        names=c("E.CO1", "E.CO2", "E.CO3", "E.CO4","E.CO5","E.CO6","E.CO7", "E.CO8", "E.NV1",
        "E.NV2", "E.NV3", "E.NV4","E.NV5","E.NV6", "E.UT1","E.UT2","E.UT3", "E.UT4",
        "E.WY1", "E.CO9"), cex.axis=0.8,
        main="Ericameria nauseosa", cex.main=1.5, col="palegoldenrod")

## ** UPDATE GROWTH FORMULA ******
boxplot((Length_cm_20220915/Length_cm_20220608) ~ Source, data=ERNA,
        xlab="Population", ylab="Plant growth", cex.lab=1.5,
        names=c("E.CO1", "E.CO2", "E.CO3", "E.CO4","E.CO5","E.CO6","E.CO7", "E.CO8", "E.NV1",
                "E.NV2", "E.NV3", "E.NV4","E.NV5","E.NV6", "E.UT1","E.UT2","E.UT3", "E.UT4",
                "E.WY1", "E.CO9"), cex.axis=0.8,
        main="Ericameria nauseosa", cex.main=1.5, col="palegoldenrod")



## BOGR ------------------
boxplot(Length_cm_20220801 ~ Source, data=BOGR,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
        "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1",
        "B.WY2", "B.WY3"), cex.axis=0.7,
        main="Bouteloua gracilis", cex.main=1.5, col="lightcyan")

## ** UPDATE GROWTH FORMULA ******
boxplot((Length_cm_20220801/Length_cm_20220627) ~ Source, data=BOGR,
        xlab="Population", ylab="Plant growth", cex.lab=1.5,
        names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
                "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1",
                "B.WY2", "B.WY3"), cex.axis=0.7,
        main="Bouteloua gracilis", cex.main=1.5, col="lightcyan")



## PEVI ---------------
boxplot(Length_cm_20220707 ~ Source, data=PEVI, ylim=c(0,10),
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        names=c("P.CO1", "P.CO2", "P.CO3", "P.CO4", "P.WY1", "P.WY2"),
        main="Penstemon virens", cex.main=1.5, col="thistle", cex.axis=1.25)
## ------------------------------------------------------------------------------------------------




## ESTIMATE VARIATION B/W POPULATIONS -------------------------------------------------------------
## BOGR
## Final size (height)
BOGR.htMod <- lmer(log(Length_cm_20220801) ~ Source + (1|Block), data = BOGR)
BOGR.ht.em <- emmeans(BOGR.htMod, specs = pairwise ~ Source)
summary(BOGR.ht.em)
plot(BOGR.ht.em)

BOGR.ht.pw <- as.data.frame(BOGR.ht.em$contrasts)
BOGR.ht.pwDiffs <- exp(BOGR.ht.pw$estimate)
max(BOGR.ht.pwDiffs)

## Plant growth
BOGR.grMod <- lmer(log(Length_cm_20220801/Length_cm_220531) ~ Source + (1|Block), data = BOGR)
BOGR.gr.em <- emmeans(BOGR.grMod, specs = pairwise ~ Source)
summary(BOGR.gr.em)
plot(BOGR.gr.em)

BOGR.gr.means <- as.data.frame(BOGR.gr.em$emmeans)
BOGR.gr.meanEsts <- exp(BOGR.gr.means$emmean)
max(BOGR.gr.meanEsts)
min(BOGR.gr.meanEsts)

BOGR.gr.pw <- as.data.frame(BOGR.gr.em$contrasts)
BOGR.gr.pwDiffs <- exp(BOGR.gr.pw$estimate)
max(BOGR.gr.pwDiffs)

## Survival
BOGR.survMod <- glmer(Survival_20220927 ~ Source + (1|Block), data = BOGR, family = binomial (link ="logit"))

## Days until flower

## Number of inflorescences
BOGR.inf.climate.lm <- glmer( num_inf_20220927 ~ Ppt_Annual_Z + min_wint_temp_Z + elev_Z + (1|Block), data = BOGR.crop, family = poisson (link ="log"))







## ESTIMATE VARIATION WITHIN POPULATIONS -------------------------------------------------------------
## ERNA ----------------
## Final size (height)


## BOGR ----------------
## Final size (height)
## Look at mean and var of CV across blocks for each source
BOGR.ht.cv <- NULL

for (vv in 1:10) {
  dat.b <- subset(BOGR, BOGR$Block==vv)
  BOGR.ht.cv <- rbind(BOGR.ht.cv, dat.b %>% group_by(Source) %>%
                        summarise(Height_CV=cv(Length_cm_20220801, na.rm=TRUE)))
}

## Test for significance in CV differences
BOGR.ht.CV.aov <- aov(Height_CV ~ Source, data=BOGR.ht.cv)
summary(BOGR.ht.CV.aov)

## Boxplot
boxplot(Height_CV ~ Source, data=BOGR.ht.cv,
        xlab="Population", ylab="CV in plant height", cex.lab=1.5,
        names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
                "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1",
                "B.WY2", "B.WY3"), cex.axis=0.7,
        main="Bouteloua gracilis", cex.main=1.5, col="lightcyan")

## Barplot
BOGR.ht.mn.se <- BOGR.ht.cv %>% group_by(Source) %>%
                 summarise(MN=mean(Height_CV, na.rm=TRUE), SE=calcSE(Height_CV))

plotCI(barplot(BOGR.ht.mn.se$MN, col="lightcyan", ylab="CV in plant height", xlab="Population",
               names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
                       "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1", "B.WY2", "B.WY3"), 
               cex.names=0.7, ylim=c(0,0.7), cex.lab=1.4, main="Bouteloua gracilis", cex.main=1.5),
               BOGR.ht.mn.se$MN, uiw=BOGR.ht.mn.se$SE, add=TRUE, pch=NA, sfrac=0)


## Plant growth    
## Look at mean and var of CV across blocks for each source
BOGR.gr.cv <- NULL

for (vv in 1:10) {
  dat.b <- subset(BOGR, BOGR$Block==vv)
  BOGR.gr.cv <- rbind(BOGR.gr.cv, dat.b %>% group_by(Source) %>%
                        summarise(Growth_CV=cv(Length_cm_20220801/Length_cm_20220627, na.rm=TRUE)))
}

## Test for significance in CV differences
BOGR.gr.CV.aov <- aov(Growth_CV ~ Source, data=BOGR.gr.cv)
summary(BOGR.gr.CV.aov)

## Boxplot
boxplot(Growth_CV ~ Source, data=BOGR.gr.cv,
        xlab="Population", ylab="CV in plant growth", cex.lab=1.5,
        names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
                "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1",
                "B.WY2", "B.WY3"), cex.axis=0.7,
        main="Bouteloua gracilis", cex.main=1.5, col="lightcyan")

## Barplot
BOGR.gr.mn.se <- BOGR.gr.cv %>% group_by(Source) %>%
            summarise(MN=mean(Growth_CV, na.rm=TRUE), SE=calcSE(Growth_CV))

plotCI(barplot(BOGR.gr.mn.se$MN, col="lightcyan", ylab="CV in plant growth", xlab="Population",
        names=c("B.AZ1","B.AZ2","B.AZ3","B.AZ4","B.AZ5","B.CO1","B.CO2","B.CO3", "B.CO4","B.CO5",
                "B.CO6","B.NM1", "B.NM2", "B.NM3", "B.NM4", "B.NM5", "B.UT1", "B.UT2","B.WY1", "B.WY2", "B.WY3"), 
        cex.names=0.7, ylim=c(0,0.56), cex.lab=1.4, main="Bouteloua gracilis", cex.main=1.5),
       BOGR.gr.mn.se$MN, uiw=BOGR.gr.mn.se$SE, add=TRUE, pch=NA, sfrac=0)
## -------------------------------------------------------------------------------






## EMMEANS 
BOGR.ht.em <- readRDS("20230318_BOGRhtPw.rds")
BOGR.ht.pw <- as.data.frame(BOGR.ht.em$contrasts)
BOGR.ht.pw$Pop1 <- gsub("\\ - .+?\\)", "", BOGR.ht.pw$contrast)
BOGR.ht.pw$Pop2 <- gsub("\\(.+?\\ - ", "", BOGR.ht.pw$contrast)
BOGR.ht.pw.pvals <- BOGR.ht.pw %>% dplyr::select(Pop1,Pop2, p.value)
BOGR.ht.pw.pvals <- as.matrix(acast(BOGR.ht.pw.pvals, Pop1~Pop2))
heatmap.2(BOGR.ht.pw.pvals, dendrogram="none", cexRow=0.7, cexCol=0.7)
## -------------------------------------------------------------------------------








## EXTRA CODE --------------------------------------------------------------------------------------------
#ARFR$Source <- factor(ARFR$Source, levels=c("ARFR-UT080-109-UINTAH-12","ARFR-CO932-316-JEFFERSON-12",
#                                            "ARFR-CO932-314-JEFFERSON-12","ARFR-NM930N-66-11",
#                                            "ARFR-WY930-44-LASANIMAS-13","ARFR-CO932-294-11","ARFR-WY040-71-10",
#                                            "ARFR-WY050-49-FREMONT-12","ARFR-AZ930-422-NAVAJO-18",
#                                            "ARFR-AZ930-423-NAVAJO-18", "ARFR-WY050-151-FREMONT-16"))

## Previous default order
#names=c("A.AZ1", "A.AZ2", "A.CO1", "A.CO2","A.CO3","A.NM1","A.UT1", 
#        "A.WY1", "A.WY2", "A.WY3","A.CO4")
#brewer.pal(11,"RdYlBu")
#col=heat.colors(11)

## Look at mean and var of CV across blocks for each source
#ARFR.ht.cv <- NULL

#for (vv in 1:10) {
#  dat.b <- subset(ARFR, ARFR$Block==vv)
# ARFR.ht.cv <- rbind(ARFR.ht.cv, dat.b %>% group_by(Source) %>%
#                summarise(Height_CV=cv(Length_cm_20220726, na.rm=TRUE)))
#}

## Test for significance in CV differences
#ARFR.ht.CV.aov <- aov(Height_CV ~ Source, data=ARFR.ht.cv)
#summary(ARFR.ht.CV.aov)

## Boxplot
#boxplot(Height_CV ~ Source, data=ARFR.ht.cv,
#        xlab="Population", ylab="CV in plant height", cex.lab=1.5,
#        names=c("A.AZ1", "A.AZ2", "A.CO1", "A.CO2","A.CO3","A.NM1","A.UT1", 
#                "A.WY1", "A.WY2", "A.WY3","A.CO4"), cex.axis=1.25,
#        main="Artemisia frigida", cex.main=1.5, col="darkseagreen3")

#ARFR.ht.mn.se <- ARFR.ht.cv %>% group_by(Source) %>%
#                 summarise(MN=mean(Height_CV, na.rm=TRUE), SE=calcSE(Height_CV))

#plotCI(barplot(ARFR.ht.mn.se$MN, col="darkseagreen3", ylab="CV in plant height", xlab="Population",
#               names=c("A.AZ1", "A.AZ2", "A.CO1", "A.CO2","A.CO3","A.NM1","A.UT1", 
#                       "A.WY1", "A.WY2", "A.WY3","A.CO4"), main="Artemisia frigida", cex.main=1.5,
#               cex.names=0.7, ylim=c(0,0.6), cex.lab=1.4),
#               ARFR.ht.mn.se$MN, uiw=ARFR.ht.mn.se$SE, add=TRUE, pch=NA, sfrac=0)


#ARFR$Replaced_YorN_20220531[ARFR$Replaced_YorN_20220531==""] <- "N" 
#If died after planting and not replaced, use subsequent surv data.
#ARFR.cl[!is.na(ARFR.cl$Replaced_YorN_20220531),]