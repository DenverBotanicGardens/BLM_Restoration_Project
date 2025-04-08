## April Goebl
## Script started 2023-01-29
## BLM Restoration project at Denver Botanic Gardens
## Analyze data from Chatfield Common Garden  




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
#setwd("C:/Users/april/Denver Botanic Gardens/Conservation - BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
#PEVI <- read.csv(file="Chatfield/20230303_ChatfieldData_PEVI.csv", sep=",", header=TRUE, dec=".")
ARFR <- read.csv(file="Chatfield/20230616_ChatfieldData2022_ARFR.csv", sep=",", header=TRUE, dec=".", na.strings="")
#ERNA <- read.csv(file="Chatfield/20230302_ChatfieldData_ERNA.csv", sep=",", header=TRUE, dec=".")
#BOGR <- read.csv(file="Chatfield/20230301_ChatfieldData_BOGR.csv", sep=",", header=TRUE, dec=".")

ARFR.SdZn <- read.csv(file="AGoebl/Seeds/20230823_ARFR_LatLong.csv", sep=",", header=TRUE, dec=".")
ARFR.ppt <- read.csv(file="AGoebl/Seeds/20230311_ARFR_pptAnnual.csv", sep=",", header=TRUE, dec=".")
ARFR.tmin <- read.csv(file="AGoebl/Seeds/20230311_ARFR_tminWinter.csv", sep=",", header=TRUE, dec=".")
ARFR.biovar <- readRDS("AGoebl/Seeds/20230814_ARFR_BiovarsAvg1980_2021")

#BOGR.SdZn <- read.csv(file="AGoebl/Seeds/20230824_BOGR_LatLongSdZn.csv", sep=",", header=TRUE, dec=".")
#BOGR.biovar <- readRDS("AGoebl/Seeds/20230825_BOGR_BiovarsAvg1980_2021")

#ERNA.SdZn <- read.csv(file="AGoebl/Seeds/Join_Features_to_SeedZones/ERNA_latLong_SdZone.csv", sep=",", header=TRUE, dec=".")
## ------------------------------------------------------------------------------------------------




## ARFR - DATA CLEAN UP ---------------------------------------------
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


## Checks 
#Some plts were dead that were selected to be harvested, check? 
ARFR.Coll <- ARFR.cl[!is.na(ARFR.cl$Harvest_20221014) | !is.na(ARFR.cl$Harvest_20221110),]
#Current data sheet only have "Harvest" marked for plts there were alive
ARFR.Coll[is.na(ARFR.Coll$AGB_MinusBag),] #Two harvested plts do not have final AGB: 1107 and 1274
ARFR.MissinfBM <- ARFR.Coll[is.na(ARFR.Coll$InfBM_Wbag),]
ARFR.MissinfBM$ID[ARFR.MissinfBM$Phenology_20220922==3] #Looks for these plts and get inf weights

## ** Add checks listed in BOGR ** 
## ----------------------------------------------------------------------------------------------




## ARFR - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ARFR.SdZn$Source <- str_replace(ARFR.SdZn$SOURCE_CODE, "4-SOS", "")
ARFR.ppt$Source <- str_replace(ARFR.ppt$Source, "4-SOS", "")
ARFR.tmin$Source <- str_replace(ARFR.tmin$Code, "4-SOS", "")
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

## Add colour columns that corresponds to seed zone and pop
#SdZn.list <- unique(ARFR.SdZn$SdZone)
ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "powderblue"        #semi-humid, cold
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen"     #semi-humid, cool
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 3 - 6", ARFR.SdZn$SdZone)] = "darkseagreen4"    #semi-humid, warm
ARFR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "darkgoldenrod1"   #semi-arid, cold
ARFR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "orange3"         #semi-arid, cool
ARFR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 6 - 12", ARFR.SdZn$SdZone)] = "tomato2"         #semi-arid, warm

## 'Order' by latitude 
pop.list <- unique(ARFR.SdZn$Source)
ARFR.SdZn$PopCol[grepl("ARFR-AZ930-423-NAVAJO-18", ARFR.SdZn$Source)] = "red4"        
ARFR.SdZn$PopCol[grepl("ARFR-AZ930-422-NAVAJO-18", ARFR.SdZn$Source)] = "red2"   
ARFR.SdZn$PopCol[grepl("ARFR-NM930N-66-11", ARFR.SdZn$Source)] = "orangered3"    
ARFR.SdZn$PopCol[grepl("ARFR-WY930-44-LASANIMAS-13", ARFR.SdZn$Source)] = "orange2"   
ARFR.SdZn$PopCol[grepl("ARFR-CO932-314-JEFFERSON-12", ARFR.SdZn$Source)] = "darkgoldenrod1"        
ARFR.SdZn$PopCol[grepl("ARFR-CO932-316-JEFFERSON-12", ARFR.SdZn$Source)] = "lightgoldenrod"        
ARFR.SdZn$PopCol[grepl("ARFR-UT080-109-UINTAH-12", ARFR.SdZn$Source)] = "darkseagreen"        
ARFR.SdZn$PopCol[grepl("ARFR-CO932-294-11", ARFR.SdZn$Source)] = "darkseagreen4"     
ARFR.SdZn$PopCol[grepl("ARFR-WY040-71-10", ARFR.SdZn$Source)] = "skyblue"    
ARFR.SdZn$PopCol[grepl("ARFR-WY050-151-FREMONT-16", ARFR.SdZn$Source)] = "steelblue"   
ARFR.SdZn$PopCol[grepl("ARFR-WY050-49-FREMONT-12", ARFR.SdZn$Source)] = "steelblue4"
## ----------------------------------------------------------------------------------------------

 
 

## ARFR - COMBINE DATA TYPES --------------------------------------------
ARFR.cl <- left_join(ARFR.cl, ARFR.SdZn, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.ppt, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.tmin, by="Source")
ARFR.cl <- left_join(ARFR.cl, ARFR.biovar, by="Source")
## ---------------------------------------------------------------




## ARFR - CONSOLIDATE PHENOLOGY DATA ------------------------------------
#ARFR.Pheno <- ARFR %>% dplyr::select(c(Source, Treatment, Block,SdZone_approx,PopAbbrev, Lat, starts_with("Phenology")))
ARFR.PlantingDate <- as.Date("2022-05-03")
ARFR.PhenoCol.List <- colnames(ARFR)[grepl("Phenology*", colnames(ARFR))] #Obtain phenology column names
ARFR.Pheno.List <- str_replace(ARFR.PhenoCol.List, "Phenology_", "")      #Obtain just date from phenology columns
ARFR.Pheno.List <- as.Date(ARFR.Pheno.List, "%Y%m%d")
ARFR.DaysToFlwr <- ARFR.Pheno.List - ARFR.PlantingDate                    #Calculate days from planting to each phenology survey 

## Loop over each phenology column & enter the num days since planting when a 3 (flower or later phenophase) first appears
ARFR.cl$DaysToFlwr <- NA
for (pp in 1:length(ARFR.PhenoCol.List)) {
  ARFR.cl$DaysToFlwr[ARFR.cl[,ARFR.PhenoCol.List[pp]]>=3 & is.na(ARFR.cl$DaysToFlwr)] <- as.integer(ARFR.DaysToFlwr)[pp]
}
## ---------------------------------------------------------------



## ARFR - ADD GROWTH RATE VARIABLES
ARFR.cl$GrwthRate_Specific <- log(ARFR.cl$Length_cm_20220726/ARFR.cl$Length_cm_20220622)
ARFR.cl$GrwthRate_Absolute <- ARFR.cl$Length_cm_20220726-ARFR.cl$Length_cm_20220622
ARFR.cl$GrwthRate_Relative <- (ARFR.cl$Length_cm_20220726-ARFR.cl$Length_cm_20220622)/ARFR.cl$Length_cm_20220622

## ** Look at early vs late growth? ** 
## ---------------------------------------------------------------






## ARFR - TEST FOR TREATMENT EFFECT -------------------------------------
hist(log(ARFR$AGB_MinusBag))
hist(ARFR$AGB_MinusBag)
#ARFR.tx.mod <- aov(ARFR$AGB_MinusBag ~ ARFR$Source + ARFR$Treatment)
ARFR.tx.mod <- lmer(log(AGB_MinusBag) ~ Source + Treatment + (1|Block), data=ARFR)
summary(ARFR.tx.mod)
ARFR.pop.mod <- lmer(log(AGB_MinusBag) ~ Source + (1|Block), data=ARFR)

models <- list(ARFR.tx.mod, ARFR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment 

hist(log(ARFR$Length_cm_20220726))
hist(ARFR$Length_cm_20220726)
#ARFR.tx.mod <- aov(log(ARFR$Length_cm_20220726) ~ ARFR$Source + ARFR$Treatment)
ARFR.tx.modlog <- lmer(log(ARFR$Length_cm_20220726) ~ Source + Treatment + (1|Block), data=ARFR)
ARFR.pop.modlog <- lmer(log(ARFR$Length_cm_20220726) ~ Source + (1|Block), data=ARFR)
ARFR.tx.mod <- lmer(ARFR$Length_cm_20220726 ~ Source + Treatment + (1|Block), data=ARFR)
ARFR.pop.mod <- lmer(ARFR$Length_cm_20220726 ~ Source + (1|Block), data=ARFR)

models <- list(ARFR.tx.mod, ARFR.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment? 
## -----------------------------------------------------------------------------------------




## ARFR - VISUALIZE RAW DATA ---------------------------------------------------------------
## Order populations for plotting 

## EITHER --- Order by seed zone -----
## Order based on aridity & temp (cool & wet to hot & dry), w/in each cat, order from N to S
ARFR.SdZn$Source <- factor(ARFR.SdZn$Source, levels=c("ARFR-WY040-71-10","ARFR-UT080-109-UINTAH-12",
                                                      "ARFR-CO932-314-JEFFERSON-12","ARFR-NM930N-66-11",
                                                      "ARFR-CO932-294-11","ARFR-CO932-316-JEFFERSON-12",
                                                      "ARFR-WY930-44-LASANIMAS-13","ARFR-WY050-49-FREMONT-12",
                                                      "ARFR-WY050-151-FREMONT-16","ARFR-AZ930-422-NAVAJO-18",
                                                      "ARFR-AZ930-423-NAVAJO-18"))
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Source),]
ARFR.cl$Source <- factor(ARFR.cl$Source, levels=ARFR.SdZn$Source)

## OR --- Order by latitude ----
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Lat),]
ARFR.cl$Source <- factor(ARFR.cl$Source, levels=ARFR.SdZn$Source)

## OR --- Order by average trait value ---
## **** ## add other traits ** 
ARFR.htByMed <- with(ARFR.cl, reorder(Source, Length_cm_20220726, median, na.rm=TRUE))
ARFR.htMeds <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(HT_LATE=median(Length_cm_20220726,na.rm=TRUE),
                                                AGB=median(AGB_MinusBag,na.rm=TRUE), REPRO_BM=median(InfBM_Wbag,na.rm=TRUE))
ARFR.SdZn <- left_join(ARFR.SdZn, ARFR.htMeds, by="Source")
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$HT_LATE),]



## Boxplots of raw data 
boxplot(Length_cm_20220726 ~ Source, data=ARFR.cl,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        cex.axis=1.1, names=ARFR.SdZn$PopAbbrev,
        main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$PopCol)

boxplot(Length_cm_20220726 ~ ARFR.htByMed, data=ARFR.cl,
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        cex.axis=1.1, names=ARFR.SdZn$PopAbbrev,
        main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

## Add time interval to growth rate calcs? **
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

boxplot(DaysToFlwr ~ Source, data=ARFR,
        xlab="Population", ylab="Days to first flower", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

boxplot(InfBM_Wbag ~ Source, data=ARFR.cl,
        xlab="Population", ylab="Reproductive biomass (g)", cex.lab=1.5, names=ARFR.SdZn$PopAbbrev,
        cex.axis=1.1, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## ---------------------------------------   




## ARFR - ESTIMATE VARIATION B/W POPULATIONS -------------------------------------------------------------
## Use Emmeans
## Final size (height and biomass)
ARFR.htMod <- lmer(log(Length_cm_20220726) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.ht <- emmeans(ARFR.htMod, specs = pairwise ~ Source)
summary(ARFR.pair.ht)
plot(ARFR.pair.ht)

ARFR.bmMod <- lmer(log(AGB_MinusBag) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.bm <- emmeans(ARFR.bmMod, specs = pairwise ~ Source)
summary(ARFR.pair.bm)
plot(ARFR.pair.bm)

## Plant growth (specific GR)
ARFR.grMod <- lmer(log(Length_cm_20220726/Length_cm_20220622) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.gr <- emmeans(ARFR.grMod, specs = pairwise ~ Source)
summary(ARFR.pair.gr)
plot(ARFR.pair.gr)

## Repro BM
ARFR.rbmMod <- lmer(log(InfBM_Wbag) ~ Source + (1|Block), data = ARFR.cl)
ARFR.pair.rbm <- emmeans(ARFR.rbmMod, specs = pairwise ~ Source)
summary(ARFR.pair.bm)
plot(ARFR.pair.rbm)

## Models show support for Source as a predictor variable
Anova(ARFR.htMod)
Anova(ARFR.grMod)
Anova(ARFR.bmMod)
Anova(ARFR.rbmMod)


## Does seed zone or lat explain variation?
ARFR.ht.szLat <- lmer(Length_cm_20220726 ~ as.factor(seed_zone) + Lat + (1|Block), data=ARFR.cl)
Anova(ARFR.ht.szLat)
summary(ARFR.ht)

ARFR.ht.sz <- lmer(Length_cm_20220726 ~ as.factor(seed_zone) + (1|Block), data=ARFR.cl)
ARFR.ht.lat <- lmer(Length_cm_20220726 ~ Lat + (1|Block), data=ARFR.cl)
models <- list(ARFR.ht.szLat, ARFR.ht.sz, ARFR.ht.lat)
mod.names <- c('Both', 'JustSz', 'JustLat')
aictab(cand.set = models, modnames = mod.names )


## Use PCA to look at clustering of individuals and if there is grouping by population
## Try PCA using different (not cov) approach
#ARFR.traits <- ARFR.cl %>% dplyr::select(c("Length_cm_20220726","Survival_20220922","DaysToFlwr",
#                                           "AGB_MinusBag","InfBM_Wbag","GrwthRate_Absolute"))
ARFR.traits <- ARFR.cl %>% dplyr::select(c("PopCol","Source","Length_cm_20220726","DaysToFlwr",
                                           "GrwthRate_Specific"))
ARFR.traitsCen <- scale(ARFR.traits[,3:ncol(ARFR.traits)], center=TRUE, scale=TRUE)
ARFR.traitsComb <- cbind(ARFR.traits[,1:2],ARFR.traitsCen)
ARFR.traitsComb.noNA <- na.omit(ARFR.traitsComb)
pca.pops <- prcomp(ARFR.traitsComb.noNA[,3:ncol(ARFR.traitsComb.noNA)]) #b/c of NAs there are few indivs w/ data for all traits. Figure out how to deal with NAs
biplot(pca.pops)
plot(x=pca.pops$x[,1], y=pca.pops$x[,2],col=ARFR.traitsComb.noNA$PopCol, pch=19, cex=0.5)
plot(x=pca.pops$x[,1], y=pca.pops$x[,3],col=ARFR.traitsComb.noNA$PopCol, pch=19, cex=0.63,ylim=c(-2.2,2.2))
plot(x=pca.pops$x[,2], y=pca.pops$x[,3],col=ARFR.traitsComb.noNA$PopCol, pch=19, cex=0.63,ylim=c(-2.2,2.3))

# ** remove outlier point ***

## Look at scree plot
pca.pops$sdev[1]**2/sum(pca.pops$sdev**2)
var.expl <- pca.pops$sdev^2/sum(pca.pops$sdev^2)
barplot(pca.pops$sdev[1:11]**2/sum(pca.pops$sdev**2))
## ---------------------------------------------------------------------------------------------------




## ARFR - ESTIMATE VARIATION WITHIN POPULATIONS -------------------------------------------------------------
## Order populations for plotting 
## EITHER --- Order by seed zone -----
## Order based on aridity & temp (cool & wet - warm & wet - cool & dry - hot & dry), w/in each cat, order from N to S
ARFR.SdZn$Source <- factor(ARFR.SdZn$Source, levels=c("ARFR-WY040-71-10","ARFR-UT080-109-UINTAH-12",
                                                      "ARFR-CO932-314-JEFFERSON-12","ARFR-NM930N-66-11",
                                                      "ARFR-CO932-294-11","ARFR-CO932-316-JEFFERSON-12",
                                                      "ARFR-WY930-44-LASANIMAS-13","ARFR-WY050-49-FREMONT-12",
                                                      "ARFR-WY050-151-FREMONT-16","ARFR-AZ930-422-NAVAJO-18",
                                                      "ARFR-AZ930-423-NAVAJO-18"))
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Source),]
ARFR.cl$Source <- factor(ARFR.cl$Source, levels=ARFR.SdZn$Source)

## OR --- Order by latitude ----
ARFR.SdZn <- ARFR.SdZn[order(ARFR.SdZn$Lat),]
ARFR.cl$Source <- factor(ARFR.cl$Source, levels=ARFR.SdZn$Source)


## Calculate the coefficient of variation 
ARFR.ht.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Height_CV=cv(Length_cm_20220726, na.rm=TRUE))
ARFR.gr.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Growth_CV=cv(log(Length_cm_20220726/Length_cm_20220622), na.rm=TRUE))
ARFR.bm.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Biomass_CV=cv(AGB_MinusBag, na.rm=TRUE))
ARFR.dtf.cv <- ARFR.cl %>% group_by(Source) %>% summarise(DaysToFlwr_CV=cv(DaysToFlwr, na.rm=TRUE))
ARFR.rbm.cv <- ARFR.cl %>% group_by(Source) %>% summarise(Rbiomass_CV=cv(InfBM_Wbag, na.rm=TRUE))
## ---------


## Barplots
barplot(ARFR.ht.cv$Height_CV, xlab="Population", ylab="CV in plant height", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.gr.cv$Growth_CV, xlab="Population", ylab="CV in plant growth", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.bm.cv$Biomass_CV, xlab="Population", ylab="CV in above-ground biomass", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.dtf.cv$DaysToFlwr_CV, xlab="Population", ylab="CV in days to first flower", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)

barplot(ARFR.rbm.cv$Rbiomass_CV, xlab="Population", ylab="CV in reproductive biomass", cex.lab=1.5, cex.axis=1.1, 
        names=ARFR.SdZn$PopAbbrev, main="Artemisia frigida", cex.main=1.5, col=ARFR.SdZn$SdZnCol)
## ------------------------------------------------------------------------------




## ARFR - LOOK AT RELATIONSHIP BETWEEN TRAIT VAR AND AVERAGE ---------------------------
ARFR.ht.mn <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Height_MN=mean(Length_cm_20220726, na.rm=TRUE))
ARFR.gr.mn <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Growth_MN=mean(log(Length_cm_20220726/Length_cm_20220622), na.rm=TRUE))
ARFR.bm.mn <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Biomass_MN=mean(AGB_MinusBag, na.rm=TRUE))
ARFR.rbm.mn <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Rbiomass_MN=mean(InfBM_Wbag, na.rm=TRUE))

ARFR.ht.sd <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Height_SD=sd(Length_cm_20220726, na.rm=TRUE))
ARFR.gr.sd <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Growth_SD=sd(log(Length_cm_20220726/Length_cm_20220622), na.rm=TRUE))
ARFR.bm.sd <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Biomass_SD=sd(AGB_MinusBag, na.rm=TRUE))
ARFR.rbm.sd <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(Rbiomass_SD=sd(InfBM_Wbag, na.rm=TRUE))

#plot(ARFR.ht.cv$Height_CV, ARFR.ht.mn$Height_MN)
#plot(ARFR.gr.cv$Growth_CV, ARFR.gr.mn$Growth_MN)
#plot(ARFR.bm.cv$Biomass_CV, ARFR.bm.mn$Biomass_MN)
#plot(ARFR.rbm.cv$Rbiomass_CV, ARFR.rbm.mn$Rbiomass_MN)

plot(ARFR.ht.sd$Height_SD, ARFR.ht.mn$Height_MN, pch=19, col=ARFR.SdZn$PopCol, cex=1.3)
plot(ARFR.gr.sd$Growth_SD, ARFR.gr.mn$Growth_MN, pch=19, col=ARFR.SdZn$PopCol, cex=1.25)
plot(ARFR.bm.sd$Biomass_SD, ARFR.bm.mn$Biomass_MN, pch=19)
plot(ARFR.rbm.sd$Rbiomass_SD, ARFR.rbm.mn$Rbiomass_MN, pch=19)
## ----------------------------------------------------------




## ARFR - EVALUATE RELATIONSHIPS B/W TRAITS AND SOURCE CLIMATE -------------------------------------------
## Visualize raw data
ARFR.summ <- ARFR.cl %>% group_by(Source) %>% 
  dplyr::summarise(AGB_MEAN=mean(AGB_MinusBag, na.rm=TRUE), AGB_SE=calcSE(AGB_MinusBag),
            HT_MEAN=mean(Length_cm_20220726, na.rm=TRUE), HT_SE=calcSE(Length_cm_20220726),
            INFBM_MEAN=mean(InfBM_Wbag, na.rm=TRUE), INFBM_SE=calcSE(InfBM_Wbag), BIO1=mean(bio1, na.rm=TRUE), 
            BIO2=mean(bio2, na.rm=TRUE), BIO4=mean(bio4, na.rm=TRUE), BIO7=mean(bio7, na.rm=TRUE), BIO12=mean(bio12, na.rm=TRUE),
            BIO9=mean(bio9, na.rm=TRUE), BIO14=mean(bio14, na.rm=TRUE), BIO15=mean(bio15, na.rm=TRUE), BIO18=mean(bio18, na.rm=TRUE),
            PPT=mean(Mean_Annual_Ppt, na.rm=TRUE), TMIN=mean(Mean_MinWinter_Temp, na.rm=TRUE))

ARFR.summ <- left_join(ARFR.summ, ARFR.SdZn, by="Source")

## Look at correlations in climate, latitude, etc.
plot(ARFR.summ$TMIN, ARFR.summ$Lat, pch=19)
summary(lm(ARFR.summ$TMIN ~ ARFR.summ$Lat))


plot(ARFR.summ$PPT, ARFR.summ$AGB_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida", 
     cex.main=1.5, xlab="Mean annual precipitation", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$PPT, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$PPT, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$TMIN, ARFR.summ$AGB_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida", 
     cex.main=1.5, xlab="Mean winter temperature", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$TMIN, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$TMIN, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$Lat, ARFR.summ$AGB_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida",
     cex.main=1.5, xlab="Latitude", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$Lat, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$Lat, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE,
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.5)

plot(ARFR.summ$Lat, ARFR.summ$HT_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.25, main="Artemisia frigida",
     cex.main=1.5, xlab="Latitude", ylab="Plant Height (cm)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$Lat, ARFR.summ$HT_MEAN-ARFR.summ$HT_SE, ARFR.summ$Lat, ARFR.summ$HT_MEAN+ARFR.summ$HT_SE,
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.5)


## Models
ARFR.cl$Ppt.sq <- (ARFR.cl$Mean_Annual_Ppt)^2
ARFR.bm.mod <- lmer(log(AGB_MinusBag) ~ Mean_Annual_Ppt + Ppt.sq + Mean_MinWinter_Temp + (1|Block), data=ARFR.cl)
summary(ARFR.bm.mod)
Anova(ARFR.bm.mod)

plot(allEffects(ARFR.bm.mod))
plot(predictorEffects(ARFR.bm.mod))
#Min winter temp looks important 



## 19 Bioclimate variables ------------------------------------------------------
## Look at correlations b/w climate variables to decide which to include in models 
cor.mat <- cor(ARFR.biovar[,2:20], method="pearson")
corrplot(cor.mat)
chart.Correlation(cor.mat, histogram=TRUE, method="pearson")

cor.p.mats <- rcorr(as.matrix(ARFR.biovar[,2:20]), type="pearson")

##Function to re-format the output (from http://www.sthda.com/english/wiki/...)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(row=rownames(cormat)[row(cormat)[ut]],
             column=rownames(cormat)[col(cormat)[ut]],
             cor=(cormat)[ut], p=pmat[ut]) }

cor.p.tbl <- flattenCorrMatrix(cor.p.mats$r, cor.p.mats$P)
cor.80 <- cor.p.tbl[cor.p.tbl$cor>0.8 | cor.p.tbl$cor<(-0.8),] #Filter by correlations > 80%
cor.80 <- cor.80[order(cor.80$row, cor.80$column),]
## List of variables to include in model if excluding vars that are > 80% correlated
#bio1, bio2, bio4, bio7, bio9, bio12, bio14, bio15, bio16

## Models
ARFR.ht.bv.mod <- lmer(Length_cm_20220726 ~ scale(bio1) + scale(bio2) + scale(bio4) + scale(bio7)
                       + scale(bio9) + scale(bio12) + scale(bio14) + scale(bio15) 
                       + (1|Block), data=ARFR.cl)
summary(ARFR.ht.bv.mod)
Anova(ARFR.ht.bv.mod)

plot(allEffects(ARFR.ht.bv.mod))
plot(predictorEffects(ARFR.ht.bv.mod))
#Bio9 looks like strongest relationship, Bio16 is not significant. 


## Look at PCA of 19 variables to possibly reduce the number of predictors 
#ARFR.biovarCen <- scale(ARFR.biovar[,2:20], center=TRUE, scale=TRUE)
#ARFR.biovarComb <- as.data.frame(cbind(ARFR.biovar[,21],ARFR.biovarCen))
pca.biovar <- prcomp(ARFR.biovar[,2:20]) 
plot(x=pca.biovar$x[,1], y=pca.biovar$x[,2], cex=1.5)
plot(x=pca.biovar$x[,1], y=pca.biovar$x[,3], pch=19)
plot(x=pca.biovar$x[,2], y=pca.biovar$x[,3], pch=19, cex=0.63)

## Add arrows on PCA plot and look at loadings 
biplot(pca.biovar)
pca.biovar$rotation

## Look at scree plot
pca.biovar$sdev[1]**2/sum(pca.biovar$sdev**2)
var.expl <- pca.biovar$sdev^2/sum(pca.biovar$sdev^2)
barplot(pca.results$sdev[1:19]**2/sum(pca.results$sdev**2))
## -----------


## Visualize raw data relationships with bioclimate vars
plot(ARFR.cl$Length_cm_20220726 ~ ARFR.cl$bio1)
plot(ARFR.cl$InfBM_Wbag ~ ARFR.cl$bio1)
plot(ARFR.cl$Length_cm_20220726 ~ ARFR.cl$bio9)
plot(ARFR.cl$InfBM_Wbag ~ ARFR.cl$bio9)
plot(ARFR.cl$Length_cm_20220726 ~ ARFR.cl$bio12)
plot(ARFR.cl$Length_cm_20220726 ~ ARFR.cl$bio4)

plot(ARFR.summ$BIO1, ARFR.summ$AGB_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.2, main="Artemisia frigida", 
     cex.main=1.5, xlab="BIO1", ylab="Above-ground biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$BIO1, ARFR.summ$AGB_MEAN-ARFR.summ$AGB_SE, ARFR.summ$BIO1, ARFR.summ$AGB_MEAN+ARFR.summ$AGB_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$BIO2, ARFR.summ$INFBM_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.2, main="Artemisia frigida", 
     cex.main=1.5, xlab="BIO2", ylab="Inf biomass (g)", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$BIO2, ARFR.summ$INFBM_MEAN-ARFR.summ$INFBM_SE, ARFR.summ$BIO2, ARFR.summ$INFBM_MEAN+ARFR.summ$INFBM_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$BIO9, ARFR.summ$HT_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.2, main="Artemisia frigida", 
     cex.main=1.5, xlab="BIO9", ylab="Plant height", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$BIO9, ARFR.summ$HT_MEAN-ARFR.summ$HT_SE, ARFR.summ$BIO9, ARFR.summ$HT_MEAN+ARFR.summ$HT_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)

plot(ARFR.summ$BIO12, ARFR.summ$HT_MEAN, col=ARFR.summ$SdZnCol, pch=19, cex=1.2, main="Artemisia frigida", 
     cex.main=1.5, xlab="BIO12", ylab="Plant height", cex.lab=1.5, cex.axis=1.1)
arrows(ARFR.summ$BIO12, ARFR.summ$HT_MEAN-ARFR.summ$HT_SE, ARFR.summ$BIO12, ARFR.summ$HT_MEAN+ARFR.summ$HT_SE, 
       angle=90, col=ARFR.summ$SdZnCol, code=3, length=0, lwd=1.6)
## -------------------------------------------------------------------------------




## ARFR - EMMEANS ** is this diff from above? **------------------------------------------
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
## -------------------------------------------------------------------
## -------------------------------------------------------------------







## ----------------------------------------------------------------------------------------------
## BOGR - LOAD DATA ------------------------------------------------------------------------------------
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

## Add colour columns that corresponds to seed zone and pop
unique(BOGR.SdZn$seed_zone)
BOGR.SdZn$SdZnColful[grepl("10 - 15 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen2"    #semi-humid, cool
BOGR.SdZn$SdZnColful[grepl("15 - 20 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen"     #semi-humid, warm
BOGR.SdZn$SdZnColful[grepl("20 - 25 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = "darkseagreen4"    #semi-humid, v.warm
BOGR.SdZn$SdZnColful[grepl("10 - 15 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "lightgoldenrod1" #semi-arid, cool
BOGR.SdZn$SdZnColful[grepl("15 - 20 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "goldenrod2"      #semi-arid, warm
BOGR.SdZn$SdZnColful[grepl("20 - 25 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orange3"         #semi-arid, v.warm
BOGR.SdZn$SdZnColful[grepl("25 - 30 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orangered"       #semi-arid, hot
BOGR.SdZn$SdZnColful[grepl("30 - 35 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = "orangered3"      #semi-arid, v.hot
BOGR.SdZn$SdZnColful[grepl("5 - 10 Deg. F. / 12 - 30", BOGR.SdZn$seed_zone)] = "pink2"           #arid, cold

colfunc <- colorRampPalette(c("grey15","grey90"))
BOGR.SdZnCols <- colfunc(length(unique(BOGR.SdZn$seed_zone)))
BOGR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[1]    #semi-humid, cool
BOGR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[2]    #semi-humid, warm
BOGR.SdZn$SdZnCol[grepl("20 - 25 Deg. F. / 3 - 6", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[3]    #semi-humid, v.warm
BOGR.SdZn$SdZnCol[grepl("10 - 15 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[4]   #semi-arid, cool
BOGR.SdZn$SdZnCol[grepl("15 - 20 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[5]   #semi-arid, warm
BOGR.SdZn$SdZnCol[grepl("20 - 25 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[6]   #semi-arid, v.warm
BOGR.SdZn$SdZnCol[grepl("25 - 30 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[7]   #semi-arid, hot
BOGR.SdZn$SdZnCol[grepl("30 - 35 Deg. F. / 6 - 12", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[8]   #semi-arid, v.hot
BOGR.SdZn$SdZnCol[grepl("5 - 10 Deg. F. / 12 - 30", BOGR.SdZn$seed_zone)] = BOGR.SdZnCols[9]   #arid, cold

## 'Order' by latitude 
unique(BOGR.SdZn$Source)
BOGR.popCols <- colfunc(length(unique(BOGR.SdZn$Source)))
BOGR.SdZn$PopCol[grepl("BOGR-AZ040-302-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[1]   
BOGR.SdZn$PopCol[grepl("BOGR-AZ040-284-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[2] 
BOGR.SdZn$PopCol[grepl("BOGR-AZ040-278-SANTACRUZ-17", BOGR.SdZn$Source)] = BOGR.popCols[3]   
BOGR.SdZn$PopCol[grepl("BOGR-AZ930-302-10", BOGR.SdZn$Source)] = BOGR.popCols[8] 
BOGR.SdZn$PopCol[grepl("BOGR-AZ930-303-10", BOGR.SdZn$Source)] = BOGR.popCols[9]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-360-FREMONT-17", BOGR.SdZn$Source)] = BOGR.popCols[12]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-362-CHAFFEE-17", BOGR.SdZn$Source)] = BOGR.popCols[13]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-307-JEFFERSON-12", BOGR.SdZn$Source)] = BOGR.popCols[14]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-274-11", BOGR.SdZn$Source)] = BOGR.popCols[15]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-363-CLEARCREEK-17", BOGR.SdZn$Source)] = BOGR.popCols[16]   
BOGR.SdZn$PopCol[grepl("BOGR-CO932-276-11", BOGR.SdZn$Source)] = BOGR.popCols[18] 
BOGR.SdZn$PopCol[grepl("BOGR-NM-930-071-08", BOGR.SdZn$Source)] = BOGR.popCols[4]   
BOGR.SdZn$PopCol[grepl("BOGR-NM930-145-10", BOGR.SdZn$Source)] = BOGR.popCols[5]   
BOGR.SdZn$PopCol[grepl("BOGR-NM080-67-CHAVES-16", BOGR.SdZn$Source)] = BOGR.popCols[6]   
BOGR.SdZn$PopCol[grepl("BOGR-NM930-140-10", BOGR.SdZn$Source)] = BOGR.popCols[7]   
BOGR.SdZn$PopCol[grepl("BOGR-NM930N-95-SANDOVAL-12", BOGR.SdZn$Source)] = BOGR.popCols[10]   
BOGR.SdZn$PopCol[grepl("BOGR-UT030-215-GARFIELD-13", BOGR.SdZn$Source)] = BOGR.popCols[11]    
BOGR.SdZn$PopCol[grepl("BOGR-UT080-163-CARBON-14", BOGR.SdZn$Source)] = BOGR.popCols[17]    
BOGR.SdZn$PopCol[grepl("BOGR-WY050-133-FREMONT-16", BOGR.SdZn$Source)] = BOGR.popCols[19]       
BOGR.SdZn$PopCol[grepl("BOGR-WY050-113-FREMONT-15", BOGR.SdZn$Source)] = BOGR.popCols[20] 
BOGR.SdZn$PopCol[grepl("BOGR-WY070-71-JOHNSON-15", BOGR.SdZn$Source)] = BOGR.popCols[21]  


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

## OR --- Order by latitude ----

## OR --- Order by average trait value ---
BOGR.htByMed <- with(BOGR.cl, reorder(Source, Length_cm_20220801, median, na.rm=TRUE))
BOGR.infByMed <- with(BOGR.cl, reorder(Source, NumInf, median, na.rm=TRUE))
#BOGR.graByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Absolute, median, na.rm=TRUE))
BOGR.grrByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Relative, median, na.rm=TRUE))
#BOGR.grsByMed <- with(BOGR.cl, reorder(Source, GrwthRate_Specific, median, na.rm=TRUE))
BOGR.dfByMed <- with(BOGR.cl, reorder(Source, DaysToFlwr, median, na.rm=TRUE))


## Calculate median trait values for each source pop
BOGR.meds <- BOGR.cl %>% group_by(Source) %>% 
  dplyr::summarise(Height_MD=median(Length_cm_20220801,na.rm=TRUE),
                   Inf_MD=median(NumInf,na.rm=TRUE),Growth_MD=median(GrwthRate_Relative,na.rm=TRUE),
                   DaysToFlwr_MD=median(DaysToFlwr,na.rm=TRUE))


## OR --- Order by average size ---
BOGR.meds <- BOGR.meds[order(BOGR.meds$Height_MD),]
BOGR.cl$Source <- factor(BOGR.cl$Source, levels=BOGR.meds$Source)




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
        xlab=NA, ylab="Plant height (cm)", cex.lab=1.25,
        cex.axis=0.9, names=BOGR.meds$PopAbbrev, las=2, col=BOGR.meds$SdZnColful)
        #main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)
legend("topleft", unique(BOGR.meds$SdZnAbbrev[order(BOGR.meds$SdZnOrder)]), 
       col=unique(BOGR.meds$SdZnColful[order(BOGR.meds$SdZnOrder)]), cex=1.95, pch=19)


## Growth rate(s) ** Add time interval to growth rate calcs? **
boxplot(GrwthRate_Relative ~ Source, data=BOGR.cl, las=2,
        xlab=NA, ylab="Plant relative growth", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.9, main=NA, col=BOGR.meds$SdZnColful)

#BOGR.meds <- BOGR.meds[order(BOGR.meds$Growth_MD),]
#boxplot(GrwthRate_Relative ~ BOGR.grrByMed, data=BOGR.cl, las=2,
#        xlab="Population", ylab="Plant relative growth", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)

#BOGR.SdZn <- BOGR.SdZn[order(BOGR.SdZn$GR_AB),]
#boxplot(GrwthRate_Absolute ~ BOGR.graByMed, data=BOGR.cl, las=2,
#        xlab="Population", ylab="Plant Absolute growth", cex.lab=1.25, names=BOGR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.SdZn$PopCol)

#BOGR.SdZn <- BOGR.SdZn[order(BOGR.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ BOGR.grsByMed, data=BOGR.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=BOGR.SdZn$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.SdZn$PopCol)

## Repro
boxplot(NumInf ~ Source, data=BOGR.cl, las=2,
        xlab=NA, ylab="Number of inflorescences", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.9, main=NA, col=BOGR.meds$SdZnColful)

#BOGR.meds <- BOGR.meds[order(BOGR.meds$Inf_MD),]
#boxplot(NumInf ~ BOGR.infByMed, data=BOGR.cl, las=2,
#        xlab="Population", ylab="Number of inflorescences", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)

## Phenology
boxplot(DaysToFlwr ~ Source, data=BOGR.cl, las=2,
        xlab=NA, ylab="Days to flower", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
        cex.axis=0.9, main=NA, col=BOGR.meds$SdZnColful)

#BOGR.meds <- BOGR.meds[order(BOGR.meds$DaysToFlwr_MD),]
#boxplot(DaysToFlwr ~ BOGR.dfByMed, data=BOGR.cl, las=2,
#        xlab="Population", ylab="Days to flower", cex.lab=1.25, names=BOGR.meds$PopAbbrev,
#        cex.axis=0.7, main="Bouteloua gracilis", cex.main=1.5, col=BOGR.meds$SdZnColful)
## ---------------------------------------------------





## BOGR - ESTIMATE VARIATION B/W POPULATIONS -------------------------------------------------------------
## Use Emmeans ----------------
## Final size (height)
BOGR.htMod <- lmer(Length_cm_20220801 ~ Source + (1|Block), data=BOGR.cl)
BOGR.pair.ht <- emmeans(BOGR.htMod, specs = pairwise ~ Source)
#summary(BOGR.pair.ht)
plot(BOGR.pair.ht)

BOGR.ht.pw <- as.data.frame(BOGR.pair.ht$contrasts)
BOGR.ht.pwDiffs <- BOGR.ht.pw$estimate
max(BOGR.ht.pwDiffs) 


## Plant growth
BOGR.grMod <- lmer(GrwthRate_Relative ~ Source + (1|Block), data=BOGR.cl)
BOGR.pair.gr <- emmeans(BOGR.grMod, specs = pairwise ~ Source)
#summary(BOGR.gr.em)
plot(BOGR.pair.gr)

BOGR.gr.means <- as.data.frame(BOGR.pair.gr$emmeans)
BOGR.gr.meanEsts <- BOGR.gr.means$emmean
max(BOGR.gr.meanEsts)
min(BOGR.gr.meanEsts)

BOGR.gr.pw <- as.data.frame(BOGR.pair.gr$contrasts)
BOGR.gr.pwDiffs <- exp(BOGR.gr.pw$estimate)
max(BOGR.gr.pwDiffs)
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

## Use PCA -----------------
## Look at clustering of individuals & if there is grouping by pop
BOGR.traits <- BOGR.cl %>% dplyr::select(c("SdZnColful","PopCol","Source","Length_cm_20220801","DaysToFlwr",
                                           "GrwthRate_Relative","NumInf"))
#BOGR.traitsCen <- scale(BOGR.traits[,4:ncol(BOGR.traits)], center=TRUE, scale=TRUE)
#BOGR.traitsComb <- cbind(BOGR.traits[,1:2], BOGR.traitsCen)
BOGR.traits.noNA <- na.omit(BOGR.traits)
#BOGR.pcaPops <- prcomp(BOGR.traitsComb.noNA[,3:ncol(BOGR.traitsComb.noNA)]) #b/c of NAs there are fewer indivs w/ data for all traits. Figure out how to deal with NAs
BOGR.pcaPops <- prcomp(BOGR.traits.noNA[4:ncol(BOGR.traits.noNA)], scale=TRUE) 
biplot(BOGR.pcaPops, pch=19)
plot(x=BOGR.pcaPops$x[,1], y=BOGR.pcaPops$x[,2],col=BOGR.traits.noNA$SdZnColful, pch=19, cex=0.8)
plot(x=BOGR.pcaPops$x[,2], y=BOGR.pcaPops$x[,3],col=BOGR.traits.noNA$SdZnColful, pch=19, cex=0.8)
plot(x=BOGR.pcaPops$x[,3], y=BOGR.pcaPops$x[,4],col=BOGR.traits.noNA$SdZnColful, pch=19, cex=0.8)

## Look at scree plot
BOGR.pcaPops$sdev[1]**2/sum(BOGR.pcaPops$sdev**2) #variation explained by PC1
BOGR.pcaPops$sdev^2/sum(BOGR.pcaPops$sdev^2)
barplot(BOGR.pcaPops$sdev[1:11]**2/sum(BOGR.pcaPops$sdev**2))

## Look at loadings
BOGR.pcaPops$rotation
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

plot(BOGR.ht.sd$Height_SD, BOGR.ht.mn$Height_MN, pch=19, col=BOGR.SdZn$PopCol, cex=1.3)
plot(BOGR.gr.sd$Growth_SD, BOGR.gr.mn$Growth_MN, pch=19, col=BOGR.SdZn$PopCol, cex=1.2)
plot(BOGR.df.sd$DaysToFlwr_SD, BOGR.df.mn$DaysToFlwr_MN, pch=19, col=BOGR.SdZn$PopCol, cex=1.2)
plot(BOGR.inf.sd$Inf_SD, BOGR.inf.mn$Inf_MN, pch=19,col=BOGR.SdZn$PopCol, cex=1.1)
## ----------------------------------------------------------





## BOGR - EVALUATE RELATIONSHIPS B/W TRAITS AND SOURCE CLIMATE -------------------------------------------
## 19 Bio-climate variables 

## Look at correlations b/w climate variables ------------------------------------------------------------
cor.mat <- cor(BOGR.biovar[,2:20], method="pearson")
corrplot(cor.mat)
#chart.Correlation(cor.mat, histogram=TRUE, method="pearson")

cor.p.mats <- rcorr(as.matrix(BOGR.biovar[,2:20]), type="pearson")
##Function to re-format the output (from http://www.sthda.com/english/wiki/...)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(row=rownames(cormat)[row(cormat)[ut]],
             column=rownames(cormat)[col(cormat)[ut]],
             cor=(cormat)[ut], p=pmat[ut]) }
 
cor.p.tbl <- flattenCorrMatrix(cor.p.mats$r, cor.p.mats$P)
cor.80 <- cor.p.tbl[cor.p.tbl$cor>0.8 | cor.p.tbl$cor<(-0.8),] #Filter by correlations > 80%
cor.80 <- cor.80[order(cor.80$row, cor.80$column),] #Use to decide which vars include in models?
## List of variables to include in model if excluding vars that are > 80% correlated
#bio1,bio2,bio3,bio5,bio8,bio12,bio13,bio14,bio15,bio17
#Exclude: bio4,bio6,bio7,bio9,bio10,bio11,bio16,bio18,bio19

## Models
#BOGR.ht.mod <- lmer(Length_cm_20220801 ~ Latitude.x + Longitude.x + (1|Block), data=BOGR.cl)
#Anova(BOGR.ht.mod)
#plot(allEffects(BOGR.ht.mod))
#plot(predictorEffects(BOGR.ht.mod))

BOGR.ht.bv.mod <- lmer(Length_cm_20220801 ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
                       + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
                       + scale(bio17) + (1|Block), data=BOGR.cl)
Anova(BOGR.ht.bv.mod)
plot(allEffects(BOGR.ht.bv.mod))
plot(predictorEffects(BOGR.ht.bv.mod))

BOGR.inf.bv.mod <- glmer(NumInf ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
                       + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
                       + scale(bio17) + (1|Block), data=BOGR.cl, family=poisson(link="log"))
Anova(BOGR.inf.bv.mod)
plot(allEffects(BOGR.inf.bv.mod))
plot(predictorEffects(BOGR.inf.bv.mod))

## ** Pick different model form? **
#BOGR.df.bv.mod <- lmer(log(DaysToFlwr) ~ scale(bio1) + scale(bio2) + scale(bio3) + scale(bio5)
#                        + scale(bio8) + scale(bio12) + scale(bio13) + scale(bio14) + scale(bio15)
#                        + scale(bio17) + (1|Block), data=BOGR.cl)
#Anova(BOGR.df.bv.mod)
#plot(allEffects(BOGR.df.bv.mod))
#plot(predictorEffects(BOGR.df.bv.mod))
## --------------------------------------------------------------


## ALTERNATIVELY - Look at PCA of 19 variables to reduce number of predictors -------------------
BOGR.biovar <- left_join(BOGR.biovar, BOGR.SdZn, by="Source")

#BOGR.biovarCen <- scale(BOGR.biovar[,2:20], center=TRUE, scale=TRUE)
#BOGR.biovarComb <- as.data.frame(cbind(BOGR.biovar[,21],BOGR.biovarCen))
BOGR.pcaBiovar <- prcomp(BOGR.biovar[,2:20], scale=TRUE)
par(pty="s")
plot(x=BOGR.pcaBiovar$x[,1], y=BOGR.pcaBiovar$x[,2], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,2], y=BOGR.pcaBiovar$x[,3], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,3], y=BOGR.pcaBiovar$x[,4], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,4], y=BOGR.pcaBiovar$x[,5], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)
plot(x=BOGR.pcaBiovar$x[,5], y=BOGR.pcaBiovar$x[,6], pch=19, cex=1.4, col=BOGR.biovar$SdZnColful)

## Look at scree plot
BOGR.pcaBiovar$sdev[1]**2/sum(BOGR.pcaBiovar$sdev**2)
BOGR.pcaBV.varExpl <- BOGR.pcaBiovar$sdev^2/sum(BOGR.pcaBiovar$sdev^2)
barplot(BOGR.pcaBiovar$sdev[1:19]**2/sum(BOGR.pcaBiovar$sdev**2))
sum(BOGR.pcaBV.varExpl[1:3]) #top 6 PC axes explain over 98% of variation

## Add arrows on PCA plot and look at loadings 
biplot(BOGR.pcaBiovar)
BOGR.pcaBiovar$rotation
#BOGR.pc1 <- BOGR.pcaBiovar$rotation[,1][order(BOGR.pcaBiovar$rotation[,1])]
#BOGR.pc2 <- BOGR.pcaBiovar$rotation[,2][order(BOGR.pcaBiovar$rotation[,2])]

## Use a cumulative loading score to select biovars to include  
BOGR.ldSum <- rowSums(abs(BOGR.pcaBiovar$rotation[,1:3]))
BOGR.ldSum <- BOGR.ldSum[order(BOGR.ldSum)]
## ----------------------------------------


## Models
## Using PC axis values
#BOGR.pcs <-as.data.frame(cbind(BOGR.biovar$Source, BOGR.pcaBiovar$x))
#colnames(BOGR.pcs) <- c("Source", colnames(BOGR.pcs[2:ncol(BOGR.pcs)]))
#BOGR.cl <- left_join(BOGR.cl, BOGR.pcs[1:5], by="Source") #Add PC1-PC4
#BOGR.cl$PC1 <- as.numeric(BOGR.cl$PC1)
#BOGR.cl$PC2 <- as.numeric(BOGR.cl$PC2)
#BOGR.cl$PC3 <- as.numeric(BOGR.cl$PC3)
#BOGR.cl$PC4 <- as.numeric(BOGR.cl$PC4)

#BOGR.ht.pc.mod <- lmer(Length_cm_20220801 ~ PC1 + PC2 + PC3 + PC4 + (1|Block), data=BOGR.cl)
#Anova(BOGR.ht.pc.mod)
#plot(allEffects(BOGR.ht.pc.mod))
#plot(predictorEffects(BOGR.ht.pc.mod))

## Use 'top' biovars as determined by PCA
 
BOGR.inf.bv.mod <- glmer(NumInf ~ scale(bio2) + scale(bio8) + scale(bio14) + scale(bio19)
                         + (1|Block), data=BOGR.cl, family=poisson(link="log"))
Anova(BOGR.inf.bv.mod)
plot(allEffects(BOGR.inf.bv.mod))
plot(predictorEffects(BOGR.inf.bv.mod))

## ** Pick different model form? Negative binomial? **
hist(BOGR.cl$DaysToFlwr)
BOGR.df.bv.mod <- lmer(log(DaysToFlwr) ~ scale(bio2) + scale(bio8) + scale(bio14) + scale(bio19)
                       + (1|Block), data=BOGR.cl)
Anova(BOGR.df.bv.mod)
plot(allEffects(BOGR.df.bv.mod))
plot(predictorEffects(BOGR.df.bv.mod))
## ---------------------------------------


## Visualize raw data relationships with bioclimate vars
plot(BOGR.cl$Length_cm_20220801 ~ BOGR.cl$bio2)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio2)
plot(BOGR.cl$DaysToFlwr ~ BOGR.cl$bio8)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio14)
plot(BOGR.cl$DaysToFlwr ~ BOGR.cl$bio19)
plot(BOGR.cl$NumInf ~ BOGR.cl$bio8)


## Calculate trait means
BOGR.summ <- BOGR.cl %>% group_by(Source) %>% 
             dplyr::summarise(DaysToFlwr_MN=mean(DaysToFlwr, na.rm=TRUE), DaysToFlwr_SE=calcSE(DaysToFlwr),
                   Height_MN=mean(Length_cm_20220801, na.rm=TRUE), Height_SE=calcSE(Length_cm_20220801),
                   Inf_MN=mean(NumInf, na.rm=TRUE), Inf_SE=calcSE(NumInf), 
                   Growth_MN=mean(GrwthRate_Relative, na.rm=TRUE), Growth_SE=calcSE(GrwthRate_Relative), 
                   BIO1=mean(bio1, na.rm=TRUE),BIO2=mean(bio2, na.rm=TRUE),BIO3=mean(bio3, na.rm=TRUE),BIO4=mean(bio4, na.rm=TRUE), 
                   BIO5=mean(bio5, na.rm=TRUE),BIO6=mean(bio6, na.rm=TRUE),BIO7=mean(bio7, na.rm=TRUE),BIO8=mean(bio8, na.rm=TRUE), 
                   BIO9=mean(bio9, na.rm=TRUE),BIO10=mean(bio10, na.rm=TRUE),BIO11=mean(bio11, na.rm=TRUE),BIO12=mean(bio12, na.rm=TRUE),
                   BIO13=mean(bio13, na.rm=TRUE),BIO14=mean(bio14, na.rm=TRUE),BIO15=mean(bio15, na.rm=TRUE),BIO16=mean(bio16, na.rm=TRUE),
                   BIO17=mean(bio17, na.rm=TRUE),BIO18=mean(bio18, na.rm=TRUE),BIO19=mean(bio19, na.rm=TRUE))
BOGR.summ <- left_join(BOGR.summ, BOGR.SdZn, by="Source")

#par(mfrow=c(2,2))
#Loop over columns with biovars
for (bb in 10:28) { 
  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Height_MN, col=BOGR.summ$PopCol, pch=19, cex=1.2, main="B. gracilis", 
     cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Height (cm)", cex.lab=1.5, cex.axis=1.1)
  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Height_MN-BOGR.summ$Height_SE, as.vector(t(BOGR.summ[,bb])), 
         BOGR.summ$Height_MN+BOGR.summ$Height_SE, 
         angle=90, col=BOGR.summ$PopCol, code=3, length=0, lwd=1.6) }

for (bb in 10:28) { 
  plot(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Growth_MN, col=BOGR.summ$PopCol, pch=19, cex=1.2, main="B. gracilis", 
       cex.main=1.5, xlab=colnames(BOGR.summ[bb]), ylab="Relative plant growth", cex.lab=1.5, cex.axis=1.1)
  arrows(as.vector(t(BOGR.summ[,bb])), BOGR.summ$Growth_MN-BOGR.summ$Growth_SE, as.vector(t(BOGR.summ[,bb])), 
         BOGR.summ$Growth_MN+BOGR.summ$Growth_SE, 
         angle=90, col=BOGR.summ$PopCol, code=3, length=0, lwd=1.6) }

for (bb in 10:28) { 
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
#Don't estimate days to mort since plants appear and disappear often within and across years
#Don't use number of days flowering/ length of flowering in 2023 as phenology data is not always consistent
#Just use data to determine number of days until first flower 
#Calculate percent flowering for each pop. I.e. how many indivs of each pop flowered. Were there some pops with no flowering? 
#Could try removing any questionable indivs, eg. any with inconsistent data entry (some highlighted in xls file)...
#If enough indivs remain, repeat analyses and check for consistency 
#Add amt repro into datasheet once fruit counts done in lab
#Add specific leaf data (leaf area and weight) into raw datasheet (csv) or combine in R




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


## PEVI ---------------
boxplot(Length_cm_20220707 ~ Source, data=PEVI, ylim=c(0,10),
        xlab="Population", ylab="Plant height (cm)", cex.lab=1.5,
        names=c("P.CO1", "P.CO2", "P.CO3", "P.CO4", "P.WY1", "P.WY2"),
        main="Penstemon virens", cex.main=1.5, col="thistle", cex.axis=1.25)
## ------------------------------------------------------------------------------------------------











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

#ARFR.traits <- ARFR.cl %>% dplyr::select(c("Length_cm_20220726","Survival_20220922","Phenology_20220922",
#                                           "AGB_MinusBag","InfBM_Wbag"))
#covMat.traits <- cov(ARFR.traits, use="pairwise.complete.obs")
#pca.results <- prcomp(covMat.traits)
#cols <- viridis(5)
#plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, col=cols, cex=1.2)
#legend("topleft", colnames(ARFR.traits), col=cols, cex=0.75, pch=19)
#ARFR.popMns <- ARFR.cl %>% group_by(Source) %>% dplyr::summarise(LEN_EARLY=mean(Length_cm_20220622,na.rm=TRUE), 
#                                                LEN_LATE=mean(Length_cm_20220726,na.rm=TRUE), SURV_LATE=mean(Survival_20220922,na.rm=TRUE),
#                                                AGB=mean(AGB_MinusBag,na.rm=TRUE), REPRO_BM=mean(InfBM_Wbag,na.rm=TRUE))
#ARFR.popMnsT <- ARFR.popMns %>% pivot_longer(cols=-1) %>% pivot_wider(names_from="Source", values_from="value") %>% 
#  rename(Trait=name) #Transpose 
#covMat.pops <- cov(as.data.frame(ARFR.popMnsT[,2:12]), use="pairwise.complete.obs")
#pca.pops <- prcomp(covMat.pops)
#cols <- viridis(11)
#plot(x=pca.pops$x[,1], y=pca.pops$x[,2],pch=19, col=cols, cex=1.2)
#legend("topleft", colnames(ARFR.popMnsT[2:12]), col=cols, cex=0.5, pch=19)

##Look at pops
#ARFR.popsByBios <- t(ARFR.justBios)
#covMat.popsByBios <- cov(ARFR.popsByBios)
#pca.popsByBios <- prcomp(covMat.popsByBios)
#plot(x=pca.popsByBios$x[,1], y=pca.popsByBios$x[,2],pch=19, cex=1.2)
#var.expl <- pca.popsByBios$sdev^2/sum(pca.popsByBios$sdev^2)
#biplot(pca.popsByBios)

# filter ARFR.biovar to keep only biovar columns
#ARFR.justBios <- ARFR.biovar %>% dplyr::select(c(starts_with("bio")))
#covMat <- cov(ARFR.justBios) #Use cov an prcomp
#pca.results <- prcomp(covMat)
#cols <- viridis(19)
#plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, col=cols, cex=1.2)
#legend("topleft", colnames(ARFR.justBios), col=cols, cex=0.75, pch=19)
#plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, col="yellow", cex=1.2)
#text(pca.results$x[,1], pca.results$x[,2], colnames(ARFR.justBios), cex=0.7)
#plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, col="yellow", cex=1.2, xlim=c(-1350,-550), ylim=c(-0,300))
#text(pca.results$x[,1], pca.results$x[,2], colnames(ARFR.justBios), cex=0.7)

#plot(x=pca.results$x[,1], y=pca.results$x[,3],pch=19, col=cols, cex=1.2)
#text(pca.results$x[,1], pca.results$x[,3], colnames(ARFR.justBios), cex=0.7)
#legend("topleft", colnames(ARFR.justBios), col=cols, cex=0.75, pch=19)

#plot(x=pca.results$x[,2], y=pca.results$x[,3],pch=19, col=cols, cex=1.2)
#legend("topleft", colnames(ARFR.justBios), col=cols, cex=0.75, pch=19)
#plot(x=pca.results$x[,2], y=pca.results$x[,3],pch=19, col="yellow", cex=1.2)
#text(pca.results$x[,2], pca.results$x[,3], colnames(ARFR.justBios), cex=0.7)

#biplot(pca.results)
## Try PCA using different (not cov) approach 

#library(ggfortify)
#autoplot(BOGR.pcaBiovar, data=BOGR.biovar, colour=BOGR.biovar$SdZnColful, loadings=TRUE,
#         loadings.label=TRUE, cex=1.2)
#autoplot(BOGR.pcaBiovar, x=2, y=3, colour=BOGR.biovar$PopCol, loadings=TRUE,
#         loadings.label=TRUE)
#autoplot(BOGR.pcaBiovar, x=2, y=3, colour=BOGR.biovar$PopCol, loadings=TRUE,
#         loadings.label=TRUE, frame=TRUE)#, shape=BOGR.biovar$PopCol)
#panel.background=element_rect(fill="white")
#library(factoextra)
#fviz_pca_biplot(BOGR.pcaBiovar, col.ind=BOGR.biovar$SdZnColful, col.var="black", repel=TRUE,
#                geom="point")
#BOGR.cvH$TraitCol <- "blue3"


## *** Need to work on this **** (move unused to ChatfieldCGdataAnalysis.R)
#Or to start, just look at all rows with inconsistent survival data. Then either mark as Remove or correct error if obvious
## CONSOLIDATE SURVIVAL DATA
#ARFR22.surv <- ARFR22.cl %>% dplyr::select(c(starts_with("Survival_")))
#for (rr in 1:nrow(ARFR.Surv)) {
#  for (cc in 1:(ncol(ARFR.Surv)-1)) {
#    if(ARFR.Surv[rr,cc]==0 & ARFR.Surv[rr,cc+1]==1) {
#     ARFR.Surv$Check[rr] <- "Y"
#    }
#  }
#}
#ARFR.Surv %>% filter_all(any_vars(.==0))
#ARFR.cl %>% filter_all(starts_with("Survival_")==0)

#ARFR22.surv <- ARFR22.surv %>% mutate(Check = ifelse(Survival_20220715 - Survival_20220622 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check2 = ifelse(Survival_20220721 - Survival_20220715 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check2=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check3 = ifelse(Survival_20220726 - Survival_20220721 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check3=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check4 = ifelse(Survival_20220804 - Survival_20220726 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check2=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check5 = ifelse(Survival_20220817 - Survival_20220804 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check5=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check6 = ifelse(Survival_20220830 - Survival_20220817 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check6=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check7 = ifelse(Survival_20220909 - Survival_20220830 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check7=="Remove?",]
#ARFR22.surv <- ARFR22.surv %>% mutate(Check8 = ifelse(Survival_20220922 - Survival_20220909 <= 0, "", "Remove?"))
#ARFR22.surv[ARFR22.surv$Check8=="Remove?",]