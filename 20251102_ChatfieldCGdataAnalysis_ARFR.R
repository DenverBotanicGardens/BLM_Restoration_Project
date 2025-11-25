## April Goebl
## Script modified 2025-11-02 
## BLM Restoration project at Denver Botanic Gardens
## Analyze ARFR data from Chatfield Common Garden  


rm(list=ls())
#dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
#library(Hmisc)
#library(tidyr)
#library(corrplot)
#library(tidyverse)
library(dplyr)
library(stringr)
library(lme4)
library(car)
library(plotrix)
#library(EnvStats)
library(effects)
#library(reshape2)
library(gplots)
library(AICcmodavg)
library(PerformanceAnalytics)
#library(vcfR)
#library(pcadapt)
calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------





## SET WORKING DIRECTORY --------------------------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
ARFR22 <- read.csv(file="Chatfield/2022_data/20251031_ChatfieldDataClean2022_ARFR.csv", sep=",", header=TRUE, dec=".")#, na.strings="")
ARFR23 <- read.csv(file="Chatfield/2023_data/20251031_ChatfieldDataClean2023_ARFR.csv", sep=",", header=TRUE, dec=".")#, na.strings="")
ARFR24 <- read.csv(file="Chatfield/2024_data/20251031_ChatfieldDataClean2024_ARFR.csv", sep=",", header=TRUE, dec=".")#, na.strings="")

## Load PCA values from Alyson's analysis
pca_vals <- read.csv(file="DNA_Seq/AlysonEmery/20251006_pcaTableFromAlyson_ARFR.csv", sep=",", header=TRUE, dec=".")
## ----------------------------------------------------------------------------------------------




## ARFR - DATA RE-FORMAT AS NEEDED --------------------------------------------------------------
str(ARFR22)
str(ARFR23)
str(ARFR24)

ARFR22$Source <- as.factor(ARFR22$Source)
ARFR23$Source <- as.factor(ARFR23$Source)
ARFR24$Source <- as.factor(ARFR24$Source)

ARFR22$Treatment <- as.factor(ARFR22$Treatment)
ARFR23$Treatment <- as.factor(ARFR23$Treatment)
ARFR24$Treatment <- as.factor(ARFR24$Treatment)
## ----------------------------------------------------------------------------------------------




## ----------------------------------------------------------------------------------------------
## ARFR - PREPARE RESPONSE VARIABLES

## CHANGE DATA TO NA BASED ON VARIOUS CONDITIONS (E.G. EXCLUDE SURV DATA IN 2023-24 IF HARVESTED IN 2022)
ARFR22.ex <- ARFR22 %>% mutate(across(c(starts_with("Survival_"),starts_with("Length_")), 
             ~case_when(ExcludeBcNotReplaced=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))
ARFR23.ex <- ARFR23 %>% mutate(across(c(starts_with("Survival_"), "Height_20230927"), 
             ~case_when(ARFR22.ExcludeBcNotReplaced=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))
ARFR23.ex <- ARFR23.ex %>% mutate(across(c(starts_with("Survival_"), "Height_20230927"), 
             ~case_when(ARFR24.ExcludeBcHarvest=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24 %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g",
             "InfBM2022smpls_HEADS_2024weigh","InfBM2022smpls_CHAFF_2024weigh"), 
            ~case_when(ARFR22.ExcludeBcNotReplaced=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24.ex %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g",
             "InfBM2022smpls_HEADS_2024weigh","InfBM2022smpls_CHAFF_2024weigh"), 
            ~case_when(ARFR23.ExcludeSurvDueToInconsistData=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))
ARFR24.ex <- ARFR24.ex %>% mutate(across(c("Survival","LeafSurfaceArea_cm2","SLA_mm2permg","DryLeafMass_g"), 
            ~case_when(ExcludeBcHarvest=="Y" ~ as.integer(NA), TRUE ~ as.numeric(.x))))





## GROWTH RATE VARIABLES -----------------------------------------------------------------------
## 'Late' growth (June to July)
ARFR22.ex$GrwthRate_Specific <- log(ARFR22.ex$Length_cm_20220726/ARFR22.ex$Length_cm_20220622)
ARFR22.ex$GrwthRate_Absolute <- ARFR22.ex$Length_cm_20220726-ARFR22.ex$Length_cm_20220622
ARFR22.ex$GrwthRate_Relative <- (ARFR22.ex$Length_cm_20220726-ARFR22.ex$Length_cm_20220622)/ARFR22.ex$Length_cm_20220622

## 'Early' growth (May to June)
## **Need to handle pre-replacement early height measurement 20220527 differently.. i.e. exclude replaced plts
ARFR22.ex$GrwthRateErly_Specific <- log(ARFR22.ex$Length_cm_20220622/ARFR22.ex$Length_cm_20220527)
ARFR22.ex$GrwthRateErly_Absolute <- ARFR22.ex$Length_cm_20220622-ARFR22.ex$Length_cm_20220527
ARFR22.ex$GrwthRateErly_Relative <- (ARFR22.ex$Length_cm_20220622-ARFR22.ex$Length_cm_20220527)/ARFR22.ex$Length_cm_20220527




## COMBINE AND ADD REPRO BM DATA --------------------------------------------------------------
## Combine flwr head and chaff/seed weights + any missed smpls from initial 2022 weights
## Change Chaff entries to zero (from NA) if no chaff weight, but heads were weighed
ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24.ex$InfBM2022smpls_HEADS_2024weigh) & is.na(ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh)] <- 0
#length(ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh[!is.na(ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh)]) 
## Add chaff and head weights together
ARFR24.ex$InfBM2022_2024updated <- ARFR24.ex$InfBM2022smpls_HEADS_2024weigh + ARFR24.ex$InfBM2022smpls_CHAFF_2024weigh
length(ARFR24.ex$InfBM2022_2024updated[!is.na(ARFR24.ex$InfBM2022_2024updated)])
## Add several indivs (524, 885, and 908) from 2022 weights that weren't available for 2024 re-weigh
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==524] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==524]
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==885] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==885]
ARFR24.ex$InfBM2022_2024updated[ARFR24.ex$ID==908] <- ARFR24.ex$InfBM2022_Wobag_g[ARFR24.ex$ID==908]



## CLEAN 2023 PLT SZ FIELD MEASUREMENTS -------------------------------------------------------
identical(ARFR23.ex$ExcludeSurvDueToInconsistData, ARFR23.ex$ExcludeSzDueToUncertainty)
nrow(ARFR23.ex[!is.na(ARFR23.ex$ExcludeSurvDueToInconsistData),])
nrow(ARFR23.ex[!is.na(ARFR23.ex$ExcludeSzDueToUncertainty),])

nrow(ARFR23.ex[!is.na(ARFR23.ex$Height_20230927>0),])
ARFR23.ex$Height_20230927[ARFR23.ex$ExcludeSzDueToUncertainty=="Y" & !is.na(ARFR23.ex$ExcludeSzDueToUncertainty)] <- NA
ARFR23.ex$Height_20230927[ARFR23.ex$ExcludeSurvDueToInconsistData=="Y" & !is.na(ARFR23.ex$ExcludeSurvDueToInconsistData)] <- NA
nrow(ARFR23[!is.na(ARFR23$Height_20230927>0),])


## CALC SURVIVAL FOR EACH YEAR HERE? **



## COMBINE RELEVANT 2022, 2023, 2024 DATA
# consider adding other lengths and early growth from 2022
ARFR22.sel <- ARFR22.ex %>% dplyr::select(c("ID", "Length_cm_20220726","GrwthRate_Relative","GrwthRateErly_Specific", 
                                            "GrwthRateErly_Absolute","GrwthRateErly_Relative","Survival_20220922"))               
ARFR23.sel <- ARFR23.ex %>% dplyr::select(c("ID","Height_20230927","Survival_20230801")) #Don't use surv 9/27 since not all blocks surveyed
ARFR.sel <- left_join(ARFR24.ex, ARFR23.sel, by="ID") 
ARFR.sel <- left_join(ARFR.sel, ARFR22.sel, by="ID") 
ARFR.sel$Source <- as.factor(ARFR.sel$Source)

## Surv checks on combined data? 
## Save?
#write.csv(ARFR.sel, "Chatfield/20251031_ChatfieldDataTraits_ARFR.csv", row.names=FALSE)
## --------------------------------------------------------------------------------------------------





## ARFR - VISUALIZE RAW DATA ---------------------------------------------------------------

## Order populations for plotting 
## Order by average size or lat or other trait(s)
#ARFR.htByMed <- with(ARFR.cl, reorder(Source, Length_cm_20220726, median, na.rm=TRUE))
#ARFR.infByMed <- with(ARFR.cl, reorder(Source, InfBM2022_2024updated, median, na.rm=TRUE))
ARFR.latByMed <- with(ARFR.sel, reorder(Source, Lat, median, na.rm=TRUE))

ARFR.meds <- ARFR.sel %>% group_by(Source) %>% 
             dplyr::summarise(Height22_MD=median(Length_cm_20220726,na.rm=TRUE), AGB22_MD=median(AGB2022_MinusBag,na.rm=TRUE),
             ReproBMrw_MD=median(InfBM2022_2024updated,na.rm=TRUE), Height23_MD=median(Height_20230927,na.rm=TRUE),
             Surv24_Count=length(na.omit(Survival)), Surv24_Sum=sum(Survival, na.rm=TRUE),
             Surv23_Count=length(na.omit(Survival_20230801)), Surv23_Sum=sum(Survival_20230801, na.rm=TRUE),
             Surv22_Count=length(na.omit(Survival_20220922)), Surv22_Sum=sum(Survival_20220922, na.rm=TRUE))
             #LeafArea_MD=median(LeafSurfaceArea_cm2, na.rm=TRUE), LeafMass_MD=median(DryLeafMass_g, na.rm=TRUE),
             #GrowthReE_MD=median(GrwthRateErly_Relative,na.rm=TRUE), SLA_MD=median(SLA_mm2permg,na.rm=TRUE),
             #GrowthSpE_MD=median(GrwthRateErly_Specific,na.rm=TRUE), GrowthAbE_MD=median(GrwthRateErly_Absolute,na.rm=TRUE))#Surv24_MN=mean(Survival, na.rm=TRUE)

#ARFR.meds <- left_join(ARFR.meds, ARFR.SdZn, by="Source") 
AddnCols <- as.data.frame(cbind(ARFR.sel$PopAbbrev,ARFR.sel$PopCol,as.numeric(ARFR.sel$Lat),as.character(ARFR.sel$Source)))
colnames(AddnCols) <- c("PopAbbrev","PopCol","Lat","Source")
ARFR.meds <- left_join(ARFR.meds, AddnCols, by="Source")
ARFR.meds <- unique(ARFR.meds)

## Estimate survival each year
surv24.pop <- ARFR.meds$Surv24_Sum/ARFR.meds$Surv24_Count
surv23.pop <- ARFR.meds$Surv23_Sum/ARFR.meds$Surv23_Count
surv22.pop <- ARFR.meds$Surv22_Sum/ARFR.meds$Surv22_Count



## Boxplots of raw data 
ARFR.meds <- ARFR.meds[order(ARFR.meds$Lat),] #Order by lat

par(mfrow=c(2,3))

## Size 2022
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Height22_MD),] #Order by median sz
boxplot(Length_cm_20220726 ~ ARFR.latByMed, data=ARFR.sel,
        ylab="Height (cm)", xlab=NA, cex.lab=1.25, horizontal=FALSE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2,
        main="FINAL SIZE 2022", cex.main=1.5, col=ARFR.meds$PopCol)


## Growth rate(s) 
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Growth_MD),]
boxplot(GrwthRateErly_Relative ~ ARFR.latByMed, data=ARFR.sel, las=2,
        xlab="Plant early relative growth", ylab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, ylim=c(-0.5,5.5),
        cex.main=1.5, col=ARFR.meds$PopCol, ylim=c(-0.4,7), main="GROWTH RATE 2022")

boxplot(GrwthRateErly_Absolute ~ ARFR.latByMed, data=ARFR.sel, las=2,
        xlab=NA, ylab="Plant early absolute growth", cex.lab=1, cex.axis=0.9, names=ARFR.meds$PopAbbrev,
        cex.main=1.5, col=ARFR.meds$PopCol, ylim=c(-0.9,20))

boxplot(GrwthRateErly_Specific ~ ARFR.latByMed, data=ARFR.sel,
        xlab="Population", ylab="Plant early specific growth", cex.lab=1.5, names=ARFR.meds$PopAbbrev,
        cex.axis=0.7, cex.main=1.5, col=ARFR.meds$PopCol)


## Repro
boxplot(InfBM2022_2024updated ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="Reproductive biomass (g)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,80),
        cex.main=1.5, col=ARFR.meds$PopCol, main="REPRODUCTIVE OUTPUT 2022")


## Size 2023
#ARFR.meds <- ARFR.meds[order(ARFR.meds$Height23_MD),] #Order by median 2023 sz
boxplot(Height_20230927 ~ ARFR.latByMed, data=ARFR.sel,
        ylab="Height (cm)", xlab=NA, cex.lab=1.25, horizontal=FALSE,
        cex.axis=0.99, names=ARFR.meds$PopAbbrev, las=2, ylim=c(15,90),
        main="FINAL SIZE 2023", cex.main=1.5, col=ARFR.meds$PopCol)


## SLA 2024
boxplot(SLA_mm2permg ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="Specific leaf area (mm2/mg)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,40),
        cex.main=1.5, col=ARFR.meds$PopCol, main="SPECIFIC LEAF AREA 2024")
boxplot(LeafSurfaceArea_cm2 ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="leaf area (cm2)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,2.3),
        cex.main=1.5, col=ARFR.meds$PopCol, main="LEAF AREA 2024")
boxplot(DryLeafMass_g ~ ARFR.latByMed, data=ARFR.sel, las=2,
        ylab="leaf mass (g)", xlab=NA, cex.lab=1.25, cex.axis=0.99, 
        names=ARFR.meds$PopAbbrev, horizontal=FALSE, ylim=c(0,0.02),
        cex.main=1.5, col=ARFR.meds$PopCol, main="LEAF MASS 2024")


## Survival 2024
barplot(surv.pop, col=ARFR.meds$PopCol, ylim=c(0,1), cex.axis=0.99, names.arg=ARFR.meds$PopAbbrev,
        las=2, ylab="Survival rate", main="SURVIVAL 2022-2024", cex.main=1.5)

## Try stacked bar plot with survival rate by year ** 
## Make stacked barplot to visualize cv for each trait and population
#ERNA.cvT <- as.matrix(rbind(as.vector(ERNA.cv$Height_CV), as.vector(ERNA.cv$Growth_CV),
#                            as.vector(ERNA.cv$DaysToFlwr_CV)))#, as.vector(ERNA.cv$GrowthE_CV)
#colnames(ERNA.cvT) <- ERNA.cv$Source

#barplot(ERNA.cvT, names=ERNA.cv$PopAbbrev, las=2, col=c("black","dodgerblue","bisque3"),
#        ylab="Coefficient of variation", cex.lab=1.3)




## Blank plot
#plot.new()
#legend("center", unique(ARFR.meds$Source[order(ARFR.meds$PopOrder, decreasing=TRUE)]), 
#       col=unique(ARFR.meds$PopCol[order(ARFR.meds$PopOrder, decreasing=TRUE)]), cex=1.1, pch=19)

## Generate legend figure with seed zone colors used in map and seed zone abbrevs 
#legend("center", unique(ARFR.meds$SdZnAbbrev[order(ARFR.meds$PopOrder, decreasing=TRUE)]), 
#       col=unique(ARFR.meds$SdZnCol[order(ARFR.meds$PopOrder, decreasing=TRUE)]), cex=2, pch=19)
## ---------------------------------------------------





## MODEL TRAIT DATA ----------------------------------

## Re-order Source as factor before running models
AddnCols.unq <- unique(AddnCols)
AddnCols.unq <- AddnCols.unq[order(AddnCols.unq$Lat),] #Order by lat

ARFR.sel$Source <- factor(ARFR.sel$Source, levels=AddnCols.unq$Source)



## Plant size

## 2022
ARFR.sel$Block <- as.factor(ARFR.sel$Block)
sz22.mod <- lmer(Length_cm_20220726 ~ Source + (1|Block), data=ARFR.sel)
summary(sz22.mod)
Anova(sz22.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sz22.mod, type="pearson")
dResid <- residuals(sz22.mod, type="deviance")
hist(pResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
hist(dResid)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(dResid)
plot(fitted(sz22.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line

## Obtain model predicted values for response variables
predForSource <- as.data.frame(AddnCols.unq$Source) #(unique(ARFR.sel$Source))
colnames(predForSource) <- "Source"
sz22.pred <- predict(sz22.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)


## 2023
sz23.mod <- lmer(Height_20230927 ~ Source + (1|Block), data=ARFR.sel)
summary(sz23.mod)
Anova(sz23.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sz23.mod, type="pearson")
hist(pResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(pResid)
plot(fitted(sz23.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line

## Obtain model predicted values for response variables
sz23.pred <- predict(sz23.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)



## SLA
hist(ARFR.sel$SLA_mm2permg)
hist(log(ARFR.sel$SLA_mm2permg))
sla.mod <- lmer(log(SLA_mm2permg) ~ Source + (1|Block), data=ARFR.sel)
summary(sla.mod)
Anova(sla.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(sla.mod, type="pearson")
dResid <- residuals(sla.mod, type="deviance")
hist(dResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(pResid)
plot(fitted(sla.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line

## Obtain model predicted values for response variables
sla.predLog <- predict(sla.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
sla.predOrigFit <- exp(sla.predLog$fit)
sla.predOrigSE <- exp(sla.predLog$se.fit)



## Reproductive biomass
hist(ARFR.sel$InfBM2022_2024updated)
hist(log(ARFR.sel$InfBM2022_2024updated))

rbm.mod <- lmer(log(InfBM2022_2024updated) ~ Source + (1|Block), data=ARFR.sel)
summary(rbm.mod)
Anova(rbm.mod)

## Check distribution of residuals to assess if model form/ family is appropriate
pResid <- residuals(rbm.mod, type="pearson")
hist(pResid)                                          #Shape should be consistent with assumed error distribution (e.g. normal)
qqnorm(pResid)                                        #Points should roughly follow the diagonal line, even at tails
qqline(pResid)
plot(fitted(rbm.mod), pResid, abline(h=0,col="red")) #Residuals should be randomly scattered around 0 line
## Try gamma with logged data? Or something for positive values with left skew **

## Obtain model predicted values for response variables
rbm.predLog <- predict(rbm.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)
rbm.predOrigFit <- exp(rbm.predLog$fit)
rbm.predOrigSE <- exp(rbm.predLog$se.fit)



## Survival
surv24.mod <- glmer(Survival ~ Source + (1|Block), data=ARFR.sel, family=binomial(link="logit"))
summary(surv24.mod)
Anova(surv24.mod)

## Obtain model predicted values for response variables
surv24.pred <- predict(surv24.mod, newdata=predForSource, type="response", re.form=~0, se.fit=TRUE)






## PLOT MODEL ESTIMATES AND SE
predForSource <- dplyr::left_join(predForSource, AddnCols, by="Source")
predForSource <- unique(predForSource)
preds <- cbind(predForSource, sz22.pred$fit, sz22.pred$se.fit, sz23.pred$fit, sz23.pred$se.fit,
               rbm.predOrigFit, rbm.predOrigSE, sla.predOrigFit, sla.predOrigSE, surv24.pred$fit, surv24.pred$se.fit)
#preds <- preds[order(preds$Lat),] #Order by lat

par(mfrow=c(1,1))
plot(NA, NA, xlab="Seed source", ylab="Height (cm)",
     main="FINAL SIZE 2022", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(13.5,48.5))
arrows(1:11, preds$`sz22.pred$fit`+preds$`sz22.pred$se.fit`, 1:11, preds$`sz22.pred$fit`-preds$`sz22.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`sz22.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.45)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)


par(mfrow=c(2,2))
plot(NA, NA, xlab="Seed source", ylab="Reproductive  biomass (g)",
     main="REPRODUCTION 2022", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(0,40))
arrows(1:11, preds$rbm.predOrigFit+preds$rbm.predOrigSE, 1:11, preds$rbm.predOrigFit-preds$rbm.predOrigSE,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$rbm.predOrigFit, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Height (cm)",
     main="FINAL SIZE 2023", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(40,68))
arrows(1:11, preds$`sz23.pred$fit`+preds$`sz23.pred$se.fit`, 1:11, preds$`sz23.pred$fit`-preds$`sz23.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`sz23.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Specific leaf area (mm2/mg)",
     main="SPECIFIC LEAF AREA 2024", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(10,16))
arrows(1:11, preds$sla.predOrigFit+preds$sla.predOrigSE, 1:11, preds$sla.predOrigFit-preds$sla.predOrigSE,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$sla.predOrigFit, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)

plot(NA, NA, xlab="Seed source", ylab="Survival rate",
     main="SURVIVAL 2022-2024", cex.lab=1.25, xaxt='n', xlim=c(1,11), ylim=c(0,1.9))
arrows(1:11, preds$`surv24.pred$fit`+preds$`surv24.pred$se.fit`, 1:11, preds$`surv24.pred$fit`-preds$`surv24.pred$se.fit`,
       angle=90, col="black", code=3, length=0, lwd=2)
points(1:11, preds$`surv24.pred$fit`, col="black", bg=preds$PopCol, pch=21, cex=1.5)
axis(side=1, at=1:11,preds$PopAbbrev, las=2, cex.axis=0.9)



## Look at correlation in trait estimates
plot(preds$`sz22.pred$fit`, preds$`sz23.pred$fit`)
plot(preds$`sz22.pred$fit`, preds$rbm.predOrigFit)
plot(preds$`sz23.pred$fit`, preds$rbm.predOrigFit)



## Add Compact Letter Display to plots to show significant comparisons

library(emmeans)
#sz23.mns <- emmeans(sz23.mod, specs="Source", type="response") #Use emmeans to get model means
library(multcomp)
sz23.glht <- glht(sz23.mod, linfct=mcp(Source="Tukey"), alternative="two.sided")
sz23.cld <- cld(sz23.glht, level=0.05) #adjust="Tukey", Letters=LETTERS,

rbm.glht <- glht(rbm.mod, linfct=mcp(Source="Tukey"), alternative="two.sided")
rbm.cld <- cld(rbm.glht, level=0.05) 

sla.glht <- glht(sla.mod, linfct=mcp(Source="Tukey"))
sla.cld <- cld(sla.glht, level=0.05) 

surv.glht <- glht(surv24.mod, linfct=mcp(Source="Tukey"))
surv.cld <- cld(surv.glht, level=0.05) 

## Try instead just looking at emmeans output and generating cld by hand or putting results in tbl**
sz23.pw <- emmeans(sz23.mod, specs = pairwise ~ Source, type="response")
rbm.pw <- emmeans(rbm.mod, specs = pairwise ~ Source, type="response")
sla.pw <- emmeans(sla.mod, specs = pairwise ~ Source, type="response")
surv.pw <- emmeans(surv24.mod, specs = pairwise ~ Source, type="response")



## Also look at coefficient of variation or sd vs mean? **





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







## Trait PCA ----------------------------------------
ARFR.traits <- ARFR.sel %>% dplyr::select(c("Survival","Length_cm_20220726","Height_20230927", "SLA_mm2permg",
                                           "InfBM2022_2024updated"))#, "DryLeafMass_g", "LeafSurfaceArea_cm2")) 
                                            #Add growth rate(s)? #,"Surv_2022","AGB2022_MinusBag",

ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Length_cm_20220726),] #Remove indivs that died early & have no data
#ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Survival),] #Remove rows without surv data
#block 1, block 2 row 6, block 8

## Look at trait correlations
ARFR.traitsCor <- cor(ARFR.traits, use="pairwise.complete.obs")
corrplot(ARFR.traitsCor)

ARFR.traitsT <- t(ARFR.traits)

## Make covariance matrix and run pca
covMat.traits <- cov(ARFR.traitsT, use="pairwise.complete.obs")
# Remove rows with all missing values
covMat.traits <- covMat.traits[rowSums(is.na(covMat.traits)) != ncol(covMat.traits), ]
test <- covMat.traits[rowSums(is.na(covMat.traits)) < ncol(covMat.traits), ]
#covMat.traits <- na.omit(covMat.traits)
#covMat.traits <- cov(ARFR.traits, use="pairwise.complete.obs")
pca.results <- prcomp(covMat.traits, center=TRUE)#, scale.=TRUE)


## Get sample list with pop ID and colors
ARFR.indivPop <- ARFR.sel %>% dplyr::select(c("Source", "ID", "HexCode_Indv"))
ARFR.indivPop$ID <- as.factor(ARFR.indivPop$ID)
indivs.traitPCA <- as.factor(colnames(ARFR.traitsT))
indivs.traitPCA <- as.data.frame(indivs.traitPCA)
colnames(indivs.traitPCA) <- "ID"
indivs.traitPCA <- left_join(indivs.traitPCA, ARFR.indivPop, by="ID")

## If some individuals have NAs in all columns, probably no data for any traits for these samples
#filter(covMat.traits, rowSums(is.na(covMat.traits)) != ncol(covMat.traits))
#covMat.traits <- covMat.traits[complete.cases(as.data.frame(covMat.traits)),] 



par(mfrow=c(1,1))
#cols <- viridis(9)
plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv, main="Trait PCA")#, ylim=c(-15000,0))
plot(x=pca.results$x[,2], y=pca.results$x[,3],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv)
#legend("topleft", colnames(ARFR.traits), col=cols, cex=0.75, pch=19)
## ** Look into why weird lines 

## Look into loadings (which traits contribute most to PC1)?
#ld <- pca.results$rotation


## Calculate PC1 mean values for each source population
trait.PCscores <- as.data.frame(cbind(pca.results$x[,1], pca.results$x[,2], as.character(indivs.traitPCA$Source)))
colnames(trait.PCscores) <- c("PC1", "PC2", "Source")
trait.PCscores$PC1 <- as.numeric(trait.PCscores$PC1)
traitPC1.mean <- trait.PCscores %>% group_by(Source) %>% summarise(PC1mean = mean(PC1), n=n())
## * Save trait PC values (trait.PCscores) as R object/rds **

## Create color gradient and assign colors based on numeric continuous PC1 mean values
# From ChatGPT
# Define a color gradient (e.g., from blue to red)
gradient_fn <- colorRamp(c("greenyellow",   "deeppink"))

# Normalize your values to [0,1] scale
vals_norm <- (traitPC1.mean$PC1mean - min(traitPC1.mean$PC1mean)) / (max(traitPC1.mean$PC1mean) - min(traitPC1.mean$PC1mean))

# Get RGB colors (as integers 0–255)
rgb_matrix <- gradient_fn(vals_norm)

# Convert to hex color strings
colors.traitPC <- rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3], maxColorValue = 255)

# Plot using colors
plot(traitPC1.mean$PC1mean, rep(1, length(traitPC1.mean$PC1mean)), col=colors.traitPC, pch=16, cex=2)

traitPC1.mean$color <- colors.traitPC
#** save mean color plot to show in supp mat as it shows separation into three groups


## Plot range of color gradient as a legend
traitPCrange <- seq(from=min(traitPC1.mean$PC1mean), to=max(traitPC1.mean$PC1mean), by=1000)
vals_normRange <- (traitPCrange - min(traitPCrange)) / (max(traitPCrange) - min(traitPCrange))
rgb_matrixRange <- gradient_fn(vals_normRange)
colors.traitPCrange <- rgb(rgb_matrixRange[,1], rgb_matrixRange[,2], rgb_matrixRange[,3], maxColorValue = 255)
plot(traitPCrange, rep(0.25, length(traitPCrange)), col=colors.traitPCrange, pch=15, cex=4)


## Make different from genomic PC gradient? 
## Try black and white?

## ---------------------------------------------------







### VCF table and PCA  --------------------------------------------------------------------------
## Get list of sample names from filtered vcf table 
genotype_mxFilt <- vcfR::extract.gt(vcf_filt, as.numeric=TRUE)
indvNames <- as.data.frame(colnames(genotype_mxFilt))
colnames(indvNames) <- "Sample"

#make column with ID using string replace 
indvNames$Temp <- str_replace(indvNames$Sample, "ARFR_", "")
indvNames$ID <- as.integer(str_replace(indvNames$Temp, "_sorted", ""))
#join by ID to get source (pop ID)
indvNames <- left_join(indvNames, ARFR24, by="ID")




## PCADAPT
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

### ** Look into what 'NA' is can clean up ***
PC1.mean <- PC1.mean[1:11,]




## Create color gradient and assign colors based on numeric continuous PC1 mean values
# From ChatGPT
# Define a color gradient (e.g., from blue to red)
#sdZnColsOld <- c("#EEB422", "#EEB422", "#B4EEB4", "#8FBC8F", "#B4EEB4", "#8FBC8F", "#B4EEB4", "#8FBC8F", 
#             "#C1FFC1", "#FFEC8B", "#ffff8b")

#gradient_fn <- colorRamp(c("#EEB422",   "#C1FFC1", "#ffff8b"))
#gradient_fn <- colorRamp(c("#EEB422",   "#C1FFC1"))
gradient_fn <- colorRamp(c("greenyellow",   "deeppink"))

# Normalize your values to [0,1] scale
vals_norm <- (PC1.mean$PC1mean - min(PC1.mean$PC1mean)) / (max(PC1.mean$PC1mean) - min(PC1.mean$PC1mean))

# Get RGB colors (as integers 0–255)
rgb_matrix <- gradient_fn(vals_norm)

# Convert to hex color strings
colors <- rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3], maxColorValue = 255)

# Plot using colors
plot(PC1.mean$PC1mean, rep(1, length(PC1.mean$PC1mean)), col=colors, pch=16, cex=2)

PC1.mean$color <- colors

## Plot range of color gradient as a legend
genPCrange <- seq(from=min(PC1.mean$PC1mean), to=max(PC1.mean$PC1mean), by=0.01)
vals_normGenRange <- (genPCrange - min(genPCrange)) / (max(genPCrange) - min(genPCrange))
rgb_matrixGenRange <- gradient_fn(vals_normGenRange)
colors.genPCrange <- rgb(rgb_matrixGenRange[,1], rgb_matrixGenRange[,2], rgb_matrixGenRange[,3], maxColorValue = 255)
plot(genPCrange, rep(0.5, length(genPCrange)), col=colors.genPCrange, pch=15, cex=4)


## ** make legend of color gradient **

## Not sure if this full works as intended... 
#colorAccording2(
#  x,
#  gradTy = "logGray",
#  nStartOmit = NULL,
#  nEndOmit = "sep",
#  revCol = FALSE,
#  alpha = 1,
#  silent = FALSE,
#  debug = FALSE,
#  callFrom = NULL
#)
#plot(1:11,PC1.mean$PC1mean,pch=16,cex=2,col=colorAccording2(PC1.mean$PC1mean))
#plot(PC1.mean$PC1mean, col=colorAccording2(PC1.mean$PC1mean), pch=19, cex=1.5)
#PC1.mean$HexCode <- colorAccording2(PC1.mean$PC1mean)
#col <- c("#00FF2EFF", "#00FFB9FF", "#FF008BFF", "#5D00FFFF", "#002EFFFF", "#FF0000FF", "#00B9FFFF")
# Create the color ramp function (from AI)
#color_function <- colorRampPalette(c("blue", "red"))
# Generate colors for the data
#PC1.mean$color <- color_function(PC1.mean$PC1mean)
## -------------------------------------------------------------------------------------------









