## April Goebl
## Script modified 2025-11-02 
## BLM Restoration project at Denver Botanic Gardens
## Analyze ARFR data from Chatfield Common Garden  


rm(list=ls())
#dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
#library(Hmisc)
#library(tidyr)
#library(car)
#library(corrplot)
#library(tidyverse)
library(dplyr)
library(stringr)
library(lme4)
library(plotrix)
library(EnvStats)
library(effects)
library(reshape2)
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
## ** From GSD **
## Emergence 
## Logistic regression 
#fitEmrg.mm <- glmer(cbind(emrg, num.planted-emrg) ~ ecotype * habitat + (1|pop) + (1|plot), 
#                    data=datEarly, family=binomial(link="logit"),
#                    glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

#fitEmrg.src <- glm(cbind(emrg, num.planted-emrg) ~ ecotype * habitat, data=datEarly,
#                   family=binomial(link="logit"))



## Seedling-to-adult survival 
## Logistic regression 
#fitSurv.mm <- glmer(cbind(surv, emrg-surv) ~ ecotype * habitat + (1|pop) + (1|plot), 
#                    data=datEarly, family=binomial(link="logit"),
#                    glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000)))

#fitSurv.src <- glm(cbind(surv, emrg-surv) ~ ecotype * habitat, data=datEarly,
#                   family=binomial(link="logit"))



## Number of seeds produced (fecundity)
## Negative binomial regression
#fitSds.mm <- glmer.nb(round(numSdsAdj) ~ ecotype * habitat + (1|pop) + (1|plot), data=datLate)

#fitSds.src <- glm.nb(round(numSdsAdj) ~ ecotype * habitat, data=datLate)  


## Model SLA and plt size with linear models. Also try growth rate? 
## Don't show SLA or growth rate but say that we models them and no sig diffs between sources? 


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
ARFR.traits <- ARFR.cl %>% dplyr::select(c("Surv_2022","Survival","Length_cm_20220726", "Height_20230927", "SLA_mm2permg",
                                           "AGB2022_MinusBag","InfBM2022_2024updated", "DryLeafMass_g", "LeafSurfaceArea_cm2")) 
                                            #Add growth rate(s)?

ARFR.traits <- ARFR.traits[!is.na(ARFR.traits$Length_cm_20220726) & ARFR.traits$Surv_2022==1,] #Remove indivs that died early & have no data
#ARFR.traits <- ARFR.traits[,2:8] #Remove survival
ARFR.traitsT <- t(ARFR.traits)

## Get sample list with pop ID and colors
ARFR.indivPop <- ARFR.cl %>% dplyr::select(c("Source", "ID", "HexCode_Indv"))
ARFR.indivPop$ID <- as.factor(ARFR.indivPop$ID)
indivs.traitPCA <- as.factor(colnames(ARFR.traitsT))
indivs.traitPCA <- as.data.frame(indivs.traitPCA)
colnames(indivs.traitPCA) <- "ID"
indivs.traitPCA <- left_join(indivs.traitPCA, ARFR.indivPop, by="ID")

## If some individuals have NAs in all columns, probably no data for any traits for these samples
#covMat.traitsNoNA <- covMat.traits %>% dplyr::drop_na(covMat.traits)
#filter(covMat.traits, rowSums(is.na(covMat.traits)) != ncol(covMat.traits))
#covMat.traits[complete.cases(as.data.frame(covMat.traits)),] #(na.omit(covMat.traits))

## Make covariance matrix and run pca
covMat.traits <- cov(ARFR.traitsT, use="pairwise.complete.obs")
pca.results <- prcomp(covMat.traits)


par(mfrow=c(1,1))
#cols <- viridis(9)
plot(x=pca.results$x[,1], y=pca.results$x[,2],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv, main="Trait PCA", ylim=c(-15000,0))
plot(x=pca.results$x[,2], y=pca.results$x[,3],pch=19, cex=1.2, col=indivs.traitPCA$HexCode_Indv)
#legend("topleft", colnames(ARFR.traits), col=cols, cex=0.75, pch=19)
## ** Look into why weird lines 
## ** Look into loadings (which traits contribute most to PC1)

## Calculate PC1 mean values for each source population
trait.PCscores <- as.data.frame(cbind(pca.results$x[,1], as.character(indivs.traitPCA$Source)))
colnames(trait.PCscores) <- c("PC1", "Source")
trait.PCscores$PC1 <- as.numeric(trait.PCscores$PC1)
traitPC1.mean <- trait.PCscores %>% group_by(Source) %>% summarise(PC1mean = mean(PC1), n=n())


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


## Plot range of color gradient as a legend
traitPCrange <- seq(from=min(traitPC1.mean$PC1mean), to=max(traitPC1.mean$PC1mean), by=1000)
vals_normRange <- (traitPCrange - min(traitPCrange)) / (max(traitPCrange) - min(traitPCrange))
rgb_matrixRange <- gradient_fn(vals_normRange)
colors.traitPCrange <- rgb(rgb_matrixRange[,1], rgb_matrixRange[,2], rgb_matrixRange[,3], maxColorValue = 255)
plot(traitPCrange, rep(0.25, length(traitPCrange)), col=colors.traitPCrange, pch=15, cex=4)


## ** ** Make different from genomic PC gradient ** 
## ** Try black and white ** 

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

## ** Look into what 'NA' is can clean up ****
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









