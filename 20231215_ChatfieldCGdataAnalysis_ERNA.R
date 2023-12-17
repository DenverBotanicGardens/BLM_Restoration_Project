## April Goebl
## Script started 2023-12-15 (modified from 20231122_ChatfieldCGdataAnalysis_BOGR)
## BLM Restoration project at Denver Botanic Gardens
## Analyze ERNA data from Chatfield Common Garden  


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
ERNA <- read.csv(file="Chatfield/20230302_ChatfieldData_ERNA.csv", sep=",", header=TRUE, dec=".")
ERNA.SdZn <- read.csv(file="AGoebl/Seeds/20231215_ERNA_LatLongSdZn_hexcodes.csv", sep=",", header=TRUE, dec=".")
ERNA.biovar <- readRDS("AGoebl/Seeds/20230825_ERNA_BiovarsAvg1980_2021")
## ----------------------------------------------------------------------------------------------




## ERNA - DATA CLEAN UP ---------------------------------------------
str(ERNA)
ERNA$Source <- as.factor(ERNA$Source)

##If OrigPltSurv_20220518 = 0 & plt not replaced (N), ignore data for this plt, i.e. future surv should be NA, not 0
ERNA.cl <- ERNA[ERNA$OrigPltSurvival_20220518==1 | (ERNA$OrigPltSurvival_20220518==0 & ERNA$Replaced_YorN=="Y"),]


## *** CHECK AND EDIT THIS FOR FUTURE ANALYSES !! ** ###
## ** FOR NOW (20231217) DON'T USE EARLY SURV OR LEN IN ANALYSES ***
#Note: Don't use OrigPltSurv_20220518 data in days to mort or other field analyses,
#this surv may not correspond to plt names in Source/Pop col (may correspond to orig planted or assigned)
#What about 6/8 surv?

#Plts that were dead and missing on 5/18 were replaced. Surv & sz was taken on 6/8. 
#If any that were dead on 5/18 were still dead 6/13, they were replaced as well.
#Ignore surv=0 on 6/8 if replaced=Y? 
#ERNA.cl[ERNA.cl$OrigPltSurvival_20220518==0 & (ERNA.cl$Survival_20220608==0 & ERNA.cl$Replaced_YorN=="Y"),]

#Length_0608 should be good to use foir early growth rate? 

#If not replaced, but died before planting, don't use subsequent surv data?
#ERNA.cl[!is.na(ERNA.cl$DateMortalityObservedPreTransplant),] #All plts that died before planting were replaced
## ****************************************************** 


## Checks 
## ** Add other checks listed in BOGR ** 
## Check that surv is only 1, 0 and maybe NA
ERNA.cl[ERNA.cl$Survival_20220608 < 0 | ERNA.cl$Survival_20220608 > 1,]
ERNA.cl[ERNA.cl$Survival_20221108 < 0 | ERNA.cl$Survival_20221108 > 1,]

#Check that surv is only integers
ERNA.cl[ERNA.cl$Survival_20220608 - floor(ERNA.cl$Survival_20220608) != 0,]
ERNA.cl[ERNA.cl$Survival_20221108 - floor(ERNA.cl$Survival_20221108) != 0,]

# ** Check that length is only numeric
# ** Check that if surv=0 for a given date, there are no phenology or height values for that date

## *** Need to work on this ****
#Check that once zero in surv on X/X or later, stays zero (if becomes 1 later, could be data entry error)
## CONSOLIDATE SURVIVAL DATA
#ERNA.Surv <- ERNA.cl %>% dplyr::select(c(starts_with("Survival_")))
#for (rr in 1:nrow(ERNA.Surv)) {
#  for (cc in 1:(ncol(ERNA.Surv)-1)) {
#    if(ERNA.Surv[rr,cc]==0 & ERNA.Surv[rr,cc+1]==1) {
#     ERNA.Surv$Check[rr] <- "Y"
#    }
#  }
#}
#ERNA.Surv %>% filter_all(any_vars(.==0))
#ERNA.cl %>% filter_all(starts_with("Survival_")==0)


#ERNA.test <- ERNA.cl[ERNA.cl$Replaced_YorN !="Y",]

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
## ----------------------------------------------------------------------------------------------




## ARFR - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ERNA.SdZn$Source <- str_replace(ERNA.SdZn$Code, "10-SOS", "")
ERNA.biovar$Source <- str_replace(ERNA.biovar$Pop, "10-SOS", "")


## Edit column names for biovariables
colnames(ERNA.biovar) <- c("Pop","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                           "bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                           "bio17","bio18","bio19","Source")




## Subset of ERNA seed zone colors (all listed in ERNA.SdZn file)
#ERNA.SdZn$HexCode[grepl("0 - 5 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "#F0FFF0"     #semi-humid, v.cold #honeydew #
#ERNA.SdZn$HexCode[grepl("20 - 25 Deg. F. / 12 - 30", ERNA.SdZn$seed_zone)] = "#CD6090" #arid, v.warm #hotpink3 
#ERNA.SdZn$HexCode[grepl("5 - 10 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "#C1FFC1"    #semi-humid, cold #darkseagreen1 
#ERNA.SdZn$HexCode[grepl("15 - 20 Deg. F. / 2 - 3", ERNA.SdZn$seed_zone)] = "#87CEFF"   #humid, warm #skyblue1 
## ---------------------------


## Add seed zone name abbreviation column (look at BOGR script).. (for legend)?

## Add column with seed zone or pop 'order' (look at BOGR script)..?
ERNA.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 2 - 3", ERNA.SdZn$seed_zone)] = "A"    #humid, warm  
ERNA.SdZn$SdZnOrder[grepl("0 - 5 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "B"      #semi-humid, v.cold 
ERNA.SdZn$SdZnOrder[grepl("5 - 10 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "C"     #semi-humid, cold
ERNA.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "D"    #semi-humid, cool
ERNA.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "E"    #semi-humid, warm
ERNA.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "F"    #semi-humid, v.warm
ERNA.SdZn$SdZnOrder[grepl("10 - 15 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "G"   #semi-arid, cool
ERNA.SdZn$SdZnOrder[grepl("15 - 20 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "H"   #semi-arid, warm
ERNA.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "I"   #semi-arid, v.warm
ERNA.SdZn$SdZnOrder[grepl("20 - 25 Deg. F. / 12 - 30", ERNA.SdZn$seed_zone)] = "J"  #arid, v.warm
## ----------------------------------------------------------------------------------------------




## ARFR - COMBINE DATA TYPES --------------------------------------------
ERNA.cl <- left_join(ERNA.cl, ERNA.SdZn, by="Source")
ERNA.cl <- left_join(ERNA.cl, ERNA.biovar, by="Source")
## ----------------------------------------------------------------------



## ERNA - CONSOLIDATE PHENOLOGY DATA AND DO CHECKS ---------------------------------------------
## Add from 20230129_ChatfieldCGdataAnalysis and BOGR script if want to include.
## Exclude for now due to lack of data in 2022 
## ---------------------------------------------------------------



## ERNA - ADD GROWTH RATE VARIABLES ------------------------------
ERNA.cl$GrwthRate_Specific <- log(ERNA.cl$Length_cm_20220915/ERNA.cl$Length_cm_20220719)
ERNA.cl$GrwthRate_Absolute <- ERNA.cl$Length_cm_20220915-ERNA.cl$Length_cm_20220719
ERNA.cl$GrwthRate_Relative <- (ERNA.cl$Length_cm_20220915-ERNA.cl$Length_cm_20220719)/ERNA.cl$Length_cm_20220719

## ** Look at early vs late growth if at least 3 height measurements are usable ** 
## ---------------------------------------------------------------



## LOOK AT REPRO BM DATA FOR 2023
## ---------------------------------------------------------------



## ADD AGB FOR 2023
## ---------------------------------------------------------------




## ERNA - TEST FOR TREATMENT EFFECT -------------------------------------------------------
## Review and update models as needed **
## Copy code from BOGR and ARFR scripts 
## -----------------------------------------------------------------------------------------




## ERNA - VISUALIZE RAW DATA ---------------------------------------------------------------

## Order populations for plotting 
## Order by average size 
ERNA.htByMed <- with(ERNA.cl, reorder(Source, Length_cm_20220915, median, na.rm=TRUE))
ERNA.meds <- ERNA.cl %>% group_by(Source) %>% 
             dplyr::summarise(Height_MD=median(Length_cm_20220915,na.rm=TRUE), GrowthAb_MD=median(GrwthRate_Absolute,na.rm=TRUE),
                   GrowthRe_MD=median(GrwthRate_Relative,na.rm=TRUE), GrowthSp_MD=median(GrwthRate_Specific,na.rm=TRUE))
ERNA.meds <- left_join(ERNA.meds, ERNA.SdZn, by="Source")


## Boxplots of raw data 
par(mfrow=c(2,2))

## Size
ERNA.meds <- ERNA.meds[order(ERNA.meds$Height_MD),] #Order by median 
boxplot(Length_cm_20220915 ~ ERNA.htByMed, data=ERNA.cl,
        xlab=NA, ylab="Height (cm)", cex.lab=1.25,
        cex.axis=0.99, names=ERNA.meds$PopAbbrev, las=2,
        main="FINAL SIZE", cex.main=1.5, col=ERNA.meds$HexCode)

## Growth rate(s) ** Add time interval to growth rate calcs? **
#ERNA.meds <- ERNA.meds[order(ERNA.meds$Growth_MD),]
boxplot(GrwthRate_Relative ~ ERNA.htByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Plant relative growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.5, col=ERNA.meds$HexCode, main="GROWTH RATE", ylim=c(-0.35,3))

#boxplot(GrwthRate_Absolute ~ ERNA.htByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant absolute growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode, main="GROWTH RATE")

#boxplot(GrwthRate_Specific ~ ERNA.htByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant specific growth", cex.lab=1, cex.axis=0.9, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode)
#ERNA.SdZn <- ERNA.SdZn[order(ERNA.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ ERNA.grsByMed, data=ERNA.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=ERNA.SdZn$PopAbbrev,
#        cex.axis=0.7, cex.main=1.5, col=ERNA.SdZn$HexCode)

## Blank plot
plot.new()
legend("center", unique(ERNA.meds$seed_zone[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), col="black",
       pt.bg=unique(ERNA.meds$HexCode[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), cex=1.5, pch=21)
## ---------------------------------------------------


