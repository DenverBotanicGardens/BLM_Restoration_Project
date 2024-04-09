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
PEVI.LatLong <- read.csv(file="AGoebl/Seeds/20230127_PEVI_LatLong.csv", sep=",", header=TRUE, dec=".")
PEVI23 <- read.csv(file="Chatfield/20240129_ChatfieldData2023_PEVI.csv", sep=",", header=TRUE, dec=".")

GravelTemps <- read.csv(file="Chatfield/ClimateData/20240409_GravelTemps2022_ChatfieldCG.csv", sep=",", header=TRUE, dec=".")
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
PEVI.LatLong$Source <- str_replace(PEVI.LatLong$SOURCE.CODE, "3-SOS", "")
PEVI.LatLong$Source <- as.factor(PEVI.LatLong$Source)

## Add Population name abbreviation column, W/in state ordered by decreasing latitude
PEVI.LatLong$PopAbbrev[grepl("PEVI-WY932-29-11", PEVI.LatLong$Source)] = "P.WY.2"   
PEVI.LatLong$PopAbbrev[grepl("PEVI-WY932-31-11", PEVI.LatLong$Source)] = "P.WY.1" 
PEVI.LatLong$PopAbbrev[grepl("PEVI-CO932-227-10", PEVI.LatLong$Source)] = "P.CO.3"   
PEVI.LatLong$PopAbbrev[grepl("PEVI-CO932-217-10", PEVI.LatLong$Source)] = "P.CO.4" 
PEVI.LatLong$PopAbbrev[grepl("PEVI-CO932-214-10", PEVI.LatLong$Source)] = "P.CO.1"   
PEVI.LatLong$PopAbbrev[grepl("PEVI-CO932-232-10", PEVI.LatLong$Source)] = "P.CO.2"   

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

## Make a flowered Yes or No column
PEVI23$FlwrYesNo <- NA 
PEVI23$FlwrYesNo[!is.na(PEVI23$DaysToFlwr)] <- 1                               #If there's a value in Days to Flwr, then enter Yes
## **If plt alive and didn't flwr, enter No **Find better way, as not sure what date(s) to use for 
PEVI23$FlwrYesNo[is.na(PEVI23$DaysToFlwr) & (PEVI23$Survival_20230808==1 | PEVI23$Survival_20230718==1 |
                                             PEVI23$Survival_20230627==1)] <- 0  
nrow(PEVI23[PEVI23$FlwrYesNo==0,])                                             #Num that survived but didn't flower 
## ---------------------------------------------------------------



## CALCULATE RAW FLOWERING RESULTS BY SOURCE AND TREATMENT
PEVI23 %>% group_by(Source) %>% dplyr::summarise(Pheno_Avg=mean(DaysToFlwr,na.rm=TRUE), Pheno_SD=sd(DaysToFlwr,na.rm=TRUE))
PEVI23 %>% group_by(Source, Treatment) %>% dplyr::summarise(Pheno_Avg=mean(DaysToFlwr,na.rm=TRUE), Pheno_SD=sd(DaysToFlwr,na.rm=TRUE))
# Most, but not all, pops in blk tx flowered earlier

PEVI23 %>% group_by(Source) %>% dplyr::summarise(FlwrYesNo_Avg=mean(FlwrYesNo,na.rm=TRUE)) #Or could try sum and/or NUM=n()
PEVI23 %>% group_by(Source, Treatment) %>% dplyr::summarise(FlwrYesNo_Avg=mean(FlwrYesNo,na.rm=TRUE))
# All pops had higher prob of repro in blk tx




## PEVI - TEST FOR SIGNIFICANCE OF TREATMENT EFFECT -------------------------------------
## Flowering phenology
hist(PEVI23$DaysToFlwr)
PEVI.tx.DaystoFlwr <- glmer.nb(DaysToFlwr ~ Source + Treatment + (1|Block), data=PEVI23)
PEVI.pop.DaystoFlwr <- glmer.nb(DaysToFlwr ~ Source + (1|Block), data=PEVI23)
models <- list(PEVI.tx.DaystoFlwr, PEVI.pop.DaystoFlwr)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# Possible support for treatment 
summary(PEVI.tx.DaystoFlwr)
Anova(PEVI.tx.DaystoFlwr)

## Probability of repro
PEVI.tx.pRepro <- glmer(FlwrYesNo ~ Source + Treatment + (1|Block), data=PEVI23, family=binomial (link="logit"))
PEVI.pop.pRepro <- glmer(FlwrYesNo ~ Source + (1|Block), data=PEVI23, family=binomial (link="logit"))
models <- list(PEVI.tx.pRepro, PEVI.pop.pRepro)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# Possible support for treatment 
summary(PEVI.tx.pRepro)
Anova(PEVI.tx.pRepro)
## ---------------------------------------------------------------





## PEVI - VISUALIZE RAW DATA -----------------------------------------------------------------------------

## Assign color based on latitude

## 2023 Flowering data - Plots means as points and SE as error bars
PEVI23.mn <- PEVI23 %>% group_by(Source) %>% dplyr::summarise(DaysToFlwr_MN=mean(DaysToFlwr,na.rm=TRUE),
                                                              DaysToFlwr_SE=calcSE(DaysToFlwr), FlwrYesNo_MN=mean(FlwrYesNo,na.rm=TRUE),
                                                              FlwrYesNo_SE=calcSE(FlwrYesNo))

pop.list <- unique(PEVI23$Source)
PEVI23.mn$PopCol <- NA
PEVI23.mn$PopCol[grepl("PEVI-CO932-217-10", PEVI23.mn$Source)] = "#FEE08B"        
PEVI23.mn$PopCol[grepl("PEVI-CO932-227-10", PEVI23.mn$Source)] = "#FFFF85"        
PEVI23.mn$PopCol[grepl("PEVI-CO932-232-10", PEVI23.mn$Source)] = "#aabe7f"        
PEVI23.mn$PopCol[grepl("PEVI-CO932-214-10", PEVI23.mn$Source)] = "#ABDDA4"     
PEVI23.mn$PopCol[grepl("PEVI-WY932-29-11", PEVI23.mn$Source)] = "#66C2A5"    
PEVI23.mn$PopCol[grepl("PEVI-WY932-31-11", PEVI23.mn$Source)] = "#3288BD"   

#PEVI23.srcCol <- PEVI23 %>% select(Source, PopCol)
PEVI23.mn <- left_join(PEVI23.mn, PEVI.LatLong, by="Source")

PEVI23.mn <- PEVI23.mn[order(PEVI23.mn$Lat),]
plot(c(1:6), PEVI23.mn$FlwrYesNo_MN, col="black", bg=PEVI23.mn$PopCol, pch=21, cex=2.5, ylab="Flowering rate",xlab="Population",
     xaxt="n",cex.lab=1.5, ylim=c(0.45,0.9))
arrows(c(1:6), PEVI23.mn$FlwrYesNo_MN+PEVI23.mn$FlwrYesNo_SE, c(1:6), PEVI23.mn$FlwrYesNo_MN-PEVI23.mn$FlwrYesNo_SE,
       angle=90, length=0, col="grey")
axis(1,at=1:6,labels=PEVI23.mn$PopAbbrev,las=1, cex.axis=1.1, col.ticks="black")


plot(c(1:6), PEVI23.mn$DaysToFlwr_MN, col="black", bg=PEVI23.mn$PopCol, pch=21, cex=2.5, ylab="Days to first flower",xlab="Population",
     xaxt="n",cex.lab=1.5, ylim=c(141,146)) #bg=PEVI23.mn$HexCode
arrows(c(1:6), PEVI23.mn$DaysToFlwr_MN+PEVI23.mn$DaysToFlwr_SE, c(1:6), PEVI23.mn$DaysToFlwr_MN-PEVI23.mn$DaysToFlwr_SE,
       angle=90, length=0, col="grey")
axis(1,at=1:6,labels=PEVI23.mn$PopAbbrev,las=1, cex.axis=1.1, col.ticks="black")

legend("center", legend=PEVI23.mn$PopAbbrev[order(PEVI23.mn$Lat, decreasing=TRUE)], col="black",
       pt.bg=PEVI23.mn$PopCol[order(PEVI23.mn$Lat, decreasing=TRUE)], cex=1.9, pch=21)




## 2023 Flowering data by source and tx - Plots means as points and SE as error bars
PEVI23.mnTx <- PEVI23 %>% group_by(Source, Treatment) %>% dplyr::summarise(DaysToFlwr_MN=mean(DaysToFlwr,na.rm=TRUE),
                                                              DaysToFlwr_SE=calcSE(DaysToFlwr), FlwrYesNo_MN=mean(FlwrYesNo,na.rm=TRUE),
                                                              FlwrYesNo_SE=calcSE(FlwrYesNo))

PEVI23.mnTx$PopCol <- NA
PEVI23.mnTx$PopCol[grepl("PEVI-CO932-217-10", PEVI23.mnTx$Source)] = "#FEE08B"        
PEVI23.mnTx$PopCol[grepl("PEVI-CO932-227-10", PEVI23.mnTx$Source)] = "#FFFF85"        
PEVI23.mnTx$PopCol[grepl("PEVI-CO932-232-10", PEVI23.mnTx$Source)] = "#aabe7f"        
PEVI23.mnTx$PopCol[grepl("PEVI-CO932-214-10", PEVI23.mnTx$Source)] = "#ABDDA4"     
PEVI23.mnTx$PopCol[grepl("PEVI-WY932-29-11", PEVI23.mnTx$Source)] = "#66C2A5"    
PEVI23.mnTx$PopCol[grepl("PEVI-WY932-31-11", PEVI23.mnTx$Source)] = "#3288BD"   

PEVI23.mnTx <- left_join(PEVI23.mnTx, PEVI.LatLong, by="Source")


PEVI23.mnTx <- PEVI23.mnTx[order(PEVI23.mnTx$Lat),]
symbs<-rep(c(22,25),6)
xSpac <- (1:12)-rep(c(0.5,1),6) #Change spacing to cluster population pairs

plot(xSpac, PEVI23.mnTx$FlwrYesNo_MN, col="black", bg=PEVI23.mnTx$PopCol, pch=symbs, cex=2.5, ylab="Flowering rate",xlab="Population",
     xaxt="n",cex.lab=1.5, ylim=c(0.35,0.95))
arrows(xSpac, PEVI23.mnTx$FlwrYesNo_MN+PEVI23.mnTx$FlwrYesNo_SE, xSpac, PEVI23.mnTx$FlwrYesNo_MN-PEVI23.mnTx$FlwrYesNo_SE,
       angle=90, length=0, col="grey")
axis(1,at=xSpac,labels=PEVI23.mnTx$PopAbbrev,las=1, cex.axis=1.1, col.ticks="white")

xSpac <- (1:12)-rep(c(0.5,1.2),6) #Change spacing to cluster population pairs
plot(xSpac, PEVI23.mnTx$DaysToFlwr_MN, col="black", bg=PEVI23.mnTx$PopCol, pch=symbs, cex=2.5, ylab="Days to first flower",xlab="Population",
     xaxt="n",cex.lab=1.5, ylim=c(140,148)) 
arrows(xSpac, PEVI23.mnTx$DaysToFlwr_MN+PEVI23.mnTx$DaysToFlwr_SE, xSpac, PEVI23.mnTx$DaysToFlwr_MN-PEVI23.mnTx$DaysToFlwr_SE,
       angle=90, length=0, col="grey")
axis(1,at=xSpac,labels=PEVI23.mnTx$PopAbbrev,las=1, cex.axis=1.1, col.ticks="white")

plot.new()
legend("center", legend=c("Warm", "Cool"), col="black",
       pt.bg="black", cex=1.9, pch=c(22,25))
## -------------------------------------------------------------------




## MAKE PLOT OF GRAVEL SURFACE TEMPS
GravTmp.mn <- GravelTemps %>% group_by(Tx) %>% dplyr::summarise(TEMP_MN=mean(Reading_DegC,na.rm=TRUE), TEMP_SE=calcSE(Reading_DegC))

#ERNA23.mn <- ERNA23.mn[order(ERNA23.mn$AliveYesNo_MN),]
barXvals<-barplot(GravTmp.mn$TEMP_MN, xlab=NA, ylab="Temperature (°C)", cex.lab=1.4, las=1, 
                  names=c("Black gravel","No gravel","White gravel"), main="Soil surface temperatures", cex.main=1.5, 
                  col=c("black","grey","white"), cex.names=1.1, ylim=c(0,110))
arrows(barXvals, GravTmp.mn$TEMP_MN+GravTmp.mn$TEMP_SE, barXvals, GravTmp.mn$TEMP_MN-GravTmp.mn$TEMP_SE,
       angle=90, length=0, col="black")
## -------------------------------------------------------------------
