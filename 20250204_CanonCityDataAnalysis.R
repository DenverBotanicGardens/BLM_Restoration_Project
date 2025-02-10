## April Goebl
## Script started 2023-11-02
## BLM Restoration project at Denver Botanic Gardens
## Analyze data from Canyon City genetic diversity experiment   




rm(list=ls())
#dev.off()


## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
library(lme4)
library(plotrix)
library(EnvStats)
library(car)
library(effects)
library(tidyverse)
library(AICcmodavg)

calcSE <- function(x){sd(x, na.rm=TRUE)/sqrt(length(x))}
## ------------------------------------------------------------------------------------------------



## SET WORKING DIRECTORY --------------------------------------------------------------------------
#setwd("C:/Users/april.goebl/OneDrive - Denver Botanic Gardens/AGoebl_NonProject_WorkRelated/BLM_Restoration_Project")
#setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/CanyonCity")
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/CanyonCity")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
#dats <- read.csv(file="20231102_CanonCityData_ARFR.csv", sep=",", header=TRUE, dec=".")
dats <- read.csv(file="20250203_CanonCityData_ARFRforAnalysis.csv", sep=",", header=TRUE, dec=".")
## ------------------------------------------------------------------------------------------------



## LOOK AT AND MODIFY STRUCTURE OF DATA -----------------------------------------------------------
str(dats)
dats$Plot.Number <- as.factor(dats$Plot.Number)
#dats$PlotTagNum <- as.factor(dats$PlotTagNum)
dats$Site <- as.factor(dats$Site)
dats$Treatment <- as.factor(dats$Treatment)
#dats$Contents <- as.factor(dats$Contents)
#Change 'Date' to date or parse out year, etc?
## ------------------------------------------------------------------------------------------------



## ADD AND REMOVE FROM DATAFRAME ------------------------------------------------------------------
unique(dats$Contents)
unique(dats$Treatment)
## Remove empty, control, and pilot categories 
#dats <- dats[dats$Contents != "empty",]
#dats <- dats[dats$Contents != "controlAndSave",]  #Or should we keep negative control? 
#dats <- dats[dats$Contents != "Fre12_1200seeds",]
#dats <- dats[dats$Contents != "Fre12_800seeds",]
#dats <- dats[dats$Contents != "422Nav_1200seeds",]
#dats <- dats[dats$Contents != "422Nav_800seeds",]
#dats <- dats[dats$Contents != "WY71_ClearPlot",]
#dats <- dats[dats$Contents != "ERNA_Rio",]
#dats <- dats[dats$Contents != "ERNA_Rou",]
#dats <- dats[dats$Contents != "ERNA_Gun",]



## Add in spatial variable designating position in 10x25 grid
dats$HalfSplit <- rep(c(rep("H1",12), rep("H2",13)),60)
#dats$QuadSplit <- rep(c(rep("H1",12), rep("H2",13)),60) ** split into 4 quads ... 



dats <- dats[dats$Treatment != "skip",]
dats <- dats[dats$Treatment != "pilot",]
#droplevels function


## Data checks
#seedling count should only be numeric
#seedling count on 20240918 should only be 1 or 0, should not be empty
#percent cover 202406 should not be empty, change NAs to 0 (below)
#other checks..?
#If seedling count>0, then percent cover should be >0


## Check if any ARFR seen in negative controls?
unique(dats$ARFR.Seedling.Count[dats$Treatment=="17"])
dats[dats$Treatment=="17" & dats$ARFR.Seedling.Count==1,]
##There is one incidence of a single ARFR seedling being counted in one survey in a single control plot
## Remove negative control rows
dats <- dats[dats$Treatment != "17",] #Negative controls



## Add column for number of seeds planted 
dats$NumSeeds <- 400

## Add column with mix type or single 
dats$SeedMix[grepl("mix_d", dats$Contents)] <- "Mix_B"  #Broad mixes
dats$SeedMix[grepl("mix_s", dats$Contents)] <- "Mix_C"  #Close mixes
dats$SeedMix[grepl("single_", dats$Contents)] <- "Single"
dats$SeedMix <- as.factor(dats$SeedMix)

## Rename mix types
dats$Contents <- str_replace(dats$Contents, "mix_s", "mix_c")
dats$Contents <- str_replace(dats$Contents, "mix_d", "mix_b")
## ------------------------------------------------------------------------------------------------






## PLOT RAW DATA AS BOX PLOTS ----------------------------------------------------------------------------------
#cols <- c(rep("darkseagreen4", 10), rep("red4",10), rep("steelblue4",12))
#par(las=2)
#boxplot(ARFR.Seedling.Count ~ Site + Contents, data=dats, col=cols, ylim=c(0,80),
#        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32,34:35,37:38,40:41,43:44,46:47))

unique(dats$Contents)
unique(dats$SeedMix)
#dats$MixCol[grepl("Single", dats$SeedMix)] = "steelblue4"    
#dats$MixCol[grepl("Mix_C", dats$SeedMix)] = "red4"     
#dats$MixCol[grepl("Mix_B", dats$SeedMix)] = "darkseagreen4"    

## 2023 emergence
dats23 <- dats[grepl("2023", dats$Date),]
dats23$Emrg <- dats23$ARFR.Seedling.Count/dats23$NumSeeds
  
#cols <- c(rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6),
#          rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6))
#boxplot((ARFR.Seedling.Count/NumSeeds) ~ Contents + Site, data=dats23, col=cols, ylim=c(0,0.14),
#        xlab=NA, ylab="Emergence Rate", xaxt="n", cex.lab=1.5)

boxplot(ARFR.Seedling.Count ~ SeedMix + Site, data=dats23)
boxplot(ARFR.Seedling.Count/NumSeeds ~ SeedMix + Site, data=dats23)
boxplot(ARFR.Seedling.Count/NumSeeds ~ SeedMix, data=dats23, ylim=c(0,0.19))
boxplot(ARFR.Seedling.Count/NumSeeds ~ Site, data=dats23, ylim=c(0,0.19), ylab="Emergence Rate",
        cex.lab=1.5)

## Plot emergence rate by site separately
## ** Plot all data points as well as boxes? 
dats23.PC <- dats23[dats23$Site=="PC",]
dats23.EP <- dats23[dats23$Site=="EP",]

boxplot(Emrg ~ SeedMix, data=dats23.PC, ylim=c(0,0.16), col=c("grey40", "grey80", "plum4"), main="Site 1")
boxplot(Emrg ~ SeedMix, data=dats23.EP, ylim=c(0,0.1), col=c("grey40", "grey80", "plum4"), main="Site 2")

par(mar=c(7,5,3,2))
cols <- c(rep("grey40", 5), rep("grey80",5), rep("plum4",6))
boxplot((ARFR.Seedling.Count/NumSeeds) ~ Contents, data=dats23.PC, col=cols, ylim=c(0,0.11),
        main="Site 1", ylab="Emergence Rate", cex.lab=1.5, las=2, xlab=NA)
boxplot((ARFR.Seedling.Count/NumSeeds) ~ Contents, data=dats23.EP, col=cols, ylim=c(0,0.075),
        main="Site 2", ylab="Emergence Rate", cex.lab=1.5, las=2, xlab=NA)



EmrgAvgs.PC <- dats23.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Emrg_MD=median(Emrg,na.rm=TRUE),
                   Emrg_MN=mean(Emrg,na.rm=TRUE),
                   Emrg_SE=calcSE(Emrg))
EmrgAvgs.PC$MixCol[grepl("single", EmrgAvgs.PC$Contents)] = "plum4"    
EmrgAvgs.PC$MixCol[grepl("mix_c", EmrgAvgs.PC$Contents)] = "grey80"     
EmrgAvgs.PC$MixCol[grepl("mix_b", EmrgAvgs.PC$Contents)] = "grey30"    

#EmrgByMed.PC <- with(dats23.PC, reorder(Contents, ARFR.Seedling.Count, median, na.rm=TRUE))
#EmrgByMed <- with(dats23, reorder(Contents, Emrg, median, na.rm=TRUE))
#EmrgAvgs.PC <- EmrgAvgs.PC[order(EmrgAvgs.PC$Emrg_MD),] #Order by median 
#boxplot(Emrg ~ EmrgByMed.PC, data=dats23.PC, ylim=c(0,0.11),
#        xlab=NA, ylab="Emergence rate", cex.lab=1.25, col=EmrgAvgs.PC$MixCol,
#        las=2)
plot.new()
legend("center", c("Single source", "Regional mix", "Broad mix"), 
       col=c("plum4","grey80","grey30"), cex=1.95, pch=19)



## Plots means at each site
EmrgAvgs.EP <- dats23.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Emrg_MD=median(Emrg,na.rm=TRUE),
                   Emrg_MN=mean(Emrg,na.rm=TRUE),
                   Emrg_SE=calcSE(Emrg))
EmrgAvgs.EP$MixCol[grepl("single", EmrgAvgs.EP$Contents)] = "plum4"    
EmrgAvgs.EP$MixCol[grepl("mix_c", EmrgAvgs.EP$Contents)] = "grey80"     
EmrgAvgs.EP$MixCol[grepl("mix_b", EmrgAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,0.07), xaxt="n", ylab="Emergence rate", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(EmrgAvgs.PC$Emrg_MN[ee],EmrgAvgs.EP$Emrg_MN[ee]), col=EmrgAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(EmrgAvgs.PC$Emrg_MN[ee],EmrgAvgs.EP$Emrg_MN[ee]), col=EmrgAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}


#EmrgAvgs <- left_join(EmrgAvgs, dats, by="Contents") 
#EmrgAvgs$PopCol <- NA
#EmrgAvgs$PopCol[grepl("mix_s1", EmrgAvgs$Contents)] = "red4"  
#EmrgAvgs$PopCol[grepl("mix_s2", EmrgAvgs$Contents)] = "red4"     
#EmrgAvgs$PopCol[grepl("mix_s3", EmrgAvgs$Contents)] = "red4"        
#EmrgAvgs$PopCol[grepl("mix_s4", EmrgAvgs$Contents)] = "red4"        
#EmrgAvgs$PopCol[grepl("mix_s5", EmrgAvgs$Contents)] = "red4"        
#EmrgAvgs$PopCol[grepl("mix_d1", EmrgAvgs$Contents)] = "darkseagreen4"         
#EmrgAvgs$PopCol[grepl("mix_d2", EmrgAvgs$Contents)] = "darkseagreen4"         
#EmrgAvgs$PopCol[grepl("mix_d3", EmrgAvgs$Contents)] = "darkseagreen4"         
#EmrgAvgs$PopCol[grepl("mix_d4", EmrgAvgs$Contents)] = "darkseagreen4"         
#EmrgAvgs$PopCol[grepl("mix_d5", EmrgAvgs$Contents)] = "darkseagreen4" 
#EmrgAvgs$PopCol[grepl("single_UT", EmrgAvgs$Contents)] = "steelblue4" 
#EmrgAvgs$PopCol[grepl("single_294", EmrgAvgs$Contents)] = "steelblue4" 
#EmrgAvgs$PopCol[grepl("single_NM", EmrgAvgs$Contents)] = "steelblue4" 
#EmrgAvgs$PopCol[grepl("single_316Jeff", EmrgAvgs$Contents)] = "steelblue4" 
#EmrgAvgs$PopCol[grepl("single_314Jeff", EmrgAvgs$Contents)] = "steelblue4" 
#EmrgAvgs$PopCol[grepl("single_LasAn", EmrgAvgs$Contents)] = "steelblue4" 

#plot(1:16, EmrgAvgs$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)


#EmrgAvgs_PC <- PC %>% group_by(Contents) %>% 
#  dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
#                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
#                   Emrg_SE=calcSE(ARFR.Seedling.Count))

#plot(1:16, EmrgAvgs_PC$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)

#EmrgAvgs_EP <- EP %>% group_by(Contents) %>% 
#  dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
#                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
#                   Emrg_SE=calcSE(ARFR.Seedling.Count))

#plot(1:16, EmrgAvgs_EP$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)




## 2024 
dats24 <- dats[grepl("2024", dats$Date),]
dats2406 <- dats24[dats24$Date!="9/18/2024 18:00",]

## Change any NAs (empty entries) to zero in percent veg cover
dats2406$ARFR.Percent.Cover[is.na(dats2406$ARFR.Percent.Cover)] <- 0


## Survival rates (2024 counts / 2023 counts)
plot(dats23$ARFR.Seedling.Count, dats2406$ARFR.Seedling.Count)

## Change zeros to NAs in 2023 counts just for surv estimates
dats23temp <- dats23
dats23temp$ARFR.Seedling.Count[dats23temp$ARFR.Seedling.Count==0] <- NA
dats2406$Surv <- dats2406$ARFR.Seedling.Count/dats23temp$ARFR.Seedling.Count

## Plot survival by site separately
dats2406.PC <- dats2406[dats2406$Site=="PC",]
dats2406.EP <- dats2406[dats2406$Site=="EP",]

par(mar=c(7,5,3,2))
cols <- c(rep("grey40", 5), rep("grey80",5), rep("plum4",6))
boxplot(Surv ~ Contents, data=dats2406.PC, col=cols, ylim=c(0,1), cex.main=1.5,
        main="Site 1", ylab="Survival rate", xlab=NA, cex.lab=1.5, las=2)
boxplot(Surv ~ Contents, data=dats2406.EP, col=cols, ylim=c(0,2), cex.main=1.5,
        main="Site 2", ylab="Survival rate", xlab=NA, cex.lab=1.5, las=2)


SurvAvgs.PC <- dats2406.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Surv_MD=median(Surv,na.rm=TRUE),
                   Surv_MN=mean(Surv,na.rm=TRUE),
                   Surv_SE=calcSE(Surv),
                   Surv_SUM=sum(Surv,na.rm=TRUE))
SurvAvgs.PC$MixCol[grepl("single", SurvAvgs.PC$Contents)] = "plum4"    
SurvAvgs.PC$MixCol[grepl("mix_c", SurvAvgs.PC$Contents)] = "grey80"     
SurvAvgs.PC$MixCol[grepl("mix_b", SurvAvgs.PC$Contents)] = "grey30"    

#SurvByMed.PC <- with(dats2406.PC, reorder(Contents, Surv, median, na.rm=TRUE))
#SurvAvgs.PC <- SurvAvgs.PC[order(SurvAvgs.PC$Surv_MD),] #Order by median 
#par(mar=c(7,5,3,2))
#boxplot(Surv ~ SurvByMed.PC, data=dats2406.PC, ylim=c(0,1),
#        xlab=NA, ylab="Survival rate", cex.lab=1.25, col=SurvAvgs.PC$MixCol,
#        las=2)

#dats2406.PCmix <- dats2406[dats2406$SeedMix!="Single",]
boxplot(Surv ~ SeedMix, data=dats2406.PC, ylim=c(0,1), col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Survival rate")
#boxplot(Surv ~ SeedMix, data=dats2406.EP, ylim=c(0,1.7), col=c("grey40", "grey80", "plum4"), main="Site 2",
#        ylab="Survival rate")


## Plots means at each site
SurvAvgs.EP <- dats2406.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Surv_MD=median(Surv,na.rm=TRUE),
                   Surv_MN=mean(Surv,na.rm=TRUE),
                   Surv_SE=calcSE(Surv))
SurvAvgs.EP$MixCol[grepl("single", SurvAvgs.EP$Contents)] = "plum4"    
SurvAvgs.EP$MixCol[grepl("mix_c", SurvAvgs.EP$Contents)] = "grey80"     
SurvAvgs.EP$MixCol[grepl("mix_b", SurvAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,0.6), xaxt="n", ylab="Survival rate", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(SurvAvgs.PC$Surv_MN[ee],SurvAvgs.EP$Surv_MN[ee]), col=SurvAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(SurvAvgs.PC$Surv_MN[ee],SurvAvgs.EP$Surv_MN[ee]), col=SurvAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}



## Counts
#cols <- c(rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6),
#          rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6))
#boxplot((ARFR.Seedling.Count) ~ Contents + Site, data=dats2406, col=cols, ylim=c(0,20),
#        xlab=NA, ylab="Number of plants per plot", xaxt="n", cex.lab=1.5)

boxplot(ARFR.Seedling.Count ~ SeedMix + Site, data=dats2406)
boxplot(ARFR.Seedling.Count ~ SeedMix, data=dats2406, ylim=c(0,26))
boxplot(ARFR.Seedling.Count ~ Site, data=dats2406, ylim=c(0,26), ylab="Number of plants/plot",
        cex.lab=1.5)
boxplot(ARFR.Seedling.Count ~ SeedMix, data=dats2406.PC, ylim=c(0,30), col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Counts")
boxplot(ARFR.Seedling.Count ~ SeedMix, data=dats2406.EP, ylim=c(0,14), col=c("grey40", "grey80", "plum4"), main="Site 2",
        ylab="Counts")



## Plot counts by site separately
cols <- c(rep("grey40", 5), rep("grey80",5), rep("plum4",6))
boxplot((ARFR.Seedling.Count) ~ Contents, data=dats2406.PC, col=cols, ylim=c(0,26), cex.main=1.5,
        main="Site 1", ylab="Number of plants/plot", xlab=NA, cex.lab=1.5, las=2)
boxplot((ARFR.Seedling.Count) ~ Contents, data=dats2406.EP, col=cols, ylim=c(0,7), cex.main=1.5,
        main="Site 2", ylab="Number of plants/plot", xlab=NA, cex.lab=1.5, las=2) #xaxt="n",


#CountByMed.PC <- with(dats2406.PC, reorder(Contents, ARFR.Seedling.Count, median, na.rm=TRUE))

CountAvgs.PC <- dats2406.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Count_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_SE=calcSE(ARFR.Seedling.Count),
                   Count_SUM=sum(ARFR.Seedling.Count,na.rm=TRUE))
CountAvgs.PC$MixCol[grepl("single", CountAvgs.PC$Contents)] = "plum4"    
CountAvgs.PC$MixCol[grepl("mix_c", CountAvgs.PC$Contents)] = "grey80"     
CountAvgs.PC$MixCol[grepl("mix_b", CountAvgs.PC$Contents)] = "grey30"    

CountAvgs.EP <- dats2406.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Count_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_SE=calcSE(ARFR.Seedling.Count),
                   Count_SUM=sum(ARFR.Seedling.Count,na.rm=TRUE))
CountAvgs.EP$MixCol[grepl("single", CountAvgs.EP$Contents)] = "plum4"    
CountAvgs.EP$MixCol[grepl("mix_c", CountAvgs.EP$Contents)] = "grey80"     
CountAvgs.EP$MixCol[grepl("mix_b", CountAvgs.EP$Contents)] = "grey30"    

#CountAvgs.PC <- CountAvgs.PC[order(CountAvgs.PC$Count_MD),] #Order by median 
#par(mar=c(7,5,3,2))
#boxplot(ARFR.Seedling.Count ~ CountByMed.PC, data=dats2406.PC, ylim=c(0,19),
#        xlab=NA, ylab="Number of plants/plot", cex.lab=1.25, col=CountAvgs.PC$MixCol,
#        las=2)


## Plots means at each site
CountAvgs.EP <- dats2406.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Count_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Count_SE=calcSE(ARFR.Seedling.Count))
CountAvgs.EP$MixCol[grepl("single", CountAvgs.EP$Contents)] = "plum4"    
CountAvgs.EP$MixCol[grepl("mix_c", CountAvgs.EP$Contents)] = "grey80"     
CountAvgs.EP$MixCol[grepl("mix_b", CountAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,9), xaxt="n", ylab="Count", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(CountAvgs.PC$Count_MN[ee],CountAvgs.EP$Count_MN[ee]), col=CountAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(CountAvgs.PC$Count_MN[ee],CountAvgs.EP$Count_MN[ee]), col=CountAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}


## Total counts across all plots
CountAvgs.PC <- CountAvgs.PC[order(CountAvgs.PC$Count_SUM),] #Order by total count
barplot(CountAvgs.PC$Count_SUM, xlab=NA, ylab="Total number of plants", cex.lab=1.25, las=2, 
        names=CountAvgs.PC$Contents, main="Site 1", cex.main=1.5, col=CountAvgs.PC$MixCol, cex.names=1)

CountAvgs.EP <- CountAvgs.EP[order(CountAvgs.EP$Count_SUM),] #Order by total count
barplot(CountAvgs.EP$Count_SUM, xlab=NA, ylab="Total number of plants", cex.lab=1.25, las=2, 
        names=CountAvgs.EP$Contents, main="Site 2", cex.main=1.5, col=CountAvgs.EP$MixCol, cex.names=1)
## ** Get avgs per mix type? ** Maybe not so meaningful but could be worth visualizing **




## Percent vegetative cover
plot(dats2406$ARFR.Seedling.Count, dats2406$ARFR.Percent.Cover) #How correlated are counts and cover?

par(mar=c(7,5,3,2))
boxplot(ARFR.Percent.Cover ~ Contents, data=dats2406.PC, col=cols, ylim=c(0,56), cex.main=1.5,
        main="Site 1", ylab="Percent A. frigida cover/plot", xlab=NA, cex.lab=1.5, las=2)
boxplot(ARFR.Percent.Cover ~ Contents, data=dats2406.EP, col=cols, ylim=c(0,20), cex.main=1.5,
        main="Site 2", ylab="Percent A. frigida cover/plot", xlab=NA, cex.lab=1.5, las=2) 


CovrAvgs.PC <- dats2406.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Covr_MD=median(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_MN=mean(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_SE=calcSE(ARFR.Percent.Cover))
CovrAvgs.PC$MixCol[grepl("single", CovrAvgs.PC$Contents)] = "plum4"    
CovrAvgs.PC$MixCol[grepl("mix_c", CovrAvgs.PC$Contents)] = "grey80"     
CovrAvgs.PC$MixCol[grepl("mix_b", CovrAvgs.PC$Contents)] = "grey30"    

#CovrByMed.PC <- with(dats2406.PC, reorder(Contents, ARFR.Percent.Cover, median, na.rm=TRUE))
#CovrAvgs.PC <- CovrAvgs.PC[order(CovrAvgs.PC$Covr_MD),] #Order by median 
#boxplot(ARFR.Percent.Cover ~ CovrByMed.PC, data=dats2406.PC, ylim=c(0,59),
#        xlab=NA, ylab="Percent A. frigida cover/plot", cex.lab=1.25, col=CovrAvgs.PC$MixCol,
#        las=2)

boxplot(ARFR.Percent.Cover ~ SeedMix + Site, data=dats2406)
boxplot(ARFR.Percent.Cover ~ Site, data=dats2406, ylim=c(0,56), ylab="Percent A. frigida cover/plot",
        cex.lab=1.5)
boxplot(ARFR.Percent.Cover ~ SeedMix, data=dats2406.PC, ylim=c(0,60), col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Percent A. frigida cover")
boxplot(ARFR.Percent.Cover ~ SeedMix, data=dats2406.EP, ylim=c(0,10), col=c("grey40", "grey80", "plum4"), main="Site 2",
        ylab="Percent A. frigida cover")

 
 
## Plots means at each site
CovrAvgs.EP <- dats2406.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Covr_MD=median(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_MN=mean(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_SE=calcSE(ARFR.Percent.Cover))
CovrAvgs.EP$MixCol[grepl("single", CovrAvgs.EP$Contents)] = "plum4"    
CovrAvgs.EP$MixCol[grepl("mix_c", CovrAvgs.EP$Contents)] = "grey80"     
CovrAvgs.EP$MixCol[grepl("mix_b", CovrAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,30), xaxt="n", ylab="Vegetative Cover", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(CovrAvgs.PC$Covr_MN[ee],CovrAvgs.EP$Covr_MN[ee]), col=CovrAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(CovrAvgs.PC$Covr_MN[ee],CovrAvgs.EP$Covr_MN[ee]), col=CovrAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}
  
  
  


## Percent reproductive cover
dats2409 <- dats24[dats24$Date=="9/18/2024 18:00",]

## Repro should be NA, not zero, if seedling=0
dats2409$ARFR.Reproductive.Percent.Cover[dats2409$ARFR.Seedling.Count==0] <- NA

plot(dats2406$ARFR.Percent.Cover,dats2409$ARFR.Reproductive.Percent.Cover) #How correlated are veg and repro cover?


## Make combined performance measure
## With reproduction
dats2409$CombPerf <- dats23$Emrg * dats2406$Surv * dats2409$ARFR.Reproductive.Percent.Cover
## With veg cover
dats2409$CombPerfCovr <- dats23$Emrg * dats2406$Surv * dats2406$ARFR.Percent.Cover

  
  
## Plot by site separately
dats2409.PC <- dats2409[dats2409$Site=="PC",]
dats2409.EP <- dats2409[dats2409$Site=="EP",]

boxplot(ARFR.Reproductive.Percent.Cover ~ Contents, data=dats2409.PC, col=cols, ylim=c(0,86), cex.main=1.5,
        main="Site 1", ylab="Percent A. frigida reproductive cover", xlab=NA, cex.lab=1.5, las=2)
## Add dummy values so all seed mix types appear in plot
dats2409.EP$ARFR.Reproductive.Percent.Cover[dats2409.EP$Plot.Number=='8'] <- 0
boxplot(ARFR.Reproductive.Percent.Cover ~ Contents, data=dats2409.EP, col=cols, ylim=c(0,16), cex.main=1.5,
        main="Site 2", ylab="Percent A. frigida reproductive cover", xlab=NA, cex.lab=1.5, las=2) 



ReproAvgs.PC <- dats2409.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Repro_MD=median(ARFR.Reproductive.Percent.Cover,na.rm=TRUE),
                   Repro_MN=mean(ARFR.Reproductive.Percent.Cover,na.rm=TRUE),
                   Repro_SE=calcSE(ARFR.Reproductive.Percent.Cover))
ReproAvgs.PC$MixCol[grepl("single", ReproAvgs.PC$Contents)] = "plum4"    
ReproAvgs.PC$MixCol[grepl("mix_c", ReproAvgs.PC$Contents)] = "grey80"     
ReproAvgs.PC$MixCol[grepl("mix_b", ReproAvgs.PC$Contents)] = "grey30"    

#ReproByMed.PC <- with(dats2409.PC, reorder(Contents, ARFR.Reproductive.Percent.Cover, median, na.rm=TRUE))
#ReproAvgs.PC <- ReproAvgs.PC[order(ReproAvgs.PC$Repro_MD),] #Order by median 
#boxplot(ARFR.Reproductive.Percent.Cover ~ ReproByMed.PC, data=dats2409.PC, ylim=c(0,86),
#        xlab=NA, ylab="Percent A. frigida reproductive cover", cex.lab=1.25, col=ReproAvgs.PC$MixCol,
#        las=2)

boxplot(ARFR.Reproductive.Percent.Cover ~ SeedMix + Site, data=dats2409)
boxplot(ARFR.Reproductive.Percent.Cover ~ SeedMix, data=dats2409.PC, col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Percent A. frigida reproductive cover")
boxplot(ARFR.Reproductive.Percent.Cover ~ SeedMix, data=dats2409.EP, ylim=c(0,9),
        ylab="Percent A. frigida reproductive cover", col=c("grey40", "grey80", "plum4"), main="Site 2")
boxplot(ARFR.Reproductive.Percent.Cover ~ Site, data=dats2409, ylim=c(0,86), 
        ylab="Percent A. frigida reproductive cover", cex.lab=1.5)



## Plot means at each site
ReproAvgs.EP <- dats2409.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Repro_MD=median(ARFR.Reproductive.Percent.Cover,na.rm=TRUE),
                   Repro_MN=mean(ARFR.Reproductive.Percent.Cover,na.rm=TRUE),
                   Repro_SE=calcSE(ARFR.Reproductive.Percent.Cover))
ReproAvgs.EP$MixCol[grepl("single", ReproAvgs.EP$Contents)] = "plum4"    
ReproAvgs.EP$MixCol[grepl("mix_c", ReproAvgs.EP$Contents)] = "grey80"     
ReproAvgs.EP$MixCol[grepl("mix_b", ReproAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,45), xaxt="n", ylab="Reproductive Cover", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(ReproAvgs.PC$Repro_MN[ee],ReproAvgs.EP$Repro_MN[ee]), col=ReproAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(ReproAvgs.PC$Repro_MN[ee],ReproAvgs.EP$Repro_MN[ee]), col=ReproAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}







## Combined Performance
## with Repro
boxplot(CombPerf ~ Contents, data=dats2409.PC, col=cols, ylim=c(0,4), cex.main=1.5,
        main="Site 1", ylab="Combined performance", xlab=NA, cex.lab=1.5, las=2)
## Add dummy values so all seed mix types appear in plot
dats2409.EP$CombPerf[dats2409.EP$Plot.Number=='8'] <- 0
boxplot(CombPerf ~ Contents, data=dats2409.EP, col=cols, ylim=c(0,0.12), cex.main=1.5,
        main="Site 2", ylab="Combined performance", xlab=NA, cex.lab=1.5, las=2) 


#CombByMed.PC <- with(dats2409.PC, reorder(Contents, CombPerf, median, na.rm=TRUE))

CombAvgs.PC <- dats2409.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Comb_MD=median(CombPerf,na.rm=TRUE),
                   Comb_MN=mean(CombPerf,na.rm=TRUE),
                   Comb_SE=calcSE(CombPerf))
CombAvgs.PC$MixCol[grepl("single", CombAvgs.PC$Contents)] = "plum4"    
CombAvgs.PC$MixCol[grepl("mix_c", CombAvgs.PC$Contents)] = "grey80"     
CombAvgs.PC$MixCol[grepl("mix_b", CombAvgs.PC$Contents)] = "grey30"    

#CombAvgs.PC <- CombAvgs.PC[order(CombAvgs.PC$Comb_MD),] #Order by median 
#boxplot(CombPerf ~ CombByMed.PC, data=dats2409.PC, ylim=c(0,8),
#        xlab=NA, ylab="Combined Performance", cex.lab=1.25, col=CombAvgs.PC$MixCol,
#        las=2)

boxplot(CombPerf ~ SeedMix + Site, data=dats2409)
boxplot(CombPerf ~ SeedMix, data=dats2409.PC, col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Combined performance")
boxplot(CombPerf ~ SeedMix, data=dats2409.EP, ylim=c(0,0.22),
        ylab="Combined performance", col=c("grey40", "grey80", "plum4"), main="Site 2")
boxplot(CombPerf ~ Site, data=dats2409, ylim=c(0,6), 
        ylab="Combined performance", cex.lab=1.5)



## Plot means at each site
CombAvgs.EP <- dats2409.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Comb_MD=median(CombPerf,na.rm=TRUE),
                   Comb_MN=mean(CombPerf,na.rm=TRUE),
                   Comb_SE=calcSE(CombPerf))
CombAvgs.EP$MixCol[grepl("single", CombAvgs.EP$Contents)] = "plum4"    
CombAvgs.EP$MixCol[grepl("mix_c", CombAvgs.EP$Contents)] = "grey80"     
CombAvgs.EP$MixCol[grepl("mix_b", CombAvgs.EP$Contents)] = "grey30"  

plot(1:2, 0:1, type="n",ylim=c(0,1.3), xaxt="n", ylab="Combined performance", xlab="Site", cex.lab=1.5)
for (ee in 1:length(unique(dats$Contents))) { 
  points(1:2, c(CombAvgs.PC$Comb_MN[ee],CombAvgs.EP$Comb_MN[ee]), col=CombAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  lines(1:2, c(CombAvgs.PC$Comb_MN[ee],CombAvgs.EP$Comb_MN[ee]), col=CombAvgs.PC$MixCol[ee], pch=19, cex=1.5)
  
}



## with veg cover
boxplot(CombPerfCovr ~ Contents, data=dats2409.PC, col=cols, ylim=c(0,3), cex.main=1.5,
        main="Site 1", ylab="Combined performance", xlab=NA, cex.lab=1.5, las=2)
## Add dummy values so all seed mix types appear in plot
#dats2409.EP$CombPerfCovr[dats2409.EP$Plot.Number=='8'] <- 0
boxplot(CombPerfCovr ~ Contents, data=dats2409.EP, col=cols, ylim=c(0,0.3), cex.main=1.5,
        main="Site 2", ylab="Combined performance", xlab=NA, cex.lab=1.5, las=2) 


#CombByMed.PC <- with(dats2409.PC, reorder(Contents, CombPerfCovr, median, na.rm=TRUE))

CombAvgs.PC <- dats2409.PC %>% group_by(Contents) %>% 
  dplyr::summarise(Comb_MD=median(CombPerfCovr,na.rm=TRUE),
                   Comb_MN=mean(CombPerfCovr,na.rm=TRUE),
                   Comb_SE=calcSE(CombPerfCovr))
CombAvgs.PC$MixCol[grepl("single", CombAvgs.PC$Contents)] = "plum4"    
CombAvgs.PC$MixCol[grepl("mix_c", CombAvgs.PC$Contents)] = "grey80"     
CombAvgs.PC$MixCol[grepl("mix_b", CombAvgs.PC$Contents)] = "grey30"    

#CombAvgs.PC <- CombAvgs.PC[order(CombAvgs.PC$Comb_MD),] #Order by median 
#boxplot(CombPerfCovr ~ CombByMed.PC, data=dats2409.PC, ylim=c(0,8),
#        xlab=NA, ylab="Combined Performance", cex.lab=1.25, col=CombAvgs.PC$MixCol,
#        las=2)

boxplot(CombPerfCovr ~ SeedMix + Site, data=dats2409)
boxplot(CombPerfCovr ~ SeedMix, data=dats2409.PC, col=c("grey40", "grey80", "plum4"), main="Site 1",
        ylab="Combined performance")
boxplot(CombPerfCovr ~ SeedMix, data=dats2409.EP, ylim=c(0,0.22),
        ylab="Combined performance", col=c("grey40", "grey80", "plum4"), main="Site 2")
boxplot(CombPerfCovr ~ Site, data=dats2409, ylim=c(0,6), 
        ylab="Combined performance", cex.lab=1.5)
## -----------------------------------------------------------------------------------------------------------------






## ** Dive into several mixes to ask if observed is above or below expected based on constituent sources? **
## -------------------------------------------------------------------------------------------
## TRY PLOTTING DEVIATION FROM EXPECTED ------------------------------------------------------
#Emrg_rates <- dats23 %>% group_by(Site, Contents) %>% reframe(RATE=ARFR.Seedling.Count/NumSeeds)

## Subset by mix/ single type, order lowest to highest emrg rate
# Site 1 (PC) only 
#Emrg_mixC1 <- subset(Emrg_rates, Contents=="mix_c1" & Site=="PC")
#Emrg_mixC1 <- Emrg_mixC1[order(Emrg_mixC1$RATE),]
#Emrg_mixC2 <- subset(Emrg_rates, Contents=="mix_c2" & Site=="PC")
#Emrg_mixC2 <- Emrg_mixC2[order(Emrg_mixC2$RATE),]
#Emrg_mixC3 <- subset(Emrg_rates, Contents=="mix_c3" & Site=="PC")
#Emrg_mixC3 <- Emrg_mixC3[order(Emrg_mixC3$RATE),]
#Emrg_mixC4 <- subset(Emrg_rates, Contents=="mix_c4" & Site=="PC")
#Emrg_mixC4 <- Emrg_mixC4[order(Emrg_mixC4$RATE),]
#Emrg_mixC5 <- subset(Emrg_rates, Contents=="mix_c5" & Site=="PC")
#Emrg_mixC5 <- Emrg_mixC5[order(Emrg_mixC5$RATE),]

#Emrg_mixB1 <- subset(Emrg_rates, Contents=="mix_b1" & Site=="PC")
#Emrg_mixB1 <- Emrg_mixB1[order(Emrg_mixB1$RATE),]
#Emrg_mixB2 <- subset(Emrg_rates, Contents=="mix_b2" & Site=="PC")
#Emrg_mixB2 <- Emrg_mixB2[order(Emrg_mixB2$RATE),]
#Emrg_mixB3 <- subset(Emrg_rates, Contents=="mix_b3" & Site=="PC")
#Emrg_mixB3 <- Emrg_mixB3[order(Emrg_mixB3$RATE),]
#Emrg_mixB4 <- subset(Emrg_rates, Contents=="mix_b4" & Site=="PC")
#Emrg_mixB4 <- Emrg_mixB4[order(Emrg_mixB4$RATE),]
#Emrg_mixB5 <- subset(Emrg_rates, Contents=="mix_b5" & Site=="PC")
#Emrg_mixB5 <- Emrg_mixB5[order(Emrg_mixB5$RATE),]

#Emrg_UT <- subset(Emrg_rates, Contents=="single_UT" & Site=="PC")
#Emrg_UT <- Emrg_UT[order(Emrg_UT$RATE),]
#Emrg_294 <- subset(Emrg_rates, Contents=="single_294" & Site=="PC")
#Emrg_294 <- Emrg_294[order(Emrg_294$RATE),]
#Emrg_NM <- subset(Emrg_rates, Contents=="single_NM" & Site=="PC")
#Emrg_NM <- Emrg_NM[order(Emrg_NM$RATE),]
#Emrg_316Jeff <- subset(Emrg_rates, Contents=="single_316Jeff" & Site=="PC")
#Emrg_316Jeff <- Emrg_316Jeff[order(Emrg_316Jeff$RATE),]
#Emrg_314Jeff <- subset(Emrg_rates, Contents=="single_314Jeff" & Site=="PC")
#Emrg_314Jeff <- Emrg_314Jeff[order(Emrg_314Jeff$RATE),]
#Emrg_LasAn <- subset(Emrg_rates, Contents=="single_LasAn" & Site=="PC")
#Emrg_LasAn <- Emrg_LasAn[order(Emrg_LasAn$RATE),]

## List what components are in each mix
#mix_c1: LasAn, 314Jeff, 316Jeff
#mix_c2: NM, 294, 316Jeff
#mix_c3: UT, 294, 314Jeff
#mix_c4: 294, 314Jeff, 316Jeff
#mix_c5: LasAn, 294, 316Jeff

EmrgAvgs.PC

## Calculate 'expected' count for each mix based on additive of components
eEmrg_c1 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_LasAn"] + 
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"] +
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3
eEmrg_c2 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_NM"] + 
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3
eEmrg_c3 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_UT"] + 
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"])/3
eEmrg_c4 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"] + 
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"])/3
eEmrg_c5 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_LasAn"] + 
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3

EmrgByPlot.PC <- dats23.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Emrg=mean(Emrg,na.rm=TRUE))

## Calculate deviation from expected for each plot of each close mix
eEmrg_c1Dev <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c1"] - eEmrg_c1
eEmrg_c2Dev <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c2"] - eEmrg_c2
eEmrg_c3Dev <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c3"] - eEmrg_c3
eEmrg_c4Dev <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c4"] - eEmrg_c4
eEmrg_c5Dev <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c5"] - eEmrg_c5

## Combine diff mixes together and add values for x-axis for plotting
eEmrgDevs <- c(eEmrg_c1Dev, eEmrg_c2Dev, eEmrg_c3Dev, eEmrg_c4Dev, eEmrg_c5Dev)
xAx <- c(rep(1,length(eEmrg_c1Dev)), rep(2,length(eEmrg_c2Dev)), rep(3,length(eEmrg_c3Dev)),
         rep(4,length(eEmrg_c4Dev)), rep(5,length(eEmrg_c5Dev)))

eEmrgDevsForPlot <- as.data.frame(cbind(eEmrgDevs, xAx))

par(pty="s")
plot(eEmrgDevsForPlot$xAx, eEmrgDevsForPlot$eEmrgDevs, pch=16, ylim=c(-0.06,0.06),
     cex=0.9, ylab="Deviation from expected\nemergence rate", col="grey45",
     xlab="Regional Mix")
abline(h=0)



## Plot means 
## Calculate deviation from expected for each plot of each regional mix
eMnEmrg_c1Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c1"] - eEmrg_c1
eMnEmrg_c2Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c2"]  - eEmrg_c2
eMnEmrg_c3Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c3"]  - eEmrg_c3
eMnEmrg_c4Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c4"]  - eEmrg_c4
eMnEmrg_c5Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c5"]  - eEmrg_c5

## Combine diff mixes together and add values for x-axis for plotting
eMnEmrgDevs <- c(eMnEmrg_c1Dev, eMnEmrg_c2Dev, eMnEmrg_c3Dev, eMnEmrg_c4Dev, eMnEmrg_c5Dev)
xAx <- c(1:5)

eMnEmrgDevsForPlot <- as.data.frame(cbind(eMnEmrgDevs, xAx))

points(eMnEmrgDevsForPlot$xAx, eMnEmrgDevsForPlot$eMnEmrgDevs, pch=16,
       cex=2, col="black") 
## ----------------------------




## BOOTSTRAP TO GENERATE EXPECTED MIX PERFORMANCE
reps <- 100
smpl.emrg <- as.data.frame(matrix(NA, reps, (5*3)))
colnames(smpl.emrg) <- c("c1.1","c1.2","c1.3","c2.1","c2.2","c2.3","c3.1","c3.2","c3.3",
                         "c4.1","c4.2","c4.3","c5.1","c5.2","c5.3")


unique(dats$Contents)

## Regional mix 1
#mix_c1: LasAn, 314Jeff, 316Jeff
emrg_c1.1 <- dats23.PC$Emrg[dats23.PC$Contents=="single_LasAn"]
emrg_c1.2 <- dats23.PC$Emrg[dats23.PC$Contents=="single_314Jeff"]
emrg_c1.3 <- dats23.PC$Emrg[dats23.PC$Contents=="single_316Jeff"]

smpl.emrg$c1.1 <- sample(emrg_c1.1[!is.na(emrg_c1.1)], reps, replace=TRUE)
smpl.emrg$c1.2 <- sample(emrg_c1.2[!is.na(emrg_c1.2)], reps, replace=TRUE)
smpl.emrg$c1.3 <- sample(emrg_c1.3[!is.na(emrg_c1.3)], reps, replace=TRUE)

emrg_c1.expBS <- smpl.emrg$c1.1/3 + smpl.emrg$c1.2/3 + smpl.emrg$c1.3/3
hist(emrg_c1.expBS, breaks=100)
emrg_c1.expBSmn <- mean(emrg_c1.expBS)

emrg_c1.devBS <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c1"] - emrg_c1.expBSmn

t.test(EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c1"], emrg_c1.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 2
#mix_c2: NM, 294, 316Jeff
emrg_c2.1 <- dats23.PC$Emrg[dats23.PC$Contents=="single_NM"]
emrg_c2.2 <- dats23.PC$Emrg[dats23.PC$Contents=="single_294"]
emrg_c2.3 <- dats23.PC$Emrg[dats23.PC$Contents=="single_316Jeff"]

smpl.emrg$c2.1 <- sample(emrg_c2.1[!is.na(emrg_c2.1)], reps, replace=TRUE)
smpl.emrg$c2.2 <- sample(emrg_c2.2[!is.na(emrg_c2.2)], reps, replace=TRUE)
smpl.emrg$c2.3 <- sample(emrg_c2.3[!is.na(emrg_c2.3)], reps, replace=TRUE)

emrg_c2.expBS <- smpl.emrg$c2.1/3 + smpl.emrg$c2.2/3 + smpl.emrg$c2.3/3
hist(emrg_c2.expBS, breaks=50)
emrg_c2.expBSmn <- mean(emrg_c2.expBS)

emrg_c2.devBS <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c2"] - emrg_c2.expBSmn

t.test(EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c2"], emrg_c2.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 3
#mix_c3: UT, 294, 314Jeff
emrg_c3.1 <- dats23.PC$Emrg[dats23.PC$Contents=="single_UT"]
emrg_c3.2 <- dats23.PC$Emrg[dats23.PC$Contents=="single_294"]
emrg_c3.3 <- dats23.PC$Emrg[dats23.PC$Contents=="single_314Jeff"]

smpl.emrg$c3.1 <- sample(emrg_c3.1[!is.na(emrg_c3.1)], reps, replace=TRUE)
smpl.emrg$c3.2 <- sample(emrg_c3.2[!is.na(emrg_c3.2)], reps, replace=TRUE)
smpl.emrg$c3.3 <- sample(emrg_c3.3[!is.na(emrg_c3.3)], reps, replace=TRUE)

emrg_c3.expBS <- smpl.emrg$c3.1/3 + smpl.emrg$c3.2/3 + smpl.emrg$c3.3/3
hist(emrg_c3.expBS, breaks=50)
emrg_c3.expBSmn <- mean(emrg_c3.expBS)

emrg_c3.devBS <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c3"] - emrg_c3.expBSmn

t.test(EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c3"], emrg_c3.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 4
#mix_c4: 294, 314Jeff, 316Jeff
emrg_c4.1 <- dats23.PC$Emrg[dats23.PC$Contents=="single_294"]
emrg_c4.2 <- dats23.PC$Emrg[dats23.PC$Contents=="single_314Jeff"]
emrg_c4.3 <- dats23.PC$Emrg[dats23.PC$Contents=="single_316Jeff"]

smpl.emrg$c4.1 <- sample(emrg_c4.1[!is.na(emrg_c4.1)], reps, replace=TRUE)
smpl.emrg$c4.2 <- sample(emrg_c4.2[!is.na(emrg_c4.2)], reps, replace=TRUE)
smpl.emrg$c4.3 <- sample(emrg_c4.3[!is.na(emrg_c4.3)], reps, replace=TRUE)

emrg_c4.expBS <- smpl.emrg$c4.1/3 + smpl.emrg$c4.2/3 + smpl.emrg$c4.3/3
hist(emrg_c4.expBS, breaks=50)
emrg_c4.expBSmn <- mean(emrg_c4.expBS)

emrg_c4.devBS <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c4"] - emrg_c4.expBSmn

t.test(EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c4"], emrg_c4.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 5
#mix_c5: LasAn, 294, 316Jeff
emrg_c5.1 <- dats23.PC$Emrg[dats23.PC$Contents=="single_LasAn"]
emrg_c5.2 <- dats23.PC$Emrg[dats23.PC$Contents=="single_294"]
emrg_c5.3 <- dats23.PC$Emrg[dats23.PC$Contents=="single_316Jeff"]

smpl.emrg$c5.1 <- sample(emrg_c5.1[!is.na(emrg_c5.1)], reps, replace=TRUE)
smpl.emrg$c5.2 <- sample(emrg_c5.2[!is.na(emrg_c5.2)], reps, replace=TRUE)
smpl.emrg$c5.3 <- sample(emrg_c5.3[!is.na(emrg_c5.3)], reps, replace=TRUE)

emrg_c5.expBS <- smpl.emrg$c5.1/3 + smpl.emrg$c5.2/3 + smpl.emrg$c5.3/3
hist(emrg_c5.expBS, breaks=50)
emrg_c5.expBSmn <- mean(emrg_c5.expBS)

emrg_c5.devBS <- EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c5"] - emrg_c5.expBSmn

t.test(EmrgByPlot.PC$Emrg[EmrgByPlot.PC$Contents=="mix_c5"], emrg_c5.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")




## Calculate 'expected' num emerged for each mix based on additive of components
#e_c1 <- ((mean(Emrg_LasAn$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
#         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
#e_c2 <- ((mean(Emrg_NM$RATE)*(400/3)) + (mean(Emrg_294$RATE)*(400/3)) 
#         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
#e_c3 <- ((mean(Emrg_UT$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
#         + (mean(Emrg_294$RATE)*(400/3)))/400
#e_c4 <- ((mean(Emrg_294$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
#         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
#e_c5 <- ((mean(Emrg_LasAn$RATE)*(400/3)) + (mean(Emrg_294$RATE)*(400/3)) 
#         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400

## Calculate deviation from expected for each plot of each close mix
#e_c1Dev <- Emrg_mixC1$RATE - e_c1
#e_c2Dev <- Emrg_mixC2$RATE - e_c2
#e_c3Dev <- Emrg_mixC3$RATE - e_c3
#e_c4Dev <- Emrg_mixC4$RATE - e_c4
#e_c5Dev <- Emrg_mixC5$RATE - e_c5

## Combine diff mixes together and add values for x-axis for plotting
#e_Devs <- c(e_c1Dev, e_c2Dev, e_c3Dev, e_c4Dev, e_c5Dev)
#xAx <- c(rep(1,length(e_c1Dev)), rep(2,length(e_c2Dev)), rep(3,length(e_c3Dev)),
#         rep(4,length(e_c4Dev)), rep(5,length(e_c5Dev)))

#e_DevsForPlot <- as.data.frame(cbind(e_Devs, xAx))

#par(pty="s")
#plot(e_DevsForPlot$xAx, e_DevsForPlot$e_Devs, pch=16, ylim=c(-0.06,0.06),
#     cex=0.9, ylab="Deviation from expected\nemergence rate", col="black",
#     xlab="Regional Mix")
#abline(h=0)


## Plots as means and error 
#EmrgAvgs.PC

#eEmrg_c1 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_LasAn"] + 
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"] +
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3
#eEmrg_c2 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_NM"] + 
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3
#eEmrg_c3 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_UT"] + 
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"])/3
#eEmrg_c4 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"] + 
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_314Jeff"])/3
#eEmrg_c5 <- (EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_LasAn"] + 
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_294"] +
#                EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="single_316Jeff"])/3

## Calculate deviation from expected for each plot of each close/ regional mix
#eEmrg_c1Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c1"] - eEmrg_c1
#eEmrg_c2Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c2"]  - eEmrg_c2
#eEmrg_c3Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c3"]  - eEmrg_c3
#eEmrg_c4Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c4"]  - eEmrg_c4
#eEmrg_c5Dev <- EmrgAvgs.PC$Emrg_MN[EmrgAvgs.PC$Contents=="mix_c5"]  - eEmrg_c5

## Combine diff mixes together and add values for x-axis for plotting
#eEmrgDevs <- c(eEmrg_c1Dev, eEmrg_c2Dev, eEmrg_c3Dev, eEmrg_c4Dev, eEmrg_c5Dev)
#xAx <- c(1:5)

#eEmrgDevsForPlot <- as.data.frame(cbind(eEmrgDevs, xAx))

#plot(eEmrgDevsForPlot$xAx, eEmrgDevsForPlot$eEmrgDevs, pch=16, ylim=c(-10,10),
#     cex=1.9, ylab="Deviation from expected\nemergence", col="black",
#     xlab="Regional Mix") 
#abline(h=0)







## 2024 --------------------------------------
#sum_data <- dats2406 %>% group_by(Site, Contents) %>% reframe(RATE=ARFR.Seedling.Count/NumSeeds)
## June

## Subset by regional mix and single type
# Site 1 (PC) only 
#mixC1 <- subset(dats2406, Contents=="mix_c1" & Site=="PC")
#mixC2 <- subset(dats2406, Contents=="mix_c2" & Site=="PC")
#mixC2 <- mixC2[order(mixC2$RATE),]
#mixC3 <- subset(dats2406, Contents=="mix_c3" & Site=="PC")
#mixC3 <- mixC3[order(mixC3$RATE),]
#mixC4 <- subset(dats2406, Contents=="mix_c4" & Site=="PC")
#mixC4 <- mixC4[order(mixC4$RATE),]
#mixC5 <- subset(dats2406, Contents=="mix_c5" & Site=="PC")
#mixC5 <- mixC5[order(mixC5$RATE),]

#singUT <- subset(dats2406, Contents=="single_UT" & Site=="PC")
#singUT <- singUT[order(singUT$RATE),]
#sing294 <- subset(dats2406, Contents=="single_294" & Site=="PC")
#sing294 <- sing294[order(sing294$RATE),]
#singNM <- subset(dats2406, Contents=="single_NM" & Site=="PC")
#singNM <- singNM[order(singNM$RATE),]
#sing316Jeff <- subset(dats2406, Contents=="single_316Jeff" & Site=="PC")
#sing316Jeff <- sing316Jeff[order(sing316Jeff$RATE),]
#sing314Jeff <- subset(dats2406, Contents=="single_314Jeff" & Site=="PC")
#sing314Jeff <- sing314Jeff[order(sing314Jeff$RATE),]
#singLasAn <- subset(dats2406, Contents=="single_LasAn" & Site=="PC")
#singLasAn <- singLasAn[order(singLasAn$RATE),]

## List what components are in each mix
#mix_c1: LasAn, 314Jeff, 316Jeff
#mix_c2: NM, 294, 316Jeff
#mix_c3: UT, 294, 314Jeff
#mix_c4: 294, 314Jeff, 316Jeff
#mix_c5: LasAn, 294, 316Jeff


## Counts 
#eCount_c1 <- (mean(singLasAn$ARFR.Seedling.Count) + mean(sing314Jeff$ARFR.Seedling.Count) 
#         + mean(sing316Jeff$ARFR.Seedling.Count))/3
#eCount_c2 <- (mean(singNM$ARFR.Seedling.Count) + mean(sing294$ARFR.Seedling.Count) 
#         + mean(sing316Jeff$ARFR.Seedling.Count))/3
#eCount_c3 <- (mean(singUT$ARFR.Seedling.Count) + mean(sing314Jeff$ARFR.Seedling.Count) 
#         + mean(sing294$ARFR.Seedling.Count))/3
#eCount_c4 <- (mean(sing294$ARFR.Seedling.Count) + mean(sing314Jeff$ARFR.Seedling.Count) 
#         + mean(sing316Jeff$ARFR.Seedling.Count))/3
#eCount_c5 <- (mean(singLasAn$ARFR.Seedling.Count) + mean(sing294$ARFR.Seedling.Count) 
#         + mean(sing316Jeff$ARFR.Seedling.Count))/3

CountAvgs.PC

## Calculate 'expected' count for each mix based on additive of components
eCount_c1 <- (CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_LasAn"] + 
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_314Jeff"] +
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_316Jeff"])/3
eCount_c2 <- (CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_NM"] + 
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_294"] +
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_316Jeff"])/3
eCount_c3 <- (CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_UT"] + 
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_294"] +
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_314Jeff"])/3
eCount_c4 <- (CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_316Jeff"] + 
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_294"] +
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_314Jeff"])/3
eCount_c5 <- (CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_LasAn"] + 
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_294"] +
               CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="single_316Jeff"])/3

CountByPlot.PC <- dats2406.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Count=mean(ARFR.Seedling.Count,na.rm=TRUE))

## Calculate deviation from expected for each plot of each close mix
eCount_c1Dev <- CountByPlot.PC$Count[CountByPlot.PC$Contents=="mix_c1"] - eCount_c1
eCount_c2Dev <- CountByPlot.PC$Count[CountByPlot.PC$Contents=="mix_c2"] - eCount_c2
eCount_c3Dev <- CountByPlot.PC$Count[CountByPlot.PC$Contents=="mix_c3"] - eCount_c3
eCount_c4Dev <- CountByPlot.PC$Count[CountByPlot.PC$Contents=="mix_c4"] - eCount_c4
eCount_c5Dev <- CountByPlot.PC$Count[CountByPlot.PC$Contents=="mix_c5"] - eCount_c5

## Combine diff mixes together and add values for x-axis for plotting
eCountDevs <- c(eCount_c1Dev, eCount_c2Dev, eCount_c3Dev, eCount_c4Dev, eCount_c5Dev)
xAx <- c(rep(1,length(eCount_c1Dev)), rep(2,length(eCount_c2Dev)), rep(3,length(eCount_c3Dev)),
         rep(4,length(eCount_c4Dev)), rep(5,length(eCount_c5Dev)))

eCountDevsForPlot <- as.data.frame(cbind(eCountDevs, xAx))

par(pty="s")
plot(eCountDevsForPlot$xAx, eCountDevsForPlot$eCountDevs, pch=16, ylim=c(-25,25),
     cex=0.9, ylab="Deviation from expected\ncounts", col="grey45",
     xlab="Regional Mix")
abline(h=0)



## Plot means 
## Calculate deviation from expected for each plot of each close/ regional mix
eMnCount_c1Dev <- CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="mix_c1"] - eCount_c1
eMnCount_c2Dev <- CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="mix_c2"]  - eCount_c2
eMnCount_c3Dev <- CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="mix_c3"]  - eCount_c3
eMnCount_c4Dev <- CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="mix_c4"]  - eCount_c4
eMnCount_c5Dev <- CountAvgs.PC$Count_MN[CountAvgs.PC$Contents=="mix_c5"]  - eCount_c5

## Combine diff mixes together and add values for x-axis for plotting
eMnCountDevs <- c(eMnCount_c1Dev, eMnCount_c2Dev, eMnCount_c3Dev, eMnCount_c4Dev, eMnCount_c5Dev)
xAx <- c(1:5)

eMnCountDevsForPlot <- as.data.frame(cbind(eMnCountDevs, xAx))

points(eMnCountDevsForPlot$xAx, eMnCountDevsForPlot$eMnCountDevs, pch=16,
     cex=2, col="black") 






## Survival
## Calculate 'expected' survival for each mix based on additive of components
#eSurv_c1 <- (mean(singLasAn$Surv) + mean(sing314Jeff$Surv) 
#              + mean(sing316Jeff$Surv))/3
#eSurv_c2 <- (mean(singNM$Surv) + mean(sing294$Surv) 
#              + mean(sing316Jeff$Surv))/3
#eSurv_c3 <- (mean(singUT$Surv, na.rm=TRUE) + mean(sing314Jeff$Surv) 
#              + mean(sing294$Surv))/3
#eSurv_c4 <- (mean(sing294$Surv) + mean(sing314Jeff$Surv) 
#              + mean(sing316Jeff$Surv))/3
#eSurv_c5 <- (mean(singLasAn$Surv) + mean(sing294$Surv) 
#              + mean(sing316Jeff$Surv))/3

SurvAvgs.PC

eSurv_c1 <- (SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_LasAn"] + 
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_314Jeff"] +
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_316Jeff"])/3
eSurv_c2 <- (SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_NM"] + 
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_294"] +
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_316Jeff"])/3
eSurv_c3 <- (SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_UT"] + 
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_294"] +
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_314Jeff"])/3
eSurv_c4 <- (SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_316Jeff"] + 
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_294"] +
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_314Jeff"])/3
eSurv_c5 <- (SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_LasAn"] + 
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_294"] +
               SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="single_316Jeff"])/3


SurvByPlot.PC <- dats2406.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Surv=mean(Surv,na.rm=TRUE))

## Calculate deviation from expected for each plot of each close mix
eSurv_c1Dev <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c1"]- eSurv_c1
eSurv_c2Dev <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c2"] - eSurv_c2
eSurv_c3Dev <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c3"] - eSurv_c3
eSurv_c4Dev <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c4"] - eSurv_c4
eSurv_c5Dev <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c5"] - eSurv_c5


## Combine diff mixes together and add values for x-axis for plotting
eSurvDevs <- c(eSurv_c1Dev, eSurv_c2Dev, eSurv_c3Dev, eSurv_c4Dev, eSurv_c5Dev)
xAx <- c(rep(1,length(eSurv_c1Dev)), rep(2,length(eSurv_c2Dev)), rep(3,length(eSurv_c3Dev)),
         rep(4,length(eSurv_c4Dev)), rep(5,length(eSurv_c5Dev)))

eSurvDevsForPlot <- as.data.frame(cbind(eSurvDevs, xAx))

plot(eSurvDevsForPlot$xAx, eSurvDevsForPlot$eSurvDevs, pch=16, ylim=c(-0.5,0.5),
     cex=0.9, ylab="Deviation from expected\nsurvival rate", col="grey45",
     xlab="Regional Mix")
abline(h=0)



## Plots as means and error 
## Calculate deviation from expected for each plot of each close/ regional mix
eMnSurv_c1Dev <- SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="mix_c1"] - eSurv_c1
eMnSurv_c2Dev <- SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="mix_c2"]  - eSurv_c2
eMnSurv_c3Dev <- SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="mix_c3"]  - eSurv_c3
eMnSurv_c4Dev <- SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="mix_c4"]  - eSurv_c4
eMnSurv_c5Dev <- SurvAvgs.PC$Surv_MN[SurvAvgs.PC$Contents=="mix_c5"]  - eSurv_c5


## Combine diff mixes together and add values for x-axis for plotting
eMnSurvDevs <- c(eMnSurv_c1Dev, eMnSurv_c2Dev, eMnSurv_c3Dev, eMnSurv_c4Dev, eMnSurv_c5Dev)
xAx <- c(1:5)

eMnSurvDevsForPlot <- as.data.frame(cbind(eMnSurvDevs, xAx))

points(eMnSurvDevsForPlot$xAx, eMnSurvDevsForPlot$eMnSurvDevs, pch=16,
     cex=2, col="black") 
## --------------------




## BOOTSTRAP TO GENERATE EXPECTED MIX PERFORMANCE
reps <- 100
smpl.surv <- as.data.frame(matrix(NA, reps, (5*3)))
colnames(smpl.surv) <- c("c1.1","c1.2","c1.3","c2.1","c2.2","c2.3","c3.1","c3.2","c3.3",
                         "c4.1","c4.2","c4.3","c5.1","c5.2","c5.3")


unique(dats$Contents)

## Regional mix 1
#mix_c1: LasAn, 314Jeff, 316Jeff
surv_c1.1 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_LasAn"]
surv_c1.2 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_314Jeff"]
surv_c1.3 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_316Jeff"]

smpl.surv$c1.1 <- sample(surv_c1.1[!is.na(surv_c1.1)], reps, replace=TRUE)
smpl.surv$c1.2 <- sample(surv_c1.2[!is.na(surv_c1.2)], reps, replace=TRUE)
smpl.surv$c1.3 <- sample(surv_c1.3[!is.na(surv_c1.3)], reps, replace=TRUE)

surv_c1.expBS <- smpl.surv$c1.1/3 + smpl.surv$c1.2/3 + smpl.surv$c1.3/3
hist(surv_c1.expBS, breaks=100)
surv_c1.expBSmn <- mean(surv_c1.expBS)

surv_c1.devBS <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c1"] - surv_c1.expBSmn

t.test(SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c1"], surv_c1.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 2
#mix_c2: NM, 294, 316Jeff
surv_c2.1 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_NM"]
surv_c2.2 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_294"]
surv_c2.3 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_316Jeff"]

smpl.surv$c2.1 <- sample(surv_c2.1[!is.na(surv_c2.1)], reps, replace=TRUE)
smpl.surv$c2.2 <- sample(surv_c2.2[!is.na(surv_c2.2)], reps, replace=TRUE)
smpl.surv$c2.3 <- sample(surv_c2.3[!is.na(surv_c2.3)], reps, replace=TRUE)

surv_c2.expBS <- smpl.surv$c2.1/3 + smpl.surv$c2.2/3 + smpl.surv$c2.3/3
hist(surv_c2.expBS, breaks=50)
surv_c2.expBSmn <- mean(surv_c2.expBS)

surv_c2.devBS <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c2"] - surv_c2.expBSmn

t.test(SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c2"], surv_c2.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 3
#mix_c3: UT, 294, 314Jeff
surv_c3.1 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_UT"]
surv_c3.2 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_294"]
surv_c3.3 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_314Jeff"]

smpl.surv$c3.1 <- sample(surv_c3.1[!is.na(surv_c3.1)], reps, replace=TRUE)
smpl.surv$c3.2 <- sample(surv_c3.2[!is.na(surv_c3.2)], reps, replace=TRUE)
smpl.surv$c3.3 <- sample(surv_c3.3[!is.na(surv_c3.3)], reps, replace=TRUE)

surv_c3.expBS <- smpl.surv$c3.1/3 + smpl.surv$c3.2/3 + smpl.surv$c3.3/3
hist(surv_c3.expBS, breaks=50)
surv_c3.expBSmn <- mean(surv_c3.expBS)

surv_c3.devBS <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c3"] - surv_c3.expBSmn

t.test(SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c3"], surv_c3.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 4
#mix_c4: 294, 314Jeff, 316Jeff
surv_c4.1 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_294"]
surv_c4.2 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_314Jeff"]
surv_c4.3 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_316Jeff"]

smpl.surv$c4.1 <- sample(surv_c4.1[!is.na(surv_c4.1)], reps, replace=TRUE)
smpl.surv$c4.2 <- sample(surv_c4.2[!is.na(surv_c4.2)], reps, replace=TRUE)
smpl.surv$c4.3 <- sample(surv_c4.3[!is.na(surv_c4.3)], reps, replace=TRUE)

surv_c4.expBS <- smpl.surv$c4.1/3 + smpl.surv$c4.2/3 + smpl.surv$c4.3/3
hist(surv_c4.expBS, breaks=50)
surv_c4.expBSmn <- mean(surv_c4.expBS)

surv_c4.devBS <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c4"] - surv_c4.expBSmn

t.test(SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c4"], surv_c4.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 5
#mix_c5: LasAn, 294, 316Jeff
surv_c5.1 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_LasAn"]
surv_c5.2 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_294"]
surv_c5.3 <- dats2406.PC$Surv[dats2406.PC$Contents=="single_316Jeff"]

smpl.surv$c5.1 <- sample(surv_c5.1[!is.na(surv_c5.1)], reps, replace=TRUE)
smpl.surv$c5.2 <- sample(surv_c5.2[!is.na(surv_c5.2)], reps, replace=TRUE)
smpl.surv$c5.3 <- sample(surv_c5.3[!is.na(surv_c5.3)], reps, replace=TRUE)

surv_c5.expBS <- smpl.surv$c5.1/3 + smpl.surv$c5.2/3 + smpl.surv$c5.3/3
hist(surv_c5.expBS, breaks=50)
surv_c5.expBSmn <- mean(surv_c5.expBS)

surv_c5.devBS <- SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c5"] - surv_c5.expBSmn

t.test(SurvByPlot.PC$Surv[SurvByPlot.PC$Contents=="mix_c5"], surv_c5.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")





## Percent Veg Cover
## Calculate 'expected' cover for each mix based on additive of components
#eCovr_c1 <- (mean(singLasAn$ARFR.Percent.Cover) + mean(sing314Jeff$ARFR.Percent.Cover) 
#             + mean(sing316Jeff$ARFR.Percent.Cover,na.rm=TRUE))/3
#eCovr_c2 <- (mean(singNM$ARFR.Percent.Cover) + mean(sing294$ARFR.Percent.Cover) 
#             + mean(sing316Jeff$ARFR.Percent.Cover,na.rm=TRUE))/3
#eCovr_c3 <- (mean(singUT$ARFR.Percent.Cover, na.rm=TRUE) + mean(sing314Jeff$ARFR.Percent.Cover) 
#             + mean(sing294$ARFR.Percent.Cover))/3
#eCovr_c4 <- (mean(sing294$ARFR.Percent.Cover) + mean(sing314Jeff$ARFR.Percent.Cover) 
#             + mean(sing316Jeff$ARFR.Percent.Cover,na.rm=TRUE))/3
#eCovr_c5 <- (mean(singLasAn$ARFR.Percent.Cover) + mean(sing294$ARFR.Percent.Cover) 
#             + mean(sing316Jeff$ARFR.Percent.Cover,na.rm=TRUE))/3

CovrAvgs.PC

eCovr_c1 <- (CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_LasAn"] + 
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_314Jeff"] +
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_316Jeff"])/3
eCovr_c2 <- (CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_NM"] + 
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_294"] +
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_316Jeff"])/3
eCovr_c3 <- (CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_UT"] + 
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_294"] +
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_314Jeff"])/3
eCovr_c4 <- (CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_316Jeff"] + 
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_294"] +
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_314Jeff"])/3
eCovr_c5 <- (CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_LasAn"] + 
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_294"] +
               CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="single_316Jeff"])/3

CovrByPlot.PC <- dats2406.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Covr=mean(ARFR.Percent.Cover,na.rm=TRUE))

## Calculate deviation from expected cover for each plot of each regional mix
eCovr_c1Dev <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c1"] - eCovr_c1
eCovr_c2Dev <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c2"] - eCovr_c2
eCovr_c3Dev <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c3"] - eCovr_c3
eCovr_c4Dev <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c4"] - eCovr_c4
eCovr_c5Dev <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c5"] - eCovr_c5


## Combine diff mixes together & add values for x-axis for plotting
eCovrDevs <- c(eCovr_c1Dev, eCovr_c2Dev, eCovr_c3Dev, eCovr_c4Dev, eCovr_c5Dev)
xAx <- c(rep(1,length(eCovr_c1Dev)), rep(2,length(eCovr_c2Dev)), rep(3,length(eCovr_c3Dev)),
         rep(4,length(eCovr_c4Dev)), rep(5,length(eCovr_c5Dev)))

eCovrDevsForPlot <- as.data.frame(cbind(eCovrDevs, xAx))

plot(eCovrDevsForPlot$xAx, eCovrDevsForPlot$eCovrDevs, pch=16, ylim=c(-45,45),
     cex=0.9, ylab="Deviation from expected\nvegetative cover (% cover)", col="grey45",
     xlab="Regional Mix") 
abline(h=0)


## Plot means 
## Calculate deviation from expected mean (avg of all plots) of each reg mix
eMnCovr_c1Dev <- CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="mix_c1"] - eCovr_c1
eMnCovr_c2Dev <- CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="mix_c2"]  - eCovr_c2
eMnCovr_c3Dev <- CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="mix_c3"]  - eCovr_c3
eMnCovr_c4Dev <- CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="mix_c4"]  - eCovr_c4
eMnCovr_c5Dev <- CovrAvgs.PC$Covr_MN[CovrAvgs.PC$Contents=="mix_c5"]  - eCovr_c5

## Combine diff mixes together and add values for x-axis for plotting
eMnCovrDevs <- c(eMnCovr_c1Dev, eMnCovr_c2Dev, eMnCovr_c3Dev, eMnCovr_c4Dev, eMnCovr_c5Dev)
xAx <- c(1:5)

eMnCovrDevsForPlot <- as.data.frame(cbind(eMnCovrDevs, xAx))

points(eMnCovrDevsForPlot$xAx, eMnCovrDevsForPlot$eMnCovrDevs, pch=16,
     cex=2, col="black") 





## Site 2 (Eco Park)
CovrAvgs.EP <- dats2406.EP %>% group_by(Contents) %>% 
  dplyr::summarise(Covr_MD=median(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_MN=mean(ARFR.Percent.Cover,na.rm=TRUE),
                   Covr_SE=calcSE(ARFR.Percent.Cover))

eCovr_c1.EP <- (CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_LasAn"] + 
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_314Jeff"] +
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_316Jeff"])/3
eCovr_c2.EP <- (CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_NM"] + 
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_294"] +
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_316Jeff"])/3
eCovr_c3.EP <- (CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_UT"] + 
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_294"] +
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_314Jeff"])/3
eCovr_c4.EP <- (CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_316Jeff"] + 
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_294"] +
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_314Jeff"])/3
eCovr_c5.EP <- (CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_LasAn"] + 
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_294"] +
               CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="single_316Jeff"])/3

CovrByPlot.EP <- dats2406.EP %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Covr=mean(ARFR.Percent.Cover,na.rm=TRUE))

## Calculate deviation from expected cover for each plot of each regional mix
eCovr_c1Dev.EP <- CovrByPlot.EP$Covr[CovrByPlot.EP$Contents=="mix_c1"] - eCovr_c1.EP
eCovr_c2Dev.EP <- CovrByPlot.EP$Covr[CovrByPlot.EP$Contents=="mix_c2"] - eCovr_c2.EP
eCovr_c3Dev.EP <- CovrByPlot.EP$Covr[CovrByPlot.EP$Contents=="mix_c3"] - eCovr_c3.EP
eCovr_c4Dev.EP <- CovrByPlot.EP$Covr[CovrByPlot.EP$Contents=="mix_c4"] - eCovr_c4.EP
eCovr_c5Dev.EP <- CovrByPlot.EP$Covr[CovrByPlot.EP$Contents=="mix_c5"] - eCovr_c5.EP


## Combine diff mixes together & add values for x-axis for plotting
eCovrDevs.EP <- c(eCovr_c1Dev.EP, eCovr_c2Dev.EP, eCovr_c3Dev.EP, eCovr_c4Dev.EP, eCovr_c5Dev.EP)
xAx <- c(rep(1,length(eCovr_c1Dev.EP)), rep(2,length(eCovr_c2Dev.EP)), rep(3,length(eCovr_c3Dev.EP)),
         rep(4,length(eCovr_c4Dev.EP)), rep(5,length(eCovr_c5Dev.EP)))

eCovrDevs.ForPlotEP <- as.data.frame(cbind(eCovrDevs.EP, xAx))

plot(eCovrDevs.ForPlotEP$xAx, eCovrDevs.ForPlotEP$eCovrDevs.EP, pch=16, ylim=c(-35,35),
     cex=0.9, ylab="Deviation from expected\nvegetative cover (% cover)", col="grey45",
     xlab="Regional Mix", main="Site 2") 
abline(h=0)


## Plots as means and error 
## Calculate deviation from expected mean (avg of all plots) of each reg mix
eMnCovr_c1Dev.EP <- CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="mix_c1"] - eCovr_c1.EP
eMnCovr_c2Dev.EP <- CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="mix_c2"]  - eCovr_c2.EP
eMnCovr_c3Dev.EP <- CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="mix_c3"]  - eCovr_c3.EP
eMnCovr_c4Dev.EP <- CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="mix_c4"]  - eCovr_c4.EP
eMnCovr_c5Dev.EP <- CovrAvgs.EP$Covr_MN[CovrAvgs.EP$Contents=="mix_c5"]  - eCovr_c5.EP

## Combine diff mixes together and add values for x-axis for plotting
eMnCovrDevs.EP <- c(eMnCovr_c1Dev.EP, eMnCovr_c2Dev.EP, eMnCovr_c3Dev.EP, eMnCovr_c4Dev.EP, eMnCovr_c5Dev.EP)
xAx <- c(1:5)

eMnCovrDevsForPlotEP <- as.data.frame(cbind(eMnCovrDevs.EP, xAx))

points(eMnCovrDevsForPlotEP$xAx, eMnCovrDevsForPlotEP$eMnCovrDevs.EP, pch=16,
       cex=2, col="black") 
## -------------------------



## BOOTSTRAP TO GENERATE EXPECTED MIX PERFORMANCE
reps <- 100
smpl.covr <- as.data.frame(matrix(NA, reps, (5*3)))
colnames(smpl.covr) <- c("c1.1","c1.2","c1.3","c2.1","c2.2","c2.3","c3.1","c3.2","c3.3",
                          "c4.1","c4.2","c4.3","c5.1","c5.2","c5.3")


unique(dats$Contents)

## Regional mix 1
#mix_c1: LasAn, 314Jeff, 316Jeff
covr_c1.1 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_LasAn"]
covr_c1.2 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_314Jeff"]
covr_c1.3 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_316Jeff"]

smpl.covr$c1.1 <- sample(covr_c1.1[!is.na(covr_c1.1)], reps, replace=TRUE)
smpl.covr$c1.2 <- sample(covr_c1.2[!is.na(covr_c1.2)], reps, replace=TRUE)
smpl.covr$c1.3 <- sample(covr_c1.3[!is.na(covr_c1.3)], reps, replace=TRUE)

covr_c1.expBS <- smpl.covr$c1.1/3 + smpl.covr$c1.2/3 + smpl.covr$c1.3/3
hist(covr_c1.expBS, breaks=100)
covr_c1.expBSmn <- mean(covr_c1.expBS)

covr_c1.devBS <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c1"] - covr_c1.expBSmn

t.test(CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c1"], covr_c1.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 2
#mix_c2: NM, 294, 316Jeff
covr_c2.1 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_NM"]
covr_c2.2 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_294"]
covr_c2.3 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_316Jeff"]

smpl.covr$c2.1 <- sample(covr_c2.1[!is.na(covr_c2.1)], reps, replace=TRUE)
smpl.covr$c2.2 <- sample(covr_c2.2[!is.na(covr_c2.2)], reps, replace=TRUE)
smpl.covr$c2.3 <- sample(covr_c2.3[!is.na(covr_c2.3)], reps, replace=TRUE)

covr_c2.expBS <- smpl.covr$c2.1/3 + smpl.covr$c2.2/3 + smpl.covr$c2.3/3
hist(covr_c2.expBS, breaks=50)
covr_c2.expBSmn <- mean(covr_c2.expBS)

covr_c2.devBS <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c2"] - covr_c2.expBSmn

t.test(CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c2"], covr_c2.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 3
#mix_c3: UT, 294, 314Jeff
covr_c3.1 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_UT"]
covr_c3.2 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_294"]
covr_c3.3 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_314Jeff"]

smpl.covr$c3.1 <- sample(covr_c3.1[!is.na(covr_c3.1)], reps, replace=TRUE)
smpl.covr$c3.2 <- sample(covr_c3.2[!is.na(covr_c3.2)], reps, replace=TRUE)
smpl.covr$c3.3 <- sample(covr_c3.3[!is.na(covr_c3.3)], reps, replace=TRUE)

covr_c3.expBS <- smpl.covr$c3.1/3 + smpl.covr$c3.2/3 + smpl.covr$c3.3/3
hist(covr_c3.expBS, breaks=50)
covr_c3.expBSmn <- mean(covr_c3.expBS)

covr_c3.devBS <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c3"] - covr_c3.expBSmn

t.test(CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c3"], covr_c3.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 4
#mix_c4: 294, 314Jeff, 316Jeff
covr_c4.1 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_294"]
covr_c4.2 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_314Jeff"]
covr_c4.3 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_316Jeff"]

smpl.covr$c4.1 <- sample(covr_c4.1[!is.na(covr_c4.1)], reps, replace=TRUE)
smpl.covr$c4.2 <- sample(covr_c4.2[!is.na(covr_c4.2)], reps, replace=TRUE)
smpl.covr$c4.3 <- sample(covr_c4.3[!is.na(covr_c4.3)], reps, replace=TRUE)

covr_c4.expBS <- smpl.covr$c4.1/3 + smpl.covr$c4.2/3 + smpl.covr$c4.3/3
hist(covr_c4.expBS, breaks=50)
covr_c4.expBSmn <- mean(covr_c4.expBS)

covr_c4.devBS <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c4"] - covr_c4.expBSmn

t.test(CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c4"], covr_c4.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 5
#mix_c5: LasAn, 294, 316Jeff
covr_c5.1 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_LasAn"]
covr_c5.2 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_294"]
covr_c5.3 <- dats2406.PC$ARFR.Percent.Cover[dats2406.PC$Contents=="single_316Jeff"]

smpl.covr$c5.1 <- sample(covr_c5.1[!is.na(covr_c5.1)], reps, replace=TRUE)
smpl.covr$c5.2 <- sample(covr_c5.2[!is.na(covr_c5.2)], reps, replace=TRUE)
smpl.covr$c5.3 <- sample(covr_c5.3[!is.na(covr_c5.3)], reps, replace=TRUE)

covr_c5.expBS <- smpl.covr$c5.1/3 + smpl.covr$c5.2/3 + smpl.covr$c5.3/3
hist(covr_c5.expBS, breaks=50)
covr_c5.expBSmn <- mean(covr_c5.expBS)

covr_c5.devBS <- CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c5"] - covr_c5.expBSmn

t.test(CovrByPlot.PC$Covr[CovrByPlot.PC$Contents=="mix_c5"], covr_c5.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")
## -------------------------------------------------------------------------------






## Repro cover
ReproAvgs.PC

## Calculate 'expected' repro for each regional mix based on additive of components
eRepro_c1 <- (ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_LasAn"] + 
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_314Jeff"] +
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_316Jeff"])/3
eRepro_c2 <- (ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_NM"] + 
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_294"] +
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_316Jeff"])/3
eRepro_c3 <- (ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_UT"] + 
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_294"] +
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_314Jeff"])/3
eRepro_c4 <- (ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_316Jeff"] + 
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_294"] +
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_314Jeff"])/3
eRepro_c5 <- (ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_LasAn"] + 
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_294"] +
                ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="single_316Jeff"])/3

ReproByPlot.PC <- dats2409.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Repro=mean(ARFR.Reproductive.Percent.Cover,na.rm=TRUE))

## Calculate deviation from expected for each plot of each close mix
eRepro_c1Dev <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c1"] - eRepro_c1
eRepro_c2Dev <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c2"] - eRepro_c2
eRepro_c3Dev <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c3"] - eRepro_c3
eRepro_c4Dev <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c4"] - eRepro_c4
eRepro_c5Dev <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c5"] - eRepro_c5

## Combine diff mixes together and add values for x-axis for plotting
eReproDevs <- c(eRepro_c1Dev, eRepro_c2Dev, eRepro_c3Dev, eRepro_c4Dev, eRepro_c5Dev)
xAx <- c(rep(1,length(eRepro_c1Dev)), rep(2,length(eRepro_c2Dev)), rep(3,length(eRepro_c3Dev)),
         rep(4,length(eRepro_c4Dev)), rep(5,length(eRepro_c5Dev)))

eReproDevsForPlot <- as.data.frame(cbind(eReproDevs, xAx))

plot(eReproDevsForPlot$xAx, eReproDevsForPlot$eReproDevs, pch=16, ylim=c(-60,60),
     cex=0.9, ylab="Deviation from expected\nreproductive cover (% cover)", col="grey45",
     xlab="Regional Mix")
abline(h=0)



## Add plot means
## Calculate deviation from expected mean (avg of all plots) of each reg mix
eMnRepro_c1Dev <- ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="mix_c1"] - eRepro_c1
eMnRepro_c2Dev <- ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="mix_c2"]  - eRepro_c2
eMnRepro_c3Dev <- ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="mix_c3"]  - eRepro_c3
eMnRepro_c4Dev <- ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="mix_c4"]  - eRepro_c4
eMnRepro_c5Dev <- ReproAvgs.PC$Repro_MN[ReproAvgs.PC$Contents=="mix_c5"]  - eRepro_c5


## Combine diff mixes together and add values for x-axis for plotting
eMnReproDevs <- c(eMnRepro_c1Dev, eMnRepro_c2Dev, eMnRepro_c3Dev, eMnRepro_c4Dev, eMnRepro_c5Dev)
xAx <- c(1:5)

eMnReproDevsForPlot <- as.data.frame(cbind(eMnReproDevs, xAx))

points(eMnReproDevsForPlot$xAx, eMnReproDevsForPlot$eMnReproDevs, pch=16,
     cex=2, col="black") 



## Site 2 (Eco Park) (too many missing values?)
#ReproAvgs.EP <- dats2409.EP %>% group_by(Contents) %>% 
#  dplyr::summarise(Repro_MD=median(ARFR.Percent.Cover,na.rm=TRUE),
#                   Repro_MN=mean(ARFR.Percent.Cover,na.rm=TRUE),
#                   Repro_SE=calcSE(ARFR.Percent.Cover))

#eRepro_c1.EP <- (ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_LasAn"] + 
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_314Jeff"] +
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_316Jeff"])/3
#eRepro_c2.EP <- (ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_NM"] + 
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_294"] +
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_316Jeff"])/3
#eRepro_c3.EP <- (ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_UT"] + 
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_294"] +
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_314Jeff"])/3
#eRepro_c4.EP <- (ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_316Jeff"] + 
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_294"] +
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_314Jeff"])/3
#eRepro_c5.EP <- (ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_LasAn"] + 
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_294"] +
#                  ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="single_316Jeff"])/3

#ReproByPlot.EP <- dats2409.EP %>% group_by(Contents, Plot.Number) %>% 
#  dplyr::summarise(Repro=mean(ARFR.Reproductive.Percent.Cover,na.rm=TRUE))

## Calculate deviation from expected cover for each plot of each regional mix
#eRepro_c1Dev.EP <- ReproByPlot.EP$Repro[ReproByPlot.EP$Contents=="mix_c1"] - eRepro_c1.EP
#eRepro_c2Dev.EP <- ReproByPlot.EP$Repro[ReproByPlot.EP$Contents=="mix_c2"] - eRepro_c2.EP
#eRepro_c3Dev.EP <- ReproByPlot.EP$Repro[ReproByPlot.EP$Contents=="mix_c3"] - eRepro_c3.EP
#eRepro_c4Dev.EP <- ReproByPlot.EP$Repro[ReproByPlot.EP$Contents=="mix_c4"] - eRepro_c4.EP
#eRepro_c5Dev.EP <- ReproByPlot.EP$Repro[ReproByPlot.EP$Contents=="mix_c5"] - eRepro_c5.EP


## Combine diff mixes together & add values for x-axis for plotting
#eReproDevs.EP <- c(eRepro_c1Dev.EP, eRepro_c2Dev.EP, eRepro_c3Dev.EP, eRepro_c4Dev.EP, eRepro_c5Dev.EP)
#xAx <- c(rep(1,length(eRepro_c1Dev.EP)), rep(2,length(eRepro_c2Dev.EP)), rep(3,length(eRepro_c3Dev.EP)),
#         rep(4,length(eRepro_c4Dev.EP)), rep(5,length(eRepro_c5Dev.EP)))

#eReproDevs.ForPlotEP <- as.data.frame(cbind(eReproDevs.EP, xAx))

#plot(eReproDevs.ForPlotEP$xAx, eReproDevs.ForPlotEP$eReproDevs.EP, pch=16, ylim=c(-35,35),
#     cex=0.9, ylab="Deviation from expected\nreproductive cover (% cover)", col="grey45",
#     xlab="Regional Mix", main="Site 2") 
#abline(h=0)


## Plots as means and error 
## Calculate deviation from expected mean (avg of all plots) of each reg mix
#eMnRepro_c1Dev.EP <- ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="mix_c1"] - eRepro_c1.EP
#eMnRepro_c2Dev.EP <- ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="mix_c2"]  - eRepro_c2.EP
#eMnRepro_c3Dev.EP <- ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="mix_c3"]  - eRepro_c3.EP
#eMnRepro_c4Dev.EP <- ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="mix_c4"]  - eRepro_c4.EP
#eMnRepro_c5Dev.EP <- ReproAvgs.EP$Repro_MN[ReproAvgs.EP$Contents=="mix_c5"]  - eRepro_c5.EP

## Combine diff mixes together and add values for x-axis for plotting
#eMnReproDevs.EP <- c(eMnRepro_c1Dev.EP, eMnRepro_c2Dev.EP, eMnRepro_c3Dev.EP, eMnRepro_c4Dev.EP, eMnRepro_c5Dev.EP)
#xAx <- c(1:5)

#eMnReproDevsForPlotEP <- as.data.frame(cbind(eMnReproDevs.EP, xAx))

#points(eMnReproDevsForPlotEP$xAx, eMnReproDevsForPlotEP$eMnReproDevs.EP, pch=16,
#       cex=2, col="black") 




## BOOTSTRAP TO GENERATE EXPECTED MIX PERFORMANCE
reps <- 100
smpl.repro <- as.data.frame(matrix(NA, reps, (5*3)))
colnames(smpl.repro) <- c("c1.1","c1.2","c1.3","c2.1","c2.2","c2.3","c3.1","c3.2","c3.3",
                      "c4.1","c4.2","c4.3","c5.1","c5.2","c5.3")


unique(dats$Contents)

## Regional mix 1
#mix_c1: LasAn, 314Jeff, 316Jeff
repro_c1.1 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_LasAn"]
repro_c1.2 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_314Jeff"]
repro_c1.3 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_316Jeff"]

smpl.repro$c1.1 <- sample(repro_c1.1[!is.na(repro_c1.1)], reps, replace=TRUE)
smpl.repro$c1.2 <- sample(repro_c1.2[!is.na(repro_c1.2)], reps, replace=TRUE)
smpl.repro$c1.3 <- sample(repro_c1.3[!is.na(repro_c1.3)], reps, replace=TRUE)

repro_c1.expBS <- smpl.repro$c1.1/3 + smpl.repro$c1.2/3 + smpl.repro$c1.3/3
hist(repro_c1.expBS, breaks=100)
repro_c1.expBSmn <- mean(repro_c1.expBS)

repro_c1.devBS <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c1"] - repro_c1.expBSmn

t.test(ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c1"], repro_c1.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 2
#mix_c2: NM, 294, 316Jeff
repro_c2.1 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_NM"]
repro_c2.2 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_294"]
repro_c2.3 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_316Jeff"]

smpl.repro$c2.1 <- sample(repro_c2.1[!is.na(repro_c2.1)], reps, replace=TRUE)
smpl.repro$c2.2 <- sample(repro_c2.2[!is.na(repro_c2.2)], reps, replace=TRUE)
smpl.repro$c2.3 <- sample(repro_c2.3[!is.na(repro_c2.3)], reps, replace=TRUE)

repro_c2.expBS <- smpl.repro$c2.1/3 + smpl.repro$c2.2/3 + smpl.repro$c2.3/3
hist(repro_c2.expBS, breaks=50)
repro_c2.expBSmn <- mean(repro_c2.expBS)

repro_c2.devBS <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c2"] - repro_c2.expBSmn

t.test(ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c2"], repro_c2.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 3
#mix_c3: UT, 294, 314Jeff
repro_c3.1 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_UT"]
repro_c3.2 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_294"]
repro_c3.3 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_314Jeff"]

smpl.repro$c3.1 <- sample(repro_c3.1[!is.na(repro_c3.1)], reps, replace=TRUE)
smpl.repro$c3.2 <- sample(repro_c3.2[!is.na(repro_c3.2)], reps, replace=TRUE)
smpl.repro$c3.3 <- sample(repro_c3.3[!is.na(repro_c3.3)], reps, replace=TRUE)

repro_c3.expBS <- smpl.repro$c3.1/3 + smpl.repro$c3.2/3 + smpl.repro$c3.3/3
hist(repro_c3.expBS, breaks=50)
repro_c3.expBSmn <- mean(repro_c3.expBS)

repro_c3.devBS <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c3"] - repro_c3.expBSmn

t.test(ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c3"], repro_c3.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 4
#mix_c4: 294, 314Jeff, 316Jeff
repro_c4.1 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_294"]
repro_c4.2 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_314Jeff"]
repro_c4.3 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_316Jeff"]

smpl.repro$c4.1 <- sample(repro_c4.1[!is.na(repro_c4.1)], reps, replace=TRUE)
smpl.repro$c4.2 <- sample(repro_c4.2[!is.na(repro_c4.2)], reps, replace=TRUE)
smpl.repro$c4.3 <- sample(repro_c4.3[!is.na(repro_c4.3)], reps, replace=TRUE)

repro_c4.expBS <- smpl.repro$c4.1/3 + smpl.repro$c4.2/3 + smpl.repro$c4.3/3
hist(repro_c4.expBS, breaks=50)
repro_c4.expBSmn <- mean(repro_c4.expBS)

repro_c4.devBS <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c4"] - repro_c4.expBSmn

t.test(ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c4"], repro_c4.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")


## Regional mix 5
#mix_c5: LasAn, 294, 316Jeff
repro_c5.1 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_LasAn"]
repro_c5.2 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_294"]
repro_c5.3 <- dats2409.PC$ARFR.Reproductive.Percent.Cover[dats2409.PC$Contents=="single_316Jeff"]

smpl.repro$c5.1 <- sample(repro_c5.1[!is.na(repro_c5.1)], reps, replace=TRUE)
smpl.repro$c5.2 <- sample(repro_c5.2[!is.na(repro_c5.2)], reps, replace=TRUE)
smpl.repro$c5.3 <- sample(repro_c5.3[!is.na(repro_c5.3)], reps, replace=TRUE)

repro_c5.expBS <- smpl.repro$c5.1/3 + smpl.repro$c5.2/3 + smpl.repro$c5.3/3
hist(repro_c5.expBS, breaks=50)
repro_c5.expBSmn <- mean(repro_c5.expBS)

repro_c5.devBS <- ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c5"] - repro_c5.expBSmn

t.test(ReproByPlot.PC$Repro[ReproByPlot.PC$Contents=="mix_c5"], repro_c5.expBS, paired="FALSE", conf.level=0.95,
       var.equal=TRUE, alternative="two.sided")






## Combined performance
CombAvgs.PC

## Calculate 'expected' repro for each regional mix based on additive of components
eComb_c1 <- (CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_LasAn"] + 
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_314Jeff"] +
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_316Jeff"])/3
eComb_c2 <- (CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_NM"] + 
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_294"] +
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_316Jeff"])/3
eComb_c3 <- (CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_UT"] + 
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_294"] +
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_314Jeff"])/3
eComb_c4 <- (CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_316Jeff"] + 
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_294"] +
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_314Jeff"])/3
eComb_c5 <- (CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_LasAn"] + 
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_294"] +
                CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="single_316Jeff"])/3

CombByPlot.PC <- dats2409.PC %>% group_by(Contents, Plot.Number) %>% 
  dplyr::summarise(Comb=mean(CombPerf,na.rm=TRUE))

## Calculate deviation from expected for each plot of each close mix
eComb_c1Dev <- CombByPlot.PC$Comb[CombByPlot.PC$Contents=="mix_c1"] - eComb_c1
eComb_c2Dev <- CombByPlot.PC$Comb[CombByPlot.PC$Contents=="mix_c2"] - eComb_c2
eComb_c3Dev <- CombByPlot.PC$Comb[CombByPlot.PC$Contents=="mix_c3"] - eComb_c3
eComb_c4Dev <- CombByPlot.PC$Comb[CombByPlot.PC$Contents=="mix_c4"] - eComb_c4
eComb_c5Dev <- CombByPlot.PC$Comb[CombByPlot.PC$Contents=="mix_c5"] - eComb_c5

## Combine diff mixes together and add values for x-axis for plotting
eCombDevs <- c(eComb_c1Dev, eComb_c2Dev, eComb_c3Dev, eComb_c4Dev, eComb_c5Dev)
xAx <- c(rep(1,length(eComb_c1Dev)), rep(2,length(eComb_c2Dev)), rep(3,length(eComb_c3Dev)),
         rep(4,length(eComb_c4Dev)), rep(5,length(eComb_c5Dev)))

eCombDevsForPlot <- as.data.frame(cbind(eCombDevs, xAx))

par(pty="s")
plot(eCombDevsForPlot$xAx, eCombDevsForPlot$eCombDevs, pch=16, ylim=c(-4,4),
     cex=0.9, ylab="Deviation from expected\ncombined performance", col="grey45",
     xlab="Regional Mix")
abline(h=0)


## Add plot means
## Calculate deviation from expected mean (avg of all plots) of each reg mix
eMnComb_c1Dev <- CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="mix_c1"] - eComb_c1
eMnComb_c2Dev <- CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="mix_c2"]  - eComb_c2
eMnComb_c3Dev <- CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="mix_c3"]  - eComb_c3
eMnComb_c4Dev <- CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="mix_c4"]  - eComb_c4
eMnComb_c5Dev <- CombAvgs.PC$Comb_MN[CombAvgs.PC$Contents=="mix_c5"]  - eComb_c5


## Combine diff mixes together and add values for x-axis for plotting
eMnCombDevs <- c(eMnComb_c1Dev, eMnComb_c2Dev, eMnComb_c3Dev, eMnComb_c4Dev, eMnComb_c5Dev)
xAx <- c(1:5)

eMnCombDevsForPlot <- as.data.frame(cbind(eMnCombDevs, xAx))

points(eMnCombDevsForPlot$xAx, eMnCombDevsForPlot$eMnCombDevs, pch=16,
       cex=2, col="black") 
## -----------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------






## Look for spatial patterns in performance data
plot(as.numeric(dats23.PC$Plot.Number), dats23.PC$Emrg)
plot(as.numeric(dats23.EP$Plot.Number), dats23.EP$Emrg)
plot(as.numeric(dats23.PC$Plot.Number), dats2406.PC$ARFR.Seedling.Count)
plot(as.numeric(dats23.EP$Plot.Number), dats2406.EP$ARFR.Seedling.Count)
plot(as.numeric(dats23.PC$Plot.Number), dats2406.PC$ARFR.Percent.Cover)
plot(as.numeric(dats23.EP$Plot.Number), dats2406.EP$ARFR.Percent.Cover)
plot(as.numeric(dats23.PC$Plot.Number), dats2409.PC$ARFR.Reproductive.Percent.Cover)
plot(as.numeric(dats23.EP$Plot.Number), dats2409.EP$ARFR.Reproductive.Percent.Cover)

boxplot(Emrg ~ HalfSplit, data=dats23.PC, ylim=c(0,0.2), 
        ylab="Emergence", cex.lab=1.5, main="Site 1")
boxplot(Emrg ~ HalfSplit, data=dats23.EP, ylim=c(0,0.1), 
        ylab="Emergence", cex.lab=1.5, main="Site 2")

boxplot(Surv ~ HalfSplit, data=dats2406.PC, ylim=c(0,1.6), 
        ylab="Survival", cex.lab=1.5, main="Site 1")
boxplot(Surv ~ HalfSplit, data=dats2406.EP, ylim=c(0,2), 
        ylab="Survival", cex.lab=1.5, main="Site 2")

boxplot(ARFR.Seedling.Count ~ HalfSplit, data=dats2406.PC, ylim=c(0,20), 
        ylab="Count", cex.lab=1.5, main="Site 1")
boxplot(ARFR.Seedling.Count ~ HalfSplit, data=dats2406.EP, ylim=c(0,10), 
        ylab="Count", cex.lab=1.5, main="Site 2")

boxplot(ARFR.Percent.Cover ~ HalfSplit, data=dats2406.PC, ylim=c(0,60), 
        ylab="Cover", cex.lab=1.5, main="Site 1")
boxplot(ARFR.Percent.Cover ~ HalfSplit, data=dats2406.EP, ylim=c(0,15), 
        ylab="Cover", cex.lab=1.5, main="Site 2")
modCovr <- lm(ARFR.Percent.Cover ~ HalfSplit, data=dats2406.PC)
summary(modCovr)
modCovr.EP <- lm(ARFR.Percent.Cover ~ HalfSplit, data=dats2406.EP)
summary(modCovr.EP)

boxplot(ARFR.Reproductive.Percent.Cover ~ HalfSplit, data=dats2409.PC, ylim=c(0,80), 
        ylab="Reproduction", cex.lab=1.5)
boxplot(ARFR.Reproductive.Percent.Cover ~ HalfSplit, data=dats2409.EP, ylim=c(0,25), 
        ylab="Reproduction", cex.lab=1.5)
modRepro.PC <- lm(ARFR.Reproductive.Percent.Cover ~ HalfSplit, data=dats2409.PC)
summary(modRepro.PC)

boxplot(CombPerf ~ HalfSplit, data=dats2409.PC, ylim=c(0,6), 
        ylab="Combined performance", cex.lab=1.5)
boxplot(CombPerf ~ HalfSplit, data=dats2409.EP, ylim=c(0,0.25), 
        ylab="Combined performance", cex.lab=1.5)
modComb.EP <- lm(CombPerf ~ HalfSplit, data=dats2409.EP)
summary(modComb.EP)



## Look at performance across environment quality ***
## Plot means at each site
#CombAvgs.EP <- dats2409.EP %>% group_by(Contents) %>% 
#  dplyr::summarise(Comb_MD=median(CombPerf,na.rm=TRUE),
#                   Comb_MN=mean(CombPerf,na.rm=TRUE),
#                   Comb_SE=calcSE(CombPerf))
#CombAvgs.EP$MixCol[grepl("single", CombAvgs.EP$Contents)] = "plum4"    
#CombAvgs.EP$MixCol[grepl("mix_c", CombAvgs.EP$Contents)] = "grey80"     
#CombAvgs.EP$MixCol[grepl("mix_b", CombAvgs.EP$Contents)] = "grey30"  

#plot(1:2, 0:1, type="n",ylim=c(0,1.3), xaxt="n", ylab="Combined performance", xlab="Site", cex.lab=1.5)
#for (ee in 1:length(unique(dats$Contents))) { 
#  points(1:2, c(CombAvgs.PC$Comb_MN[ee],CombAvgs.EP$Comb_MN[ee]), col=CombAvgs.PC$MixCol[ee], pch=19, cex=1.5)
#  lines(1:2, c(CombAvgs.PC$Comb_MN[ee],CombAvgs.EP$Comb_MN[ee]), col=CombAvgs.PC$MixCol[ee], pch=19, cex=1.5)
 #}
## ****
## ---------------------------------------------------------------------







## LOOK AT COEFFICIENT OF VARIATION ------------------------------------
## Calculate the coefficient of variation 
dats23PC.cv <- dats23.PC %>% group_by(Contents) %>% summarise(Emrg_CV=cv(Emrg, na.rm=TRUE))
dats23EP.cv <- dats23.EP %>% group_by(Contents) %>% summarise(Emrg_CV=cv(Emrg, na.rm=TRUE))
dats23.cv <- dats23 %>% group_by(Contents) %>% summarise(Emrg_CV=cv(Emrg, na.rm=TRUE))

dats23PC.cv$MixCol[grepl("single", dats23PC.cv$Contents)] = "plum4"    
dats23PC.cv$MixCol[grepl("mix_c", dats23PC.cv$Contents)] = "grey80"    
dats23PC.cv$MixCol[grepl("mix_b", dats23PC.cv$Contents)] = "grey40"    

dats23.cv$MixCol[grepl("single", dats23.cv$Contents)] = "plum4"    
dats23.cv$MixCol[grepl("mix_c", dats23.cv$Contents)] = "grey80"    
dats23.cv$MixCol[grepl("mix_b", dats23.cv$Contents)] = "grey40"    


par(mar=c(7,5,3,2))
barplot(dats23PC.cv$Emrg_CV, xlab=NA, ylab="Coefficient of variation in emergence", cex.lab=1.25, las=2, 
        names=dats23PC.cv$Contents, main="Site 1", cex.main=1.5, col=dats23PC.cv$MixCol, cex.names=1)

par(mar=c(7,5,3,2))
barplot(dats23.cv$Emrg_CV, xlab=NA, ylab="Coefficient of variation in emergence", cex.lab=1.25, las=2, 
        names=dats23.cv$Contents, main="Both sites", cex.main=1.5, col=dats23.cv$MixCol, cex.names=1)



## 2024 June
dats2406PC.cv <- dats2406.PC %>% group_by(Contents) %>% summarise(Surv_CV=cv(Surv, na.rm=TRUE),
                                                                        Counts_CV=cv(ARFR.Seedling.Count, na.rm=TRUE),
                                                                        Covr_CV=cv(ARFR.Percent.Cover, na.rm=TRUE))
dats2406EP.cv <- dats2406.EP %>% group_by(Contents) %>% summarise(Surv_CV=cv(Surv, na.rm=TRUE),
                                                                  Counts_CV=cv(ARFR.Seedling.Count, na.rm=TRUE),
                                                                  Covr_CV=cv(ARFR.Percent.Cover, na.rm=TRUE))
dats2406.cv <- dats2406 %>% group_by(Contents) %>% summarise(Surv_CV=cv(Surv, na.rm=TRUE),
                                                                  Counts_CV=cv(ARFR.Seedling.Count, na.rm=TRUE),
                                                                  Covr_CV=cv(ARFR.Percent.Cover, na.rm=TRUE))

dats2406PC.cv$MixCol[grepl("single", dats2406PC.cv$Contents)] = "plum4"    
dats2406PC.cv$MixCol[grepl("mix_c", dats2406PC.cv$Contents)] = "grey80"    
dats2406PC.cv$MixCol[grepl("mix_b", dats2406PC.cv$Contents)] = "grey40"    

dats2406EP.cv$MixCol[grepl("single", dats2406EP.cv$Contents)] = "plum4"    
dats2406EP.cv$MixCol[grepl("mix_c", dats2406EP.cv$Contents)] = "grey80"    
dats2406EP.cv$MixCol[grepl("mix_b", dats2406EP.cv$Contents)] = "grey40"    

dats2406.cv$MixCol[grepl("single", dats2406.cv$Contents)] = "plum4"    
dats2406.cv$MixCol[grepl("mix_c", dats2406.cv$Contents)] = "grey80"    
dats2406.cv$MixCol[grepl("mix_b", dats2406.cv$Contents)] = "grey40"    

barplot(dats2406PC.cv$Surv_CV, xlab=NA, ylab="Coefficient of variation in survival", cex.lab=1.25, las=2, 
        names=dats2406PC.cv$Contents, main="Site 1", cex.main=1.5, col=dats2406PC.cv$MixCol, cex.names=1)

barplot(dats2406PC.cv$Counts_CV, xlab=NA, ylab="Coefficient of variation in counts", cex.lab=1.25, las=2, 
        names=dats2406PC.cv$Contents, main="Site 1", cex.main=1.5, col=dats2406PC.cv$MixCol, cex.names=1)

barplot(dats2406PC.cv$Covr_CV, xlab=NA, ylab="Coefficient of variation in cover", cex.lab=1.25, las=2, 
        names=dats2406PC.cv$Contents, main="Site 1", cex.main=1.5, col=dats2406PC.cv$MixCol, cex.names=1)

barplot(dats2406.cv$Counts_CV, xlab=NA, ylab="Coefficient of variation in counts", cex.lab=1.25, las=2, 
        names=dats2406.cv$Contents, main="Both sites", cex.main=1.5, col=dats2406.cv$MixCol, cex.names=1)

barplot(dats2406.cv$Covr_CV, xlab=NA, ylab="Coefficient of variation in cover", cex.lab=1.25, las=2, 
        names=dats2406.cv$Contents, main="Both sites", cex.main=1.5, col=dats2406.cv$MixCol, cex.names=1)


#barplot(dats2406EP.cv$Surv_CV, xlab=NA, ylab="Coefficient of variation uin survival", cex.lab=1.25, las=2, 
#        names=dats2406EP.cv$Contents, main="Site 2", cex.main=1.5, col=dats2406EP.cv$MixCol, cex.names=1)

barplot(dats2406EP.cv$Counts_CV, xlab=NA, ylab="Coefficient of variation in counts", cex.lab=1.25, las=2, 
        names=dats2406EP.cv$Contents, main="Site 2", cex.main=1.5, col=dats2406EP.cv$MixCol, cex.names=1)

barplot(dats2406EP.cv$Covr_CV, xlab=NA, ylab="Coefficient of variation in percent cover", cex.lab=1.25, las=2, 
        names=dats2406EP.cv$Contents, main="Site 2", cex.main=1.5, col=dats2406EP.cv$MixCol, cex.names=1)




## 2024 September
dats2409PC.cv <- dats2409.PC %>% group_by(Contents) %>% summarise(Repro_CV=cv(ARFR.Reproductive.Percent.Cover, na.rm=TRUE))
dats2409EP.cv <- dats2409.EP %>% group_by(Contents) %>% summarise(Repro_CV=cv(ARFR.Reproductive.Percent.Cover, na.rm=TRUE))
dats2409.cv <- dats2409 %>% group_by(Contents) %>% summarise(Repro_CV=cv(ARFR.Reproductive.Percent.Cover, na.rm=TRUE))

dats2409PC.cv$MixCol[grepl("single", dats2409PC.cv$Contents)] = "plum4"    
dats2409PC.cv$MixCol[grepl("mix_c", dats2409PC.cv$Contents)] = "grey80"    
dats2409PC.cv$MixCol[grepl("mix_b", dats2409PC.cv$Contents)] = "grey40"    

dats2409EP.cv$MixCol[grepl("single", dats2409EP.cv$Contents)] = "plum4"    
dats2409EP.cv$MixCol[grepl("mix_c", dats2409EP.cv$Contents)] = "grey80"    
dats2409EP.cv$MixCol[grepl("mix_b", dats2409EP.cv$Contents)] = "grey40"    

dats2409.cv$MixCol[grepl("single", dats2409.cv$Contents)] = "plum4"    
dats2409.cv$MixCol[grepl("mix_c", dats2409.cv$Contents)] = "grey80"    
dats2409.cv$MixCol[grepl("mix_b", dats2409.cv$Contents)] = "grey40"    

barplot(dats2409PC.cv$Repro_CV, xlab=NA, ylab="Coefficient of variation in reproduction", cex.lab=1.25, las=2, 
        names=dats2409PC.cv$Contents, main="Site 1", cex.main=1.5, col=dats2409PC.cv$MixCol, cex.names=1)

#barplot(dats2409EP.cv$Repro_CV, xlab=NA, ylab="Coefficient of variation in reproduction", cex.lab=1.25, las=2, 
#        names=dats2409EP.cv$Contents, main="Site 2", cex.main=1.5, col=dats2409EP.cv$MixCol, cex.names=1)

barplot(dats2409.cv$Repro_CV, xlab=NA, ylab="Coefficient of variation in reproduction", cex.lab=1.25, las=2, 
        names=dats2409.cv$Contents, main="Both sites", cex.main=1.5, col=dats2409.cv$MixCol, cex.names=1)
## ---------










## PLOT RAW DATA AS MIX VS SINGLE ----------------------------------------------------------------------------------
Emrg_rates <- dats %>% group_by(Site, Contents) %>% summarise(RATE=ARFR.Seedling.Count/NumSeeds)

## Subset by mix/ single type, order lowest to highest emrg rate
# Site PC only 
Emrg_mixS1 <- subset(Emrg_rates, Contents=="mix_s1" & Site=="PC")
Emrg_mixS1 <- Emrg_mixS1[order(Emrg_mixS1$RATE),]
Emrg_mixS2 <- subset(Emrg_rates, Contents=="mix_s2" & Site=="PC")
Emrg_mixS2 <- Emrg_mixS2[order(Emrg_mixS2$RATE),]
Emrg_mixS3 <- subset(Emrg_rates, Contents=="mix_s3" & Site=="PC")
Emrg_mixS3 <- Emrg_mixS3[order(Emrg_mixS3$RATE),]
Emrg_mixS4 <- subset(Emrg_rates, Contents=="mix_s4" & Site=="PC")
Emrg_mixS4 <- Emrg_mixS4[order(Emrg_mixS4$RATE),]
Emrg_mixS5 <- subset(Emrg_rates, Contents=="mix_s5" & Site=="PC")
Emrg_mixS5 <- Emrg_mixS5[order(Emrg_mixS5$RATE),]

Emrg_UT <- subset(Emrg_rates, Contents=="single_UT" & Site=="PC")
Emrg_UT <- Emrg_UT[order(Emrg_UT$RATE),]
Emrg_294 <- subset(Emrg_rates, Contents=="single_294" & Site=="PC")
Emrg_294 <- Emrg_294[order(Emrg_294$RATE),]
Emrg_NM <- subset(Emrg_rates, Contents=="single_NM" & Site=="PC")
Emrg_NM <- Emrg_NM[order(Emrg_NM$RATE),]
Emrg_316Jeff <- subset(Emrg_rates, Contents=="single_316Jeff" & Site=="PC")
Emrg_316Jeff <- Emrg_316Jeff[order(Emrg_316Jeff$RATE),]
Emrg_314Jeff <- subset(Emrg_rates, Contents=="single_314Jeff" & Site=="PC")
Emrg_314Jeff <- Emrg_314Jeff[order(Emrg_314Jeff$RATE),]
Emrg_LasAn <- subset(Emrg_rates, Contents=="single_LasAn" & Site=="PC")
Emrg_LasAn <- Emrg_LasAn[order(Emrg_LasAn$RATE),]

## List what components are in each mix
#mix_s1: LasAn, 314Jeff, 316Jeff
#mix_s2: NM, 294, 316Jeff
#mix_s3: UT, 294, 314Jeff
#mix_s4: 294, 314Jeff, 316Jeff
#mix_s5: LasAn, 294, 316Jeff


## 

## Plot
par(pty="s")

plot((subset(Emrg_rates, Contents=="mix_s1" & Site=="PC"))$RATE,
     (subset(Emrg_rates, Contents=="single_316Jeff" & Site=="PC"))$RATE)

mean(Emrg_316Jeff$RATE)
mean(Emrg_314Jeff$RATE)
mean(Emrg_LasAn$RATE)
mean(Emrg_NM$RATE)
mean(Emrg_UT$RATE)
mean(Emrg_294$RATE)


plot(Emrg_316Jeff$RATE, Emrg_mixS1$RATE, ylim=c(0,0.11), xlim=c(0,0.11))
abline(0,1)

#plot(Emrg_316Jeff$RATE, Emrg_mixS2$RATE, ylim=c(0,0.11), xlim=c(0,0.11))
#abline(0,1)

#plot(Emrg_UT$RATE, Emrg_mixS3$RATE, ylim=c(0,0.11), xlim=c(0,0.11))
#abline(0,1)

plot(Emrg_316Jeff$RATE, Emrg_mixS4$RATE, ylim=c(0,0.12), xlim=c(0,0.12))
abline(0,1)

plot(Emrg_316Jeff$RATE, Emrg_mixS5$RATE, ylim=c(0,0.12), xlim=c(0,0.12))
abline(0,1)

#How to deal with diff number of plots per mix/ single? For avgs of singles that make up mixes, add NAs.
#or randomly drop plots so num remaining equals min number. Do the same when num plots per mix != num plots per single.
#Even with this, I don't know if these comparisons are appropriate. And if ordering by ascending val is ok.
## --------------------------------------------------------------------------------------------






## MODEL DATA -------------------------------------------------------------------------------
mod1 <- glmer(cbind(ARFR.Seedling.Count, NumSeeds-ARFR.Seedling.Count) ~ Contents + (1|Site),
              data=dats, family=binomial(link="logit"))
summary(mod1) #Not enough levels to use site as a random effect


mod2 <- glm(cbind(ARFR.Seedling.Count, NumSeeds-ARFR.Seedling.Count) ~ Contents + Site,
            data=dats, family=binomial(link="logit"))
summary(mod2)

dats <- dats %>% mutate(SeedMix = relevel(SeedMix, ref="Single")) #Make Single source the ref level
mod3 <- glm(cbind(ARFR.Seedling.Count, NumSeeds-ARFR.Seedling.Count) ~ SeedMix + Site,
            data=dats, family=binomial(link="logit"))
summary(mod3)
Anova(mod3)


#Look at each site separately
PC <- dats %>% dplyr::filter(Site=="PC")
mod4 <- glm(cbind(ARFR.Seedling.Count, NumSeeds-ARFR.Seedling.Count) ~ SeedMix,
            data=PC, family=binomial(link="logit"))
summary(mod4)
Anova(mod4)

EP <- dats %>% dplyr::filter(Site=="EP")
mod5 <- glm(cbind(ARFR.Seedling.Count, NumSeeds-ARFR.Seedling.Count) ~ SeedMix,
            data=EP, family=binomial(link="logit"))
summary(mod5)
Anova(mod5)
