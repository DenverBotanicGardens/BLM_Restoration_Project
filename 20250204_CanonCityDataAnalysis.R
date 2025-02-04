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

dats <- dats[dats$Treatment != "skip",]
dats <- dats[dats$Treatment != "pilot",]
#droplevels function

## Add column for number of seeds planted
dats$NumSeeds <- 400

## Add column with mix type or single 
dats$SeedMix[grepl("mix_d", dats$Contents)] <- "Mix_B"  #Broad mixes
dats$SeedMix[grepl("mix_s", dats$Contents)] <- "Mix_C"  #Close mixes
dats$SeedMix[grepl("single_", dats$Contents)] <- "Single"
dats$SeedMix <- as.factor(dats$SeedMix)


## Data checks
#seedling count should only be numeric
#seedling count on 20240918 should only be 1 or 0
#other checks..?

## Check if any ARFR seen in negative controls?
unique(dats$ARFR.Seedling.Count[dats$Treatment=="17"])
dats[dats$Treatment=="17" & dats$ARFR.Seedling.Count==1,]
##There is one incidence of a single ARFR seedling being counted in one survey in a single control plot
## Remove negative control rows
dats <- dats[dats$Treatment != "17",] #Negative controls
## ------------------------------------------------------------------------------------------------




## PLOT RAW DATA AS BOX PLOTS ----------------------------------------------------------------------------------
#cols <- c(rep("darkseagreen4", 10), rep("red4",10), rep("steelblue4",12))
#par(las=2)
#boxplot(ARFR.Seedling.Count ~ Site + Contents, data=dats, col=cols, ylim=c(0,80),
#        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32,34:35,37:38,40:41,43:44,46:47))

## 2023 emergence
dats23 <- dats[grepl("2023", dats$Date),]
cols <- c(rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6),
          rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6))
boxplot((ARFR.Seedling.Count/NumSeeds) ~ Contents + Site, data=dats23, col=cols, ylim=c(0,0.14),
        xlab=NA, ylab="Emergence Rate", xaxt="n", cex.lab=1.5)


boxplot(ARFR.Seedling.Count ~ SeedMix + Site, data=dats23)


EmrgByMed <- with(dats23, reorder(Contents, ARFR.Seedling.Count, median, na.rm=TRUE))
boxplot(ARFR.Seedling.Count ~ EmrgByMed, data=dats23,
        xlab="Seed composition", ylab="Emergence rate", cex.lab=1.25)


EmrgAvgs <- dats23 %>% group_by(Contents) %>% 
  dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_SE=calcSE(ARFR.Seedling.Count))

#EmrgAvgs <- left_join(EmrgAvgs, dats, by="Contents") 
EmrgAvgs$PopCol <- NA
EmrgAvgs$PopCol[grepl("mix_s1", EmrgAvgs$Contents)] = "red4"  
EmrgAvgs$PopCol[grepl("mix_s2", EmrgAvgs$Contents)] = "red4"     
EmrgAvgs$PopCol[grepl("mix_s3", EmrgAvgs$Contents)] = "red4"        
EmrgAvgs$PopCol[grepl("mix_s4", EmrgAvgs$Contents)] = "red4"        
EmrgAvgs$PopCol[grepl("mix_s5", EmrgAvgs$Contents)] = "red4"        
EmrgAvgs$PopCol[grepl("mix_d1", EmrgAvgs$Contents)] = "darkseagreen4"         
EmrgAvgs$PopCol[grepl("mix_d2", EmrgAvgs$Contents)] = "darkseagreen4"         
EmrgAvgs$PopCol[grepl("mix_d3", EmrgAvgs$Contents)] = "darkseagreen4"         
EmrgAvgs$PopCol[grepl("mix_d4", EmrgAvgs$Contents)] = "darkseagreen4"         
EmrgAvgs$PopCol[grepl("mix_d5", EmrgAvgs$Contents)] = "darkseagreen4" 
EmrgAvgs$PopCol[grepl("single_UT", EmrgAvgs$Contents)] = "steelblue4" 
EmrgAvgs$PopCol[grepl("single_294", EmrgAvgs$Contents)] = "steelblue4" 
EmrgAvgs$PopCol[grepl("single_NM", EmrgAvgs$Contents)] = "steelblue4" 
EmrgAvgs$PopCol[grepl("single_316Jeff", EmrgAvgs$Contents)] = "steelblue4" 
EmrgAvgs$PopCol[grepl("single_314Jeff", EmrgAvgs$Contents)] = "steelblue4" 
EmrgAvgs$PopCol[grepl("single_LasAn", EmrgAvgs$Contents)] = "steelblue4" 

plot(1:16, EmrgAvgs$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)


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

## Counts
cols <- c(rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6),
          rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6))
boxplot((ARFR.Seedling.Count) ~ Contents + Site, data=dats2406, col=cols, ylim=c(0,20),
        xlab=NA, ylab="Number of plants per plot", xaxt="n", cex.lab=1.5)

boxplot(ARFR.Seedling.Count ~ SeedMix + Site, data=dats2406)

## ** Look at total counts ** 
## ** Look at sites separately **

## Percent vegetative cover
boxplot((ARFR.Percent.Cover) ~ Contents + Site, data=dats2406, col=cols, ylim=c(0,60),
        xlab=NA, ylab="Percent cover/plot", xaxt="n", cex.lab=1.5)

boxplot(ARFR.Percent.Cover ~ SeedMix + Site, data=dats2406)


## Percent reproductive cover
dats2409 <- dats24[dats24$Date=="9/18/2024 18:00",]
boxplot((ARFR.Reproductive.Percent.Cover) ~ Contents + Site, data=dats2409, col=cols, ylim=c(0,80),
        xlab=NA, ylab="Percent reproductive cover/plot", xaxt="n", cex.lab=1.5)

boxplot(ARFR.Reproductive.Percent.Cover ~ SeedMix + Site, data=dats2409)








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





## TRY PLOTTING DEVIATION FROM EXPECTED ------------------------------------------------------
## Calculate 'expected' num emerged for each mix based on additive of components
e_s1 <- ((mean(Emrg_LasAn$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
e_s2 <- ((mean(Emrg_NM$RATE)*(400/3)) + (mean(Emrg_294$RATE)*(400/3)) 
         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
e_s3 <- ((mean(Emrg_UT$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
         + (mean(Emrg_294$RATE)*(400/3)))/400
e_s4 <- ((mean(Emrg_294$RATE)*(400/3)) + (mean(Emrg_314Jeff$RATE)*(400/3)) 
         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400
e_s5 <- ((mean(Emrg_LasAn$RATE)*(400/3)) + (mean(Emrg_294$RATE)*(400/3)) 
         + (mean(Emrg_316Jeff$RATE)*(400/3)))/400

## Calculate deviation from expected for each plot of each close mix
e_s1Dev <- Emrg_mixS1$RATE - e_s1
e_s2Dev <- Emrg_mixS2$RATE - e_s2
e_s3Dev <- Emrg_mixS3$RATE - e_s3
e_s4Dev <- Emrg_mixS4$RATE - e_s4
e_s5Dev <- Emrg_mixS5$RATE - e_s5

plot(rep(1,length(e_s5Dev)),e_s5Dev,)
abline(h=0)

## Combine diff mixes together and add values for x-axis for plotting
e_Devs <- c(e_s1Dev, e_s2Dev, e_s3Dev, e_s4Dev, e_s5Dev)
xAx <- c(rep(1,length(e_s1Dev)), rep(2,length(e_s2Dev)), rep(3,length(e_s3Dev)),
         rep(4,length(e_s4Dev)), rep(5,length(e_s5Dev)))

e_DevsForPlot <- as.data.frame(cbind(e_Devs, xAx))

plot(e_DevsForPlot$xAx, e_DevsForPlot$e_Devs, pch=16, ylim=c(-0.06,0.06),
     cex=0.9, ylab="Deviation from expected\nemergence rate", col="black",
     xlab="Close Mix")
abline(h=0)








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
