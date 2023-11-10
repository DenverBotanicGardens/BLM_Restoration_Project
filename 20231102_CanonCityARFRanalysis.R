## April Goebl
## Script started 2023-11-02
## BLM Restoration project at Denver Botanic Gardens
## Analyze data from Canyon City genetic diversity experiment   




rm(list=ls())
dev.off()


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
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/CanyonCity")
## ------------------------------------------------------------------------------------------------



## LOAD DATA --------------------------------------------------------------------------------------
dats <- read.csv(file="20231102_CanonCityData_ARFR.csv", sep=",", header=TRUE, dec=".")
## ------------------------------------------------------------------------------------------------



## LOOK AT AND MODIFY STRUCTURE OF DATA -----------------------------------------------------------
str(dats)
dats$Plot.Number <- as.factor(dats$Plot.Number)
dats$PlotTagNum <- as.factor(dats$PlotTagNum)
dats$Site <- as.factor(dats$Site)
dats$Treatment <- as.factor(dats$Treatment)
dats$Contents <- as.factor(dats$Contents)
## ------------------------------------------------------------------------------------------------



## ADD AND REMOVE FROM DATAFRAME ------------------------------------------------------------------
unique(dats$Contents)
#Remove empty, control, and pilot categories 
dats <- dats[dats$Contents != "empty",]
dats <- dats[dats$Contents != "controlAndSave",]  #Or should we keep negative control? 
dats <- dats[dats$Contents != "Fre12_1200seeds",]
dats <- dats[dats$Contents != "Fre12_800seeds",]
dats <- dats[dats$Contents != "422Nav_1200seeds",]
dats <- dats[dats$Contents != "422Nav_800seeds",]
dats <- dats[dats$Contents != "WY71_ClearPlot",]
dats <- dats[dats$Contents != "ERNA_Rio",]
dats <- dats[dats$Contents != "ERNA_Rou",]
#droplevels function

#Add column for number of seeds planted
dats$NumSeeds <- 400

#Add column with mix type or single 
dats$SeedMix[grepl("mix_d", dats$Contents)] <- "Mix_D"
dats$SeedMix[grepl("mix_s", dats$Contents)] <- "Mix_S"
dats$SeedMix[grepl("single_", dats$Contents)] <- "Single"
dats$SeedMix <- as.factor(dats$SeedMix)
## ------------------------------------------------------------------------------------------------




## PLOT RAW DATA ----------------------------------------------------------------------------------
cols <- c(rep("darkseagreen4", 10), rep("red4",10), rep("steelblue4",12))
par(las=2)
boxplot(ARFR.Seedling.Count ~ Site + Contents, data=dats, col=cols, ylim=c(0,80),
        at=c(1:2,4:5,7:8,10:11,13:14,16:17,19:20,22:23,25:26,28:29,31:32,34:35,37:38,40:41,43:44,46:47))

cols <- c(rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6),
          rep("darkseagreen4", 5), rep("red4",5), rep("steelblue4",6))
boxplot(ARFR.Seedling.Count ~ Contents + Site, data=dats, col=cols, ylim=c(0,75))


boxplot(ARFR.Seedling.Count ~ SeedMix + Site, data=dats)


EmrgByMed <- with(dats, reorder(Contents, ARFR.Seedling.Count, median, na.rm=TRUE))
boxplot(ARFR.Seedling.Count ~ EmrgByMed, data=dats,
        xlab="Seed composition", ylab="Emergence rate", cex.lab=1.25)


EmrgAvgs <- dats %>% group_by(Contents) %>% 
            dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_SE=calcSE(ARFR.Seedling.Count))

#EmrgAvgs <- left_join(EmrgAvgs, dats, by="Contents") 
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


EmrgAvgs_PC <- PC %>% group_by(Contents) %>% 
  dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_SE=calcSE(ARFR.Seedling.Count))

plot(1:16, EmrgAvgs_PC$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)



EmrgAvgs_EP <- EP %>% group_by(Contents) %>% 
  dplyr::summarise(Emrg_MD=median(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_MN=mean(ARFR.Seedling.Count,na.rm=TRUE),
                   Emrg_SE=calcSE(ARFR.Seedling.Count))

plot(1:16, EmrgAvgs_EP$Emrg_MN, col=EmrgAvgs$PopCol, pch=19, cex=1.25)





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
