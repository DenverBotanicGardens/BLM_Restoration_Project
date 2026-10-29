## April Goebl
## Script started 2026-10-19 (modified from 20231215_ChatfieldCGdataAnalysis_ERNAarchive)
## BLM Restoration project at Denver Botanic Gardens
## Analyze ERNA (Ericameria nauseosa) data from Chatfield Farms Common Garden Experiment 


rm(list=ls())
#dev.off()


## LOAD PACKAGES AND FUNCTIONS -------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)
## -----------------------------------------------------------------------------



## SET WORKING DIRECTORY -------------------------------------------------------
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland")
## -----------------------------------------------------------------------------




## ADD AND COMBINE FILES FOR GWAS/ GEA -----------------------------------------
ERNApheno <- read.csv(file="AGoebl/ProjectWriteupsAndSlides/ERNA_ms/20241018_ChatfieldPhenotypes_2022.2023_ERNA.csv", sep=",", header=TRUE, dec=".")
ERNAfam <- read.csv(file="AGoebl/ProjectWriteupsAndSlides/ERNA_ms/ERNA_biallelic_individ60miss_hetexcess_depthmin3max15_sitemiss0.8_mq40.fam", sep="\t", header=FALSE)

## Add new ID column to match fam file
ERNApheno$IDnew <- paste0("AG_ERNA_", ERNApheno$ID)

## Add desired phenotype for relevant indivs to fam file
colnames(ERNAfam) <- c("IDnew","V2","V3","V4","V5","V6")
ERNAfam <- left_join(ERNAfam, (ERNApheno %>% select(IDnew,DaysToFlwr2023)), by="IDnew")
## Remove column with all -9
ERNAfam <- ERNAfam %>% select(1:5,7)
colnames(ERNAfam) <- NULL

## Save updated fam file
write.table(ERNAfam, file="AGoebl/ProjectWriteupsAndSlides/ERNA_ms/ERNA_biallelic_individ60miss_hetexcess_depthmin3max15_sitemiss0.8_mq40_DAYSTOFLWR.fam", 
            sep="\t", row.names=FALSE, quote=FALSE)
## -----------------------------------------------------------------------------


## LOOK AT FILES FROM GWAS/ GEA ------------------------------------------------
ERNAassoc <- read.csv(file="AGoebl/ProjectWriteupsAndSlides/ERNA_ms/ERNA_biallelic_individ60miss_hetexcess_depthmin3max15_sitemiss0.8_mq40_DaysToFlwr.ulmm.assoc.txt", sep="\t", header=TRUE)
hist(ERNAassoc$p_wald)
hist(-log10(ERNAassoc$p_wald))
max(-log10(ERNAassoc$p_wald))
length(unique(ERNAassoc$chr))

## Look at lowest p-vals
ERNAassoc.ord <- ERNAassoc[order(ERNAassoc$p_wald),] 
ERNAassoc.loP <- ERNAassoc[ERNAassoc$p_wald<=0.001,]
ERNAassoc.loP <- ERNAassoc[-log10(ERNAassoc$p_wald)>=4,]

## Plot p-values along genomic coordiates
chr.list <- unique(ERNAassoc.loP$chr)

par(mfrow=c(5,2))
for (cc in chr.list) {
  chr.subst <- subset(ERNAassoc, chr == cc)
  
  plot(chr.subst$ps, -log10(chr.subst$p_wald), pch=19,cex=0.8,main=cc,
       xlab="Position along scafffold", ylab="-log10(p-val)",ylim=c(0,7))
}
## -----------------------------------------------------------------------------







## LOAD DATA --------------------------------------------------------------------------------------
ERNA22 <- read.csv(file="Chatfield/2022_data/20230302_ChatfieldData2022_ERNA.csv", sep=",", header=TRUE, dec=".")
ERNA.SdZn <- read.csv(file="AGoebl/Seeds/20231215_ERNA_LatLongSdZn_hexcodes.csv", sep=",", header=TRUE, dec=".")
ERNA.biovar <- readRDS("AGoebl/Seeds/20240131_ERNA_BiovarsAvg1980_2021")
ERNA23 <- read.csv(file="Chatfield/2023_data/20241017_ChatfieldData2023_ERNA.csv", sep=",", header=TRUE, dec=".")
ERNA.bioVarCG <- readRDS("AGoebl/Seeds/20240207_Chatfield_Biovars2022")
ERNA.sla <- read.csv(file="Chatfield/2023_data/20241017_ChatfieldSLAdata2023_ERNA.csv", sep=",", header=TRUE, dec=".")

## ADD 2024, 2025, and 2026 data 
ERNA24 <- read.csv(file="Chatfield/2024_data/CommonGarden/ERNA/20250321_ChatfieldCGsurveyData2024_ERNA_amgUpdates.csv", sep=",", header=TRUE, dec=".")
ERNA25 <- read.csv(file="Chatfield/2025_data/20251009_ERNA_plant_1_amgUpdates.csv", sep=",", header=TRUE, dec=".")
ERNA26 <- read.csv(file="Chatfield/2026_data/20261001_ERNA_plant_1_amgUpdates.csv", sep=",", header=TRUE, dec=".")
## ----------------------------------------------------------------------------------------------




## ERNA - DATA CLEAN UP --------------------------------------------------------
str(ERNA22)
str(ERNA23)
str(ERNA24)
str(ERNA25)


ERNA22$Source <- as.factor(ERNA22$Source)
ERNA22$Block <- as.factor(ERNA22$Block)
ERNA22$ID <- as.factor(ERNA22$ID)
ERNA23$Source <- as.factor(ERNA23$Source)
ERNA23$Block <- as.factor(ERNA23$Block)
ERNA23$ID <- as.factor(ERNA23$ID)
ERNA23$Survival_20231013 <- as.integer(as.character(ERNA23$Survival_20231013))
ERNA23$Flowering_20231013 <- as.integer(as.character(ERNA23$Flowering_20231013))
ERNA24$ID <- as.factor(ERNA24$plantID)
ERNA25$ID <- as.factor(ERNA25$plantID)




## 2022 -----------------------------
## Note: Seedlings planted 4/20-4/28, first surv 5/18, first replacements 5/18-5/19 (replaced if dead & missing & source available)
## (not replaced if appeared dead-no green tissue- but not missing, since these might not yet be fully dead.)
## Second surv/sz 6/8-6/10, second replacements 6/13 (replaced if dead on both 5/18 and 6/8, not replaced if only dead 6/8).

## Uncertainty around replacements and not tracking which individuals were replaced on which dates,
## led to confusion around which data could be used reliably for individuals that died early and were either replaced or not.
## Given this, ignore plants that died (likely from transplant shock) after initial planting and were either replaced or not. 

## Look at how many individuals fall into early death, which of these had tissue collected for DNA, etc. 
length(ERNA22$Replaced_YorN[ERNA22$Replaced_YorN!=""])
nrow(ERNA22[ERNA22$Replaced_YorN!="" & ERNA22$TissueCollected=="Y",])
ERNA22[ERNA22$Replaced_YorN!="" & ERNA22$TissueCollected=="Y",]
ERNA22[ERNA22$Replaced_YorN=="N" & ERNA22$TissueCollected=="Y",]
length(ERNA22$TissueCollected[ERNA22$TissueCollected=="Y"])

ERNA22$ExcludeBcReplaceYorN[ERNA22$Replaced_YorN!=""] <- "Y"
testRm <- ERNA22[ERNA22$Replaced_YorN!="",] 
testKeep <- ERNA22[ERNA22$Replaced_YorN=="",] 

## Check that surv is only 1, 0 and maybe NA
ERNA22[(ERNA22$Survival_20220608 < 0 | ERNA22$Survival_20220608 > 1) & !is.na(ERNA22$Survival_20220608),]
ERNA22[(ERNA22$Survival_20220719 < 0 | ERNA22$Survival_20220719 > 1) & !is.na(ERNA22$Survival_20220719),]
ERNA22[(ERNA22$Survival_20220909 < 0 | ERNA22$Survival_20220909 > 1) & !is.na(ERNA22$Survival_20220909),]
ERNA22[(ERNA22$Survival_20221108 < 0 | ERNA22$Survival_20221108 > 1) & !is.na(ERNA22$Survival_20221108),]
## **check other dates? **

## Check that surv is only integers
ERNA22[ERNA22$Survival_20220608 - floor(ERNA22$Survival_20220608) != 0,]
ERNA22[ERNA22$Survival_20220719 - floor(ERNA22$Survival_20220719) != 0,] 
ERNA22[ERNA22$Survival_20221108 - floor(ERNA22$Survival_20221108) != 0,]

## Most indivs didn't flower in the first year, ignore phenology for 2022

## Check that length is only numeric
is.numeric(ERNA22$Length_cm_20220608)
is.numeric(ERNA22$Length_cm_20220719)
is.numeric(ERNA22$Length_cm_20220915)

## Check that if surv=0 for a given date, there are no phenology or height values for that date 
ERNA22$Length_cm_20220608[ERNA22$Survival_20220608 == 0]
ERNA22$Length_cm_20220719[ERNA22$Survival_20220719 == 0]
ERNA22$Phenology_20220909[ERNA22$Survival_20220909 == 0]
ERNA22$Length_cm_20220915[ERNA22$Survival_20220915 == 0]
ERNA22$Phenology_20221108[ERNA22$Survival_20221108 == 0]
## ----------------------------------




## 2023 -----------------------------
## Look at and address data from individuals following harvesting 
## If 'H', 'harvested', 'harvest', 'h.', 'coll AGB' in notes (e.g. 'H' in Notes_20231013) subsequent data should be NA
setdiff(ERNA23$Notes_20230811,ERNA23$Notes_20230811.1) #Check difference in these two cols
unique(ERNA23$Notes_20230607[ERNA23$Notes_20230607!=""])
unique(ERNA23$Notes_20230705[ERNA23$Notes_20230705!=""])
unique(ERNA23$Notes_20230728[ERNA23$Notes_20230728!=""])
unique(ERNA23$Notes_20230811[ERNA23$Notes_20230811!=""])     #'coll AGB'
unique(ERNA23$Notes_20230811.1[ERNA23$Notes_20230811.1!=""]) #'coll AGB'
unique(ERNA23$Notes_20230901 [ERNA23$Notes_20230901!=""])    #'harvest','harvest/seed','dying, H','harvested','h. no inf dying','seed, h.','h. dry','seeding h.','seed H.' 
unique(ERNA23$Notes_20230914[ERNA23$Notes_20230914!=""])     #'H'
unique(ERNA23$Notes_20231013[ERNA23$Notes_20231013!=""])     #'H','Not harvested','no data collected; H','H? Looks like sm cut stems' 

## How to deal with data from harvested indivs, after harvesting? 
## Make subsequent (after harvesting for a given indiv) data NA only after harvest date
## Also add a column with Harvested Yes/No (any date in 2023) to use in filtering of data from later years

ERNA23[grepl("coll", ERNA23$Notes_20230811.1),] 
#ERNA23[grepl("coll", ERNA23$Notes_20230811),] #Use Notes column 20230811.1 as it contains more than 20230811
#ERNA23$Survival_20230901[grepl("coll", ERNA23$Notes_20230811.1)] <- NA
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20230901, Survival_20230914, Survival_20231013), ~ if_else(Notes_20230811.1=="coll AGB", NA, .x)))
ERNA23 <- ERNA23 %>% mutate(across(c(Flowering_20230901, Flowering_20230914, Flowering_20231013), ~ if_else(Notes_20230811.1=="coll AGB", NA, .x)))

ERNA23[grepl("harvest", ERNA23$Notes_20230901),]          
ERNA23[grepl("h\\.", ERNA23$Notes_20230901),]          
ERNA23[grepl("H", ERNA23$Notes_20230901),]      
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20230914, Survival_20231013, Flowering_20230914, Flowering_20231013), 
                                   ~ if_else(str_starts(ERNA23$Notes_20230901,"harvest"), NA, .x)))
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20230914, Survival_20231013, Flowering_20230914, Flowering_20231013), 
                                   ~ if_else(str_detect(ERNA23$Notes_20230901,"h\\."), NA, .x)))
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20230914, Survival_20231013, Flowering_20230914, Flowering_20231013), 
                                   ~ if_else(str_detect(ERNA23$Notes_20230901,"H"), NA, .x)))

ERNA23[grepl("H", ERNA23$Notes_20230914),]  
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20231013, Flowering_20231013), ~ if_else(Notes_20230914=="H", NA, .x)))

ERNA23[grepl("H", ERNA23$Notes_20231013),] 
ERNA23 <- ERNA23 %>% mutate(across(c(Survival_20231013, Flowering_20231013), ~ if_else(Notes_20231013=="H", NA, .x)))


## Add a column with Harvested Yes/No (any date in 2023) to use in filtering of data from later years **
#ERNA23 <- ERNA23 %>% mutate(Harvested_2023=if_else(
#                     str_starts(Notes_20230811.1,"coll AGB") | str_starts(Notes_20230901,"harvest") |
#                     str_detect(Notes_20230901,"h\\.") | str_detect(Notes_20230901,"H") |
#                     str_detect(Notes_20230914,"H") | str_detect(Notes_20231013,"H"), "Y", "N"))
#ERNA23$Harvested_2023[!is.na(ERNA23test5$AboveGroundBiomass_g)] #Some harvested plants didn't have note added 
ERNA23 <- ERNA23 %>% mutate(Harvested_2023=if_else(!is.na(AboveGroundBiomass_g), "Y", "N"))
testKeep <- ERNA23[ERNA23$Harvested_2023=="N",] 


## Check that surv is only 1, 0 and maybe NA
ERNA23[(ERNA23$Survival_20230607 < 0 | ERNA23$Survival_20230607 > 1) & !is.na(ERNA23$Survival_20230607),]
ERNA23[(ERNA23$Survival_20230620 < 0 | ERNA23$Survival_20230620 > 1) & !is.na(ERNA23$Survival_20230620),]
ERNA23[(ERNA23$Survival_20230627 < 0 | ERNA23$Survival_20230627 > 1) & !is.na(ERNA23$Survival_20230627),]
ERNA23[(ERNA23$Survival_20230705 < 0 | ERNA23$Survival_20230705 > 1) & !is.na(ERNA23$Survival_20230705),]
ERNA23[(ERNA23$Survival_20230728 < 0 | ERNA23$Survival_20230728 > 1) & !is.na(ERNA23$Survival_20230728),]
ERNA23[(ERNA23$Survival_20230811 < 0 | ERNA23$Survival_20230811 > 1) & !is.na(ERNA23$Survival_20230811),]
ERNA23[(ERNA23$Survival_20230901 < 0 | ERNA23$Survival_20230901 > 1) & !is.na(ERNA23$Survival_20230901),]
ERNA23[(ERNA23$Survival_20230914 < 0 | ERNA23$Survival_20230914 > 1) & !is.na(ERNA23$Survival_20230914),]
ERNA23[(ERNA23$Survival_20231013 < 0 | ERNA23$Survival_20231013 > 1) & !is.na(ERNA23$Survival_20231013),]

#Check that surv is only integers
ERNA23[ERNA23$Survival_20230607 - floor(ERNA23$Survival_20230607) != 0,]
ERNA23[ERNA23$Survival_20230620 - floor(ERNA23$Survival_20230620) != 0,]

## Check that biomass is only numeric
is.numeric(ERNA23$AboveGroundBiomass_g)
is.numeric(ERNA23$ReproBiomass_g)

## Check that if surv=0 for a given date, there are no phenology values for that date
ERNA23$Flowering_20230607[ERNA23$Survival_20230607 == 0]
ERNA23$Flowering_20230620[ERNA23$Survival_20230620 == 0]
ERNA23$Flowering_20230627[ERNA23$Survival_20230627 == 0]
ERNA23$Flowering_20230705[ERNA23$Survival_20230705 == 0] <- NA
ERNA23$Flowering_20230728[ERNA23$Survival_20230728 == 0]
ERNA23$Flowering_20230811[ERNA23$Survival_20230811 == 0]
ERNA23$Flowering_20230901[ERNA23$Survival_20230901 == 0]
ERNA23$Flowering_20230914[ERNA23$Survival_20230914 == 0]
ERNA23$Flowering_20231013[ERNA23$Survival_20231013 == 0]

#If plt alive and 1 was not entered in flowering col, enter 0 (not NA)
ERNA23$Flowering_20230607[ERNA23$Survival_20230607==1 & is.na(ERNA23$Flowering_20230607)] <- 0
ERNA23$Flowering_20230620[ERNA23$Survival_20230620==1 & is.na(ERNA23$Flowering_20230620)] <- 0
ERNA23$Flowering_20230705[ERNA23$Survival_20230705==1 & is.na(ERNA23$Flowering_20230705)] <- 0
ERNA23$Flowering_20230901[ERNA23$Survival_20230901==1 & is.na(ERNA23$Flowering_20230901)] <- 0 
ERNA23$Flowering_20230914[ERNA23$Survival_20230914==1 & is.na(ERNA23$Flowering_20230914)] <- 0
ERNA23$Flowering_20231013[ERNA23$Survival_20231013==1 & is.na(ERNA23$Flowering_20231013)] <- 0

## ** Think about for repro, if didn't flower, should repro bm be zero or NA? ** 

## Rename columns with 'Flowering' to 'Phenology' to match other years
names(ERNA23) <- sub("Flowering_", "Phenology_", names(ERNA23))
## ---------------------------------------------------------------

## For all years, add checks on integers only, or 1/0/NA only, and numeric only for size, etc? (before combining)
## *****



## 2024 -----------------------------
setdiff(ERNA24$plantID,ERNA24$ID) #Check difference in these two cols
anyDuplicated(ERNA24)
sum(duplicated(ERNA24))

## Make an updated date column: remove time, save as date data format? 
ERNA24 <- ERNA24 %>% mutate(Date=as.Date(sub(" .*","", date), format="%m/%d/%Y"))

## Combine date columns separated by just a day or two, as these should represent one collection point in time
ERNA24 <- ERNA24 %>% mutate(Date = case_when(Date == as.Date("2024-06-27") ~ as.Date("2024-06-26"),
                                             Date == as.Date("2024-07-03") ~ as.Date("2024-07-01"),
                                             TRUE ~ Date), Date=format(Date, "%Y-%m-%d"))

## Change from long to wide format to make comparable to other years (modified from ChatGPT)
ERNA24wide <- select(ERNA24, -c("GlobalID","plantID","plantID_db","ParentGlobalID","CreationDate","Creator","EditDate","Editor")) %>% 
                pivot_wider(id_cols=c(ID,Plant.Position.In.Row), names_from=Date, values_from=c(survival,phenology,notes,lengthInCm),
                            names_glue="{.value}_{format(as.Date(Date), '%Y%m%d')}")

dim(ERNA24wide)
n_distinct(ERNA24wide$ID)
anyDuplicated(ERNA24wide$ID)

## Rename columns to match trait_date format from other years
names(ERNA24wide) <- sub("^survival_", "Survival_", names(ERNA24wide))
names(ERNA24wide) <- sub("^phenology_", "Phenology_", names(ERNA24wide))
names(ERNA24wide) <- sub("^lengthInCm_", "Length_", names(ERNA24wide))
names(ERNA24wide) <- sub("^notes_", "Notes_", names(ERNA24wide))

## Check that surv is only 1, 0 and maybe NA **
#ERNA25[(ERNA25$Survival_20251009 < 0 | ERNA25$Survival_20251009 > 1) & !is.na(ERNA25$Survival_20251009),]

## Check that phenology is only 1, 0 and maybe NA **
#ERNA25[(ERNA25$Phenology_20251009 < 0 | ERNA25$Phenology_20251009 > 1) & !is.na(ERNA25$Phenology_20251009),]

## Check that size is only numeric **
## Note: Size measured on a subset of plts on 10/10 (end of 2024 season)

## Check that if surv=0 for a given date, there are no phenology values for that date **
#ERNA25$Phenology_20251009[ERNA25$Survival_20251009 == 0]

## If plt alive and 1 was not entered in flowering col, enter 0 (not NA) **
#ERNA23$Flowering_20230607[ERNA23$Survival_20230607==1 & is.na(ERNA23$Flowering_20230607)] <- 0
## ---------------------------------------------------------------





## 2025 -----------------------------
ERNA25 <- ERNA25 %>% rename(Survival_20251009 = survival, Phenology_20251009 = phenology)

## No size data was collected in 2025

## Check that surv is only 1, 0 and maybe NA
ERNA25[(ERNA25$Survival_20251009 < 0 | ERNA25$Survival_20251009 > 1) & !is.na(ERNA25$Survival_20251009),]

#Check that phenology is only 1, 0 and maybe NA
ERNA25[(ERNA25$Phenology_20251009 < 0 | ERNA25$Phenology_20251009 > 1) & !is.na(ERNA25$Phenology_20251009),]

## Check that if surv=0 for a given date, there are no phenology values for that date
ERNA25$Phenology_20251009[ERNA25$Survival_20251009 == 0]

#If plt alive and 1 was not entered in flowering col, enter 0 (not NA)


## ** Look into who was removed for SER transplants (before combining?)



## 2026 -----------------------------
## ** rename phenology and survival columns to include year
## ** size in 2026
## ** Look into who was removed for SER transplants (before combining?)




## ** After combining data (here probably),
## don't use individuals that have "Y" or "N" in 2022 Replaced_YorN column b/c of confusion around timing and what it means for subsequent data **

## * Confirm who was harvested 2023-2026 (after combining)
## * Check for consistency in survival and flowering across the year (do this for 22-26)
# ** Check that phenology only increases or stays the same; once flowering=1, it shouldn't go back to 0 
# ** Check that once zero in surv on X/X or later, stays zero (if becomes 1 later, could be data entry error, or not depending on species)




## Checks 

## *** Need to work on this **** Look at ARFR and Katie's methods for PEVI ** 
## CONSOLIDATE SURVIVAL DATA
#ERNA.Surv <- ERNA22.cl %>% dplyr::select(c(starts_with("Survival_")))
#for (rr in 1:nrow(ERNA.Surv)) {
#  for (cc in 1:(ncol(ERNA.Surv)-1)) {
#    if(ERNA.Surv[rr,cc]==0 & ERNA.Surv[rr,cc+1]==1) {
#     ERNA.Surv$Check[rr] <- "Y"
#    }
#  }
#}
#ERNA.Surv %>% filter_all(any_vars(.==0))
#ERNA22.cl %>% filter_all(starts_with("Survival_")==0)
## ----------------------------------------------------------------------------------------------

## **********************************************************************************************






## ERNA - DATA MODS ------------------------------------
## Add Source column where source name format matches Source in main data frame
ERNA.SdZn$Source <- str_replace(ERNA.SdZn$Code, "10-SOS", "")
ERNA.biovar$Source <- str_replace(ERNA.biovar$Pop, "10-SOS", "")


## Edit column names for biovariables
colnames(ERNA.biovar) <- c("Pop","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                           "bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                           "bio17","bio18","bio19","Source")


## Add seed zone name abbreviation column
ERNA.SdZn$SdZnAbbrev[grepl("15 - 20 Deg. F. / 2 - 3", ERNA.SdZn$seed_zone)] = "humid, warm"  
ERNA.SdZn$SdZnAbbrev[grepl("0 - 5 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "semi-humid, v.cold" 
ERNA.SdZn$SdZnAbbrev[grepl("5 - 10 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "semi-humid, cold"
ERNA.SdZn$SdZnAbbrev[grepl("10 - 15 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "semi-humid, cool"
ERNA.SdZn$SdZnAbbrev[grepl("15 - 20 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "semi-humid, warm"
ERNA.SdZn$SdZnAbbrev[grepl("20 - 25 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "semi-humid, v.warm"
ERNA.SdZn$SdZnAbbrev[grepl("10 - 15 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "semi-arid, cool"
ERNA.SdZn$SdZnAbbrev[grepl("15 - 20 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "semi-arid, warm"
ERNA.SdZn$SdZnAbbrev[grepl("20 - 25 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "semi-arid, v.warm"
ERNA.SdZn$SdZnAbbrev[grepl("20 - 25 Deg. F. / 12 - 30", ERNA.SdZn$seed_zone)] = "arid, v.warm"

## Add column with seed zone or pop 'order'
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

ERNA.SdZn$SdZnOrderNum[grepl("15 - 20 Deg. F. / 2 - 3", ERNA.SdZn$seed_zone)] = "1"    #humid, warm  
ERNA.SdZn$SdZnOrderNum[grepl("0 - 5 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "2"      #semi-humid, v.cold 
ERNA.SdZn$SdZnOrderNum[grepl("5 - 10 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "3"     #semi-humid, cold
ERNA.SdZn$SdZnOrderNum[grepl("10 - 15 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "4"    #semi-humid, cool
ERNA.SdZn$SdZnOrderNum[grepl("15 - 20 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "5"    #semi-humid, warm
ERNA.SdZn$SdZnOrderNum[grepl("20 - 25 Deg. F. / 3 - 6", ERNA.SdZn$seed_zone)] = "6"    #semi-humid, v.warm
ERNA.SdZn$SdZnOrderNum[grepl("10 - 15 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "7"   #semi-arid, cool
ERNA.SdZn$SdZnOrderNum[grepl("15 - 20 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "8"   #semi-arid, warm
ERNA.SdZn$SdZnOrderNum[grepl("20 - 25 Deg. F. / 6 - 12", ERNA.SdZn$seed_zone)] = "9"   #semi-arid, v.warm
ERNA.SdZn$SdZnOrderNum[grepl("20 - 25 Deg. F. / 12 - 30", ERNA.SdZn$seed_zone)] = "10"  #arid, v.warm
## ----------------------------------------------------------------------------------------------




## ERNA - COMBINE DATA TYPES --------------------------------------------
ERNA.biovar <- left_join(ERNA.biovar, ERNA.SdZn, by="Source")
ERNA22.cl <- left_join(ERNA22.cl, ERNA.biovar, by="Source")
#ERNA23 <- left_join(ERNA23, ERNA.SdZn, by="Source")
ERNA23 <- left_join(ERNA23, ERNA.biovar, by="Source")
## ----------------------------------------------------------------------

## *********************



## ERNA 2023 - CREATE FLOWERING/ PHENOLOGY VARIABLE(S) ---------------------------------------------
## For 2023, estimate days to 1st flwr & if flowered at all based on '1' in any pheno survey (not just last as not always consistent)

#ERNA.Pheno <- ERNA23 %>% dplyr::select(starts_with("Flowering"))
#min(ERNA.Pheno, na.rm=TRUE)
#max(ERNA.Pheno, na.rm=TRUE)

## Calculate days to first flower 
ERNA.StartDate <- as.Date("2023-01-01")
ERNA.PhenoCol.List <- colnames(ERNA23)[grepl("Flowering*", colnames(ERNA23))]   #Obtain phenology column names
ERNA.Pheno.List <- str_replace(ERNA.PhenoCol.List, "Flowering_", "")            #Obtain just date from phenology columns
ERNA.Pheno.List <- as.Date(ERNA.Pheno.List, "%Y%m%d")
ERNA.DaysToFlwr <- ERNA.Pheno.List - ERNA.StartDate                             #Calculate days from Jan 1 to each phenology survey 

## Loop over each phenology column & enter the num days since Jan 1 when a 1 (bud or later repro stage) first appears
ERNA23$DaysToFlwr <- NA
for (pp in 1:length(ERNA.PhenoCol.List)) {
  ERNA23$DaysToFlwr[ERNA23[,ERNA.PhenoCol.List[pp]]==1 & is.na(ERNA23$DaysToFlwr)] <- as.integer(ERNA.DaysToFlwr)[pp]
}
ERNA23 %>% group_by(Source) %>% dplyr::summarise(Pheno_Avg=mean(DaysToFlwr,na.rm=TRUE))

## Make a flowered Yes or No column
ERNA23$FlwrYesNo <- NA 
ERNA23$FlwrYesNo[!is.na(ERNA23$DaysToFlwr)] <- 1                               #If there's a value in Days to Flwr, then enter Yes
ERNA23$FlwrYesNo[is.na(ERNA23$DaysToFlwr) & ERNA23$Survival_20231013==1] <- 0  #If plt alive and didn't flw, enter No
nrow(ERNA23[ERNA23$FlwrYesNo==0,])                                             #Num that survived but didn't flower 
#ERNA23 %>% group_by(Source) %>% dplyr::summarise(FlwrYesNo_Avg=mean(FlwrYesNo,na.rm=TRUE)) #Or could try sum and/or NUM=n()
## ---------------------------------------------------------------
## ---------------------------------------------------------------

## *** ESTIMATE PHENOLOGY FROM 2024 AS WELL ****



## ERNA - ADD GROWTH RATE VARIABLES ------------------------------
ERNA22.cl$GrwthRate_Specific <- log(ERNA22.cl$Length_cm_20220915/ERNA22.cl$Length_cm_20220719)
ERNA22.cl$GrwthRate_Absolute <- ERNA22.cl$Length_cm_20220915-ERNA22.cl$Length_cm_20220719
ERNA22.cl$GrwthRate_Relative <- (ERNA22.cl$Length_cm_20220915-ERNA22.cl$Length_cm_20220719)/ERNA22.cl$Length_cm_20220719

## ** Look at early vs late growth if at least 3 height measurements are usable ** 
ERNA22.cl$GrwthRateE_Specific <- log(ERNA22.cl$Length_cm_20220719/ERNA22.cl$Length_cm_20220608)
ERNA22.cl$GrwthRateE_Absolute <- ERNA22.cl$Length_cm_20220719-ERNA22.cl$Length_cm_20220608
ERNA22.cl$GrwthRateE_Relative <- (ERNA22.cl$Length_cm_20220719-ERNA22.cl$Length_cm_20220608)/ERNA22.cl$Length_cm_20220608

plot(ERNA22.cl$GrwthRate_Relative, ERNA22.cl$GrwthRateE_Relative)
## ---------------------------------------------------------------



## LOOK AT AGB FOR 2023
## Add any mods? 
ERNA23$AboveGroundBiomass_g
ERNA23$AboveGroundBiomass_g[!is.na(ERNA23$AboveGroundBiomass_g)] #Number of plts with bm data
length(ERNA23$AboveGroundBiomass_g[!is.na(ERNA23$AboveGroundBiomass_g)])
hist(ERNA23$AboveGroundBiomass_g)
str(ERNA23$AboveGroundBiomass_g)
## ---------------------------------------------------------------



## LOOK AT REPRO BM DATA FOR 2023
ERNA23$ReproBiomass_g
ERNA23$ReproBiomass_g[!is.na(ERNA23$ReproBiomass_g)] #Number of plts with bm data
length(ERNA23$ReproBiomass_g[!is.na(ERNA23$ReproBiomass_g)])
hist(ERNA23$ReproBiomass_g)
str(ERNA23$ReproBiomass_g)
## ---------------------------------------------------------------



## ERNA - ESTIMATE SURVIVAL 
## For 2023, estimate survival based on if alive at end of season (i.e. 20231013 survey)
ERNA23$AliveYesNo <- 0
ERNA23$AliveYesNo[ERNA23$Survival_20231013==1] <- 1
# ** Remove plants that died early in 2022 (e.g. from transplant shock) and were not replaced **

## 2022
ERNA22.cl$AliveYesNo <- 0
ERNA22.cl$AliveYesNo[ERNA22.cl$Survival_20221108==1] <- 1
ERNA22.cl$AliveYesNo[(ERNA22.cl$OrigPltSurvival_20220518==0) & (ERNA22.cl$Replaced_YorN=="N" | ERNA22.cl$Replaced_YorN=="")] <- NA
#ERNA.cl$AliveYesNo[(ERNA.cl$OrigPltSurvival_20220518==0 | ERNA.cl$Survival_20220608==0) & (ERNA.cl$Replaced_YorN=="N" | ERNA.cl$Replaced_YorN=="")] <- NA
#ERNA.cl %>% group_by(Source) %>% dplyr::summarise(AliveYesNo_Avg=mean(AliveYesNo,na.rm=TRUE))

## ** ESTIMATE SURVIVAL FOR 2023-2026 **
## ** Models with all years together and year as a fixed effect? **
## ** Try CoxHaz model for time to death estimates ** 
## ---------------------------------------------------------------



## ADD SOME 2023 DATA TO 2022 DATAFRAME TO ALLOW FOR COMPARISONS AND FILTERING 
ERNA23.sel <- ERNA23 %>% dplyr::select(c("ID","SentForSequencing_20230803","DaysToFlwr","AliveYesNo","FlwrYesNo",
                                         "AboveGroundBiomass_g","ReproBiomass_g")) 
ERNA.cl <- left_join(ERNA22.cl, ERNA23.sel, by="ID") 


## Add in SLA data
ERNA.cl <- left_join(ERNA.cl, ERNA.sla, by="ID")
str(ERNA.cl$SLA_mm2permg)
ERNA.cl$SLA_mm2permg <- as.numeric(as.character(ERNA.cl$SLA_mm2permg))
## ---------------------------------------------------------------






## MAKE FILE WITH JUST COLUMNS/ PHENOTYPES OF INTEREST FOR ANALYSIS WITH GENOTYPE DATA --------------
ERNA.clSel <- ERNA.cl %>% dplyr::select(c("Source","ID","SentForSequencing_20230803","Length_cm_20220608",
                                          "Length_cm_20220719","Length_cm_20220915","DaysToFlwr","SLA_mm2permg",
                                          "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
                                          "bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19",
                                          "Longitude","Latitude","HexCode","seed_zone","SdZnAbbrev","PopAbbrev")) 

write.csv(ERNA.clSel, "Chatfield/20241018_ChatfieldPhenotypes_ERNA.csv", row.names=FALSE)
## --------------------------------------------------------------------------------------------------

## ***********************************************






## ERNA - TEST FOR WATER TREATMENT EFFECT in 2022 ---------------------------------------------------
hist(log(ERNA.cl$Length_cm_20220915))
hist(ERNA.cl$Length_cm_20220915)
ERNA.tx.mod <- lmer(Length_cm_20220915 ~ Source + Treatment + (1|Block), data=ERNA.cl)
ERNA.pop.mod <- lmer(Length_cm_20220915 ~ Source + (1|Block), data=ERNA.cl)
models <- list(ERNA.tx.mod, ERNA.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No strong support for treatment (delta AICc very small)

hist(log(ERNA.cl$GrwthRate_Specific))
hist(ERNA.cl$GrwthRate_Specific)
hist(ERNA.cl$GrwthRate_Absolute)
hist(ERNA.cl$GrwthRate_Relative)
hist(log(ERNA.cl$GrwthRate_Relative))
ERNA.tx.mod <- lmer(GrwthRate_Specific ~ Source + Treatment + (1|Block), data=ERNA.cl)
ERNA.pop.mod <- lmer(GrwthRate_Specific ~ Source + (1|Block), data=ERNA.cl)
ERNA.tx.mod <- lmer(GrwthRate_Absolute ~ Source + Treatment + (1|Block), data=ERNA.cl)
ERNA.pop.mod <- lmer(GrwthRate_Absolute ~ Source + (1|Block), data=ERNA.cl)
ERNA.tx.mod <- lmer(log(GrwthRate_Relative) ~ Source + Treatment + (1|Block), data=ERNA.cl)
ERNA.pop.mod <- lmer(log(GrwthRate_Relative) ~ Source + (1|Block), data=ERNA.cl)
models <- list(ERNA.tx.mod, ERNA.pop.mod)
mod.names <- c('IncldTx', 'JustPop')
aictab(cand.set = models, modnames = mod.names )
# No support for treatment
## -----------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------





## ERNA - VISUALIZE RAW DATA ---------------------------------------------------------------

## Order populations for plotting 
## 2022 and 2023 - Order by average size or longitude or seed zone
#ERNA.htByMed <- with(ERNA.cl, reorder(Source, Length_cm_20220915, median, na.rm=TRUE))
ERNA.longByMed <- with(ERNA.cl, reorder(Source, Longitude, median, na.rm=TRUE))
str(ERNA.cl$SdZnOrderNum)
ERNA.cl$SdZnOrderNum <- as.integer(ERNA.cl$SdZnOrderNum)
ERNA.sdznByMed <- with(ERNA.cl, reorder(Source, as.integer(SdZnOrderNum), median, na.rm=TRUE))


ERNA.meds <- ERNA.cl %>% group_by(Source) %>% 
  dplyr::summarise(Height22_MD=median(Length_cm_20220915,na.rm=TRUE), GrowthAb_MD=median(GrwthRate_Absolute,na.rm=TRUE),
                   GrowthRe_MD=median(GrwthRate_Relative,na.rm=TRUE), GrowthSp_MD=median(GrwthRate_Specific,na.rm=TRUE),
                   GrowthAbE_MD=median(GrwthRateE_Absolute,na.rm=TRUE), GrowthReE_MD=median(GrwthRateE_Relative,na.rm=TRUE), 
                   GrowthSpE_MD=median(GrwthRateE_Specific,na.rm=TRUE), DaysToFlwr23_MD=median(DaysToFlwr,na.rm=TRUE),
                   AGB23_MD=median(AboveGroundBiomass_g,na.rm=TRUE), ReproBM23_MD=median(ReproBiomass_g,na.rm=TRUE))
ERNA.meds <- left_join(ERNA.meds, ERNA.SdZn, by="Source")


## 2023 - Order by time to flower 
#ERNA.dfByMed <- with(ERNA23, reorder(Source, DaysToFlwr, median, na.rm=TRUE))
#ERNA23.meds <- ERNA23 %>% group_by(Source) %>% 
#               dplyr::summarise(DaysToFlwr_MD=median(DaysToFlwr,na.rm=TRUE))
#ERNA23.meds <- left_join(ERNA23.meds, ERNA.SdZn, by="Source")


## Boxplots of raw data 
#ERNA.meds <- ERNA.meds[order(ERNA.meds$Longitude),] #Order by longitude
ERNA.meds$SdZnOrderNum <- as.integer(ERNA.meds$SdZnOrderNum)
str(ERNA.meds$SdZnOrderNum)
ERNA.meds <- ERNA.meds[order(ERNA.meds$SdZnOrderNum),] #Order by seed zone


par(mfrow=c(4,2))

## Size 2022
#ERNA.meds <- ERNA.meds[order(ERNA.meds$Height22_MD),] #Order by median 
boxplot(Length_cm_20220915 ~ ERNA.sdznByMed, data=ERNA.cl,
        xlab=NA, ylab="Height (cm)", cex.lab=1.25,
        cex.axis=0.99, names=ERNA.meds$PopAbbrev, las=2,
        main="FINAL SIZE 2022", cex.main=1.1, col=ERNA.meds$HexCode)

## Growth rate(s) 2022 ** Add time interval to growth rate calcs? **
#boxplot(GrwthRateE_Relative ~ ERNA.longByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant relative growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode, main="EARLY GROWTH RATE", ylim=c(-0.25,5.6))

boxplot(GrwthRateE_Absolute ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Plant absolute growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="EARLY SEASON GROWTH RATE 2022", ylim=c(-10,35))

#boxplot(GrwthRateE_Specific ~ ERNA.longByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant specific growth", cex.lab=1, cex.axis=0.9, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode, main="EARLY GROWTH RATE")


#ERNA.meds <- ERNA.meds[order(ERNA.meds$Growth_MD),]
#boxplot(GrwthRate_Relative ~ ERNA.longByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant relative growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode, main="GROWTH RATE", ylim=c(-0.35,3))

boxplot(GrwthRate_Absolute ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Plant absolute growth", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="LATER SEASON GROWTH RATE 2022")

#boxplot(GrwthRate_Specific ~ ERNA.longByMed, data=ERNA.cl, las=2,
#        xlab=NA, ylab="Plant specific growth", cex.lab=1, cex.axis=0.9, names=ERNA.meds$PopAbbrev,
#        cex.main=1.5, col=ERNA.meds$HexCode)

#plot.new()

#ERNA.SdZn <- ERNA.SdZn[order(ERNA.SdZn$GR_SP),]
#boxplot(GrwthRate_Specific ~ ERNA.grsByMed, data=ERNA.cl,
#        xlab="Population", ylab="Plant specific growth", cex.lab=1.5, names=ERNA.SdZn$PopAbbrev,
#        cex.axis=0.7, cex.main=1.5, col=ERNA.SdZn$HexCode)
#par(mfrow=c(1,1))


## Flowering time 2023
boxplot(DaysToFlwr ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Days to first flower", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="PHENOLOGY 2023")

## ABG 2023
boxplot(AboveGroundBiomass_g ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Above-ground biomass (g)", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="PLANT BIOMASS 2023", ylim=c(0,350))

## Repro BM 2023
boxplot(ReproBiomass_g ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Reproductive biomass (g)", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="REPRODUCTIVE OUTPUT 2023")

## SLA 2023
boxplot(SLA_mm2permg ~ ERNA.sdznByMed, data=ERNA.cl, las=2,
        xlab=NA, ylab="Specific leaf area (mm2/mg)", cex.lab=1.25, cex.axis=0.99, names=ERNA.meds$PopAbbrev,
        cex.main=1.1, col=ERNA.meds$HexCode, main="SPECIFIC LEAF AREA 2023", ylim=c(2,12))


## Blank plot
plot.new()
par(mfrow=c(1,1))
legend("center", unique(ERNA.meds$SdZnAbbrev[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), col="black",
       pt.bg=unique(ERNA.meds$HexCode[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), cex=1.65, pch=21)
## -------------------------


## 2023
#par(mfrow=c(1,2))

## Days to Flwr
#ERNA23.meds <- ERNA23.meds[order(ERNA23.meds$DaysToFlwr_MD),] #Order by median 
#boxplot(DaysToFlwr ~ ERNA.dfByMed, data=ERNA23,
#        xlab=NA, ylab="Days to first flower", cex.lab=1.25,
#        cex.axis=0.99, names=ERNA23.meds$PopAbbrev, las=2,
#        main="PHENOLOGY", cex.main=1.5, col=ERNA23.meds$HexCode)

## Blank plot
#plot.new()
#legend("center", unique(ERNA.meds$seed_zone[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), col="black",
#       pt.bg=unique(ERNA.meds$HexCode[order(ERNA.meds$SdZnOrder, decreasing=FALSE)]), cex=1.24, pch=21)
## -------------------------





## 2023 - Plots means as points and SE as error bars
ERNA23.mn <- ERNA23 %>% group_by(Source) %>% dplyr::summarise(AliveYesNo_MN=mean(AliveYesNo,na.rm=TRUE),
                                                              AliveYesNo_SE=calcSE(AliveYesNo), DaysToFlwr_MN=mean(DaysToFlwr,na.rm=TRUE),
                                                              DaysToFlwr_SE=calcSE(DaysToFlwr), FlwrYesNo_MN=mean(FlwrYesNo,na.rm=TRUE),
                                                              FlwrYesNo_SE=calcSE(FlwrYesNo), AliveYesNo_SD=sd(AliveYesNo,na.rm=TRUE)) #Or could try sum and/or NUM=n()

ERNA23.mn <- left_join(ERNA23.mn, ERNA.biovar, by="Source")

ERNA23.mn <- ERNA23.mn[order(ERNA23.mn$SdZnOrder),]
plot(c(1:20), ERNA23.mn$AliveYesNo_MN, col="black", pch=21, cex=1.5, ylab="Survival rate",xlab="Population",
     ylim=c(0.45,0.95), bg=ERNA23.mn$HexCode, xaxt="n",cex.lab=1.2)
arrows(c(1:20), ERNA23.mn$AliveYesNo_MN+ERNA23.mn$AliveYesNo_SE, c(1:20), ERNA23.mn$AliveYesNo_MN-ERNA23.mn$AliveYesNo_SE,
       angle=90, length=0, col="grey")

plot(c(1:20), ERNA23.mn$FlwrYesNo_MN, col="black", pch=21, cex=1.5, ylab="Flowering rate",xlab="Population",
     ylim=c(0.65,1.01), bg=ERNA23.mn$HexCode, xaxt="n",cex.lab=1.2)
arrows(c(1:20), ERNA23.mn$FlwrYesNo_MN+ERNA23.mn$FlwrYesNo_SE, c(1:20), ERNA23.mn$FlwrYesNo_MN-ERNA23.mn$FlwrYesNo_SE,
       angle=90, length=0, col="grey")

plot(c(1:20), ERNA23.mn$DaysToFlwr_MN, col="black", pch=21, cex=1.5, ylab="Days to first flower",xlab="Population",
     ylim=c(185,242), bg=ERNA23.mn$HexCode, xaxt="n",cex.lab=1.2)
arrows(c(1:20), ERNA23.mn$DaysToFlwr_MN+ERNA23.mn$DaysToFlwr_SE, c(1:20), ERNA23.mn$DaysToFlwr_MN-ERNA23.mn$DaysToFlwr_SE,
       angle=90, length=0, col="grey")

legend("bottonright", unique(ERNA23.mn$SdZnAbbrev[order(ERNA23.mn$SdZnOrder)]), 
       col=unique(BOGR.meds$SdZnColful[order(BOGR.meds$SdZnOrder)]), cex=1.95, pch=19)


## *****
## Filter 2023 survival based on 2022 data (i.e. did plt die early due to transplant shock) an re-plot means and sE ***
## *****


## 2023 - Bar plots of data (color by seed zone)
## Order populations for plotting BY VALUE
par(mfrow=c(1,1))
ERNA23.mn <- ERNA23.mn[order(ERNA23.mn$AliveYesNo_MN),]
barXvals<-barplot(ERNA23.mn$AliveYesNo_MN, xlab=NA, ylab="Survival rate", cex.lab=1.3, las=2, 
                  names=ERNA23.mn$PopAbbrev, main="Ericameria nauseosa", cex.main=1.5, col=ERNA23.mn$HexCode, cex.names=0.9)
arrows(barXvals, ERNA23.mn$AliveYesNo_MN+ERNA23.mn$AliveYesNo_SE, barXvals, ERNA23.mn$AliveYesNo_MN-ERNA23.mn$AliveYesNo_SE,
       angle=90, length=0, col="black")

ERNA23.mn <- ERNA23.mn[order(ERNA23.mn$FlwrYesNo_MN),]
barplot(ERNA23.mn$FlwrYesNo_MN, xlab=NA, ylab="Flowering rate", cex.lab=1.3, las=2, 
        names=ERNA23.mn$PopAbbrev, main="Ericameria nauseosa", cex.main=1.5, col=ERNA23.mn$HexCode, cex.names=0.9)
arrows(barXvals, ERNA23.mn$FlwrYesNo_MN+ERNA23.mn$FlwrYesNo_SE, barXvals, ERNA23.mn$FlwrYesNo_MN-ERNA23.mn$FlwrYesNo_SE,
       angle=90, length=0, col="black")

ERNA23.mn <- ERNA23.mn[order(ERNA23.mn$DaysToFlwr_MN),]
barplot(ERNA23.mn$DaysToFlwr_MN, xlab=NA, ylab="Days to first flower", cex.lab=1.3, las=2, ylim=c(0,250),
        names=ERNA23.mn$PopAbbrev, main="Ericameria nauseosa", cex.main=1.5, col=ERNA23.mn$HexCode, cex.names=0.9)
arrows(barXvals, ERNA23.mn$DaysToFlwr_MN+ERNA23.mn$DaysToFlwr_SE, barXvals, ERNA23.mn$DaysToFlwr_MN-ERNA23.mn$DaysToFlwr_SE,
       angle=90, length=0, col="black")
## ---------------------------------------------------




## ERNA - LOOK AT RELATIONSHIP BETWEEN TRAIT VAR AND AVERAGE ---------------------------
#par(mfrow=c(1,2))

#ERNA.ht.mn <- ERNA.cl %>% group_by(Source) %>% dplyr::summarise(Height_MN=mean(Length_cm_20220915, na.rm=TRUE))
#ERNA.gr.mn <- ERNA.cl %>% group_by(Source) %>% dplyr::summarise(Growth_MN=mean(GrwthRate_Relative, na.rm=TRUE))

#ERNA.ht.sd <- ERNA.cl %>% group_by(Source) %>% dplyr::summarise(Height_SD=sd(Length_cm_20220915, na.rm=TRUE))
#ERNA.gr.sd <- ERNA.cl %>% group_by(Source) %>% dplyr::summarise(Growth_SD=sd(GrwthRate_Relative, na.rm=TRUE))

#plot(ERNA.ht.sd$Height_SD, ERNA.ht.mn$Height_MN, pch=19, col="black", cex=1.3, 
#     xlab="Standard deviation of plant height", ylab="Mean plant height")
#plot(ERNA.gr.sd$Growth_SD, ERNA.gr.mn$Growth_MN, pch=19, col="black", cex=1.3,
#     xlab="Standard deviation of relative growth rate", ylab="Mean relative growth rate")
## ----------------------------------------------------------------------------------------




## ERNA - LOOK AT CORRELATION BETWEEN TRAITS ----------------------------------------------
#par(mfrow=c(1,1))
#plot(ERNA.ht.mn$Height_MN, ERNA.gr.mn$Growth_MN, pch=19, col="black", cex=1.3,
#     xlab="Mean plant height", ylab="Mean relative growth rate", cex.lab=1.2)

plot(ERNA23.mn$AliveYesNo_MN, ERNA23.mn$FlwrYesNo_MN, col="black", pch=21, cex=1.5,
     bg=ERNA23.mn$HexCode)

## Plot early vs late growth -----------
plot(ERNA.meds$GrowthRe_MD, ERNA.meds$GrowthReE_MD)
plot(ERNA.meds$GrowthAb_MD, ERNA.meds$GrowthAbE_MD)
plot(ERNA.meds$GrowthSp_MD, ERNA.meds$GrowthSpE_MD)
## ----------------------------------------------------------------------------------------




## ERNA - ESTIMATE VARIATION WITHIN POPULATIONS ------------------------------------------------------
## Calculate the coefficient of variation 
ERNA.cv <- ERNA.cl %>% group_by(Source) %>% summarise(Height_CV=cv(Length_cm_20220915, na.rm=TRUE),
                                                      Growth_CV=cv(GrwthRate_Relative, na.rm=TRUE),
                                                      DaysToFlwr_CV=cv(DaysToFlwr, na.rm=TRUE))
#GrowthE_CV=cv(GrwthRateE_Relative, na.rm=TRUE))
ERNA.cv$Sum_CV <- rowSums(ERNA.cv[2:5]) #Order by sum of cv
ERNA.cv <- ERNA.cv[order(ERNA.cv$Sum_CV),]
ERNA.cv <- left_join(ERNA.cv, ERNA.SdZn, by="Source")
## ----------------------------


## Make a vertical stripchart showing the CVs for all traits for each population
par(mfrow=c(1,2))

ERNA.cvH <- as.data.frame(cbind(ERNA.cv$Source, ERNA.cv$Height_CV))
ERNA.cvH$Trait <- "Height_CV"
ERNA.cvG <- as.data.frame(cbind(ERNA.cv$Source, ERNA.cv$Growth_CV))
ERNA.cvG$Trait <- "Growth_CV"
#ERNA.cvGE <- as.data.frame(cbind(ERNA.cv$Source, ERNA.cv$GrowthE_CV))
#ERNA.cvGE$Trait <- "GrowthE_CV"
ERNA.cvP <- as.data.frame(cbind(ERNA.cv$Source, ERNA.cv$DaysToFlwr_CV))
ERNA.cvP$Trait <- "DaysToFlwr_CV"
ERNA.cvAll <- rbind(ERNA.cvH, ERNA.cvG, ERNA.cvP) #ERNA.cvGE, 
colnames(ERNA.cvAll) <- c("Source","CV","Trait")

stripchart(as.numeric(CV) ~ Source, data=ERNA.cvAll, vertical=TRUE, pch=19, group.names=ERNA.cv$PopAbbrev, 
           las=2, ylab="Coefficient of variation", cex.lab=1.4, cex=1.25)
ERNA.Gcol <- ERNA.cvAll$Trait == "Growth_CV"
stripchart(as.numeric(CV) ~ Source, data=ERNA.cvAll[ERNA.Gcol,], col="dodgerblue", vertical=TRUE, pch=19,cex=1.25, add=TRUE)
#ERNA.GEcol <- ERNA.cvAll$Trait == "GrowthE_CV"
#stripchart(as.numeric(CV) ~ Source, data=ERNA.cvAll[ERNA.GEcol,], col="blueviolet", vertical=TRUE, pch=19,cex=1.25, add=TRUE)
ERNA.Pcol <- ERNA.cvAll$Trait == "DaysToFlwr_CV"
stripchart(as.numeric(CV) ~ Source, data=ERNA.cvAll[ERNA.Pcol,], col="bisque3", vertical=TRUE, pch=19,cex=1.25, add=TRUE)

## Make stacked barplot to visualize cv for each trait and population
ERNA.cvT <- as.matrix(rbind(as.vector(ERNA.cv$Height_CV), as.vector(ERNA.cv$Growth_CV),
                            as.vector(ERNA.cv$DaysToFlwr_CV)))#, as.vector(ERNA.cv$GrowthE_CV)
colnames(ERNA.cvT) <- ERNA.cv$Source

barplot(ERNA.cvT, names=ERNA.cv$PopAbbrev, las=2, col=c("black","dodgerblue","bisque3"),
        ylab="Coefficient of variation", cex.lab=1.3)
plot.new()
legend("center", c("Flowering phenology","Growth rate","Plant height"), 
       col=c("bisque3","dodgerblue","black"), cex=1.7, pch=15, bty="n")
## ------------------------------------------------------------------------------






## ERNA - EVALUATE RELATIONSHIPS B/W TRAITS AND SOURCE CLIMATE ------------------------------------------
## Look at PCA of 19 bioclim variables to reduce number of predictors 
ERNA.pcaBiovar <- prcomp(ERNA.biovar[,2:20], scale=TRUE)
par(pty="s")
par(mfrow=c(1,1))
plot(x=ERNA.pcaBiovar$x[,1], y=ERNA.pcaBiovar$x[,2], pch=19, cex=1.4, col=ERNA.biovar$HexCode)
plot(x=ERNA.pcaBiovar$x[,2], y=ERNA.pcaBiovar$x[,3], pch=19, cex=1.4, col=ERNA.biovar$HexCode)
plot(x=ERNA.pcaBiovar$x[,3], y=ERNA.pcaBiovar$x[,4], pch=19, cex=1.4, col=ERNA.biovar$HexCode)

## Look at scree plot
ERNA.pcaBiovar$sdev[1]**2/sum(ERNA.pcaBiovar$sdev**2)
ERNA.pcaBV.varExpl <- ERNA.pcaBiovar$sdev^2/sum(ERNA.pcaBiovar$sdev^2)
barplot(ERNA.pcaBiovar$sdev[1:19]**2/sum(ERNA.pcaBiovar$sdev**2))
sum(ERNA.pcaBV.varExpl[1:2]) #top 2 PC axes explain over 63% of variation

## Add arrows on PCA plot and look at loadings 
biplot(ERNA.pcaBiovar)
ERNA.pcaBiovar$rotation #loadings
ERNA.pcaBiovar$rotation[,1:2] #loadings
ERNA.pc1 <- ERNA.pcaBiovar$rotation[,1][order(abs(ERNA.pcaBiovar$rotation[,1]))]
ERNA.pc2 <- ERNA.pcaBiovar$rotation[,2][order(abs(ERNA.pcaBiovar$rotation[,2]))]


## Models
## Use 'top' biovars as determined by PCA
## The following bioclim vars selected due to high loadings in orthogonal directions along PC1 and PC2 
##BIO17 (or 14), BIO11 (or 1), BIO19 (or 16 or 13), BIO7 (or 8 or 4)

bv.names <- c("Temperature Annual Range (BIO7)",
              "Mean Temp. Coldest Quarter (Deg. C)","Precipitation Driest Quarter (BIO17)",
              "Precipitation Coldest Quarter (BIO19)")

ERNA.ht.bv.mod <- lmer(Length_cm_20220915 ~ scale(bio7) + scale(bio11) + scale(bio17) + scale(bio19) 
                       + (1|Block), data=ERNA22.cl)
summary(ERNA.ht.bv.mod)
Anova(ERNA.ht.bv.mod)
plot(allEffects(ERNA.ht.bv.mod))
plot(predictorEffects(ERNA.ht.bv.mod))


#Select model form for growth **


## Is negative binomial an appropriate model form?
hist(ERNA23$DaysToFlwr)
ERNA.df.bv.mod <- glmer.nb(DaysToFlwr ~ scale(bio7) + scale(bio11) + scale(bio17) + scale(bio19)
                           + (1|Block), data=ERNA23)
Anova(ERNA.df.bv.mod)
summary(ERNA.df.bv.mod)
plot(allEffects(ERNA.df.bv.mod))
plot(predictorEffects(ERNA.df.bv.mod))


## Add survival model **
## ---------------------------------------


## Visualize raw data relationships with bioclimate vars
plot(ERNA.cl$Length_cm_20220915 ~ ERNA.cl$bio17) #Is there a curved relationship?
abline(v=ERNA.biovar$bio17[20], col="red")
plot(ERNA.cl$Length_cm_20220915 ~ ERNA.cl$bio19)
abline(v=ERNA.biovar$bio19[20], col="red")
plot(ERNA.cl$Length_cm_20220915 ~ ERNA.cl$bio11)
abline(v=ERNA.biovar$bio11[20], col="red")
plot(ERNA23$DaysToFlwr ~ ERNA23$bio11)
abline(v=ERNA.biovar$bio11[20], col="red")



## Calculate trait means
ERNA.mn <- ERNA22.cl %>% group_by(Source) %>% 
  dplyr::summarise(Height_MN=mean(Length_cm_20220915,na.rm=TRUE), Height_SE=calcSE(Length_cm_20220915),
                   GrowthRe_MN=mean(GrwthRate_Relative,na.rm=TRUE), GrowthRe_SE=calcSE(GrwthRate_Relative),
                   GrowthReE_MN=mean(GrwthRateE_Relative,na.rm=TRUE),GrowthReE_SE=calcSE(GrwthRateE_Relative))
ERNA.mn <- left_join(ERNA.mn, ERNA.biovar, by="Source")



## Plot trait means and sE and add model lines to predictors with model support
## Size
Anova(ERNA.ht.bv.mod)
eff.ht.bio11 <- as.data.frame(predictorEffects(ERNA.ht.bv.mod)[2]) #Extract values for plotting model lines
eff.ht.bio17 <- as.data.frame(predictorEffects(ERNA.ht.bv.mod)[3])
eff.ht.bio19 <- as.data.frame(predictorEffects(ERNA.ht.bv.mod)[4])

plot(as.vector(t(ERNA.mn[,19])), ERNA.mn$Height_MN, bg=ERNA.mn$HexCode, pch=21, col="black", cex=1.5, main=NA, 
     cex.main=1.5, xlab=bv.names[2], ylab="Height (cm)", cex.lab=1.5, cex.axis=1.1, ylim=c(25,60))
arrows(as.vector(t(ERNA.mn[,19])), ERNA.mn$Height_MN-ERNA.mn$Height_SE, as.vector(t(ERNA.mn[,19])), 
       ERNA.mn$Height_MN+ERNA.mn$Height_SE, angle=90, col="grey", code=3, length=0, lwd=1.6) #Error bars
lines(eff.ht.bio11$bio11[,1], eff.ht.bio11$bio11[,2],lwd=1,col="black") #Model fit
abline(v=ERNA.biovar$bio11[20], col="red")

plot(as.vector(t(ERNA.mn[,25])), ERNA.mn$Height_MN, bg=ERNA.mn$HexCode, pch=21, col="black", cex=1.5, main=NA, 
     cex.main=1.5, xlab=bv.names[3], ylab="Height (cm)", cex.lab=1.5, cex.axis=1.1, ylim=c(25,60))
arrows(as.vector(t(ERNA.mn[,25])), ERNA.mn$Height_MN-ERNA.mn$Height_SE, as.vector(t(ERNA.mn[,25])), 
       ERNA.mn$Height_MN+ERNA.mn$Height_SE, angle=90, col="grey", code=3, length=0, lwd=1.6) #Error bars
lines(eff.ht.bio17$bio17[,1], eff.ht.bio17$bio17[,2],lwd=1,col="black") #Model fit
abline(v=ERNA.biovar$bio17[20], col="red")

plot(as.vector(t(ERNA.mn[,27])), ERNA.mn$Height_MN, bg=ERNA.mn$HexCode, pch=21, col="black", cex=1.5, main=NA, 
     cex.main=1.5, xlab=bv.names[4], ylab="Height (cm)", cex.lab=1.5, cex.axis=1.1, ylim=c(25,60))
arrows(as.vector(t(ERNA.mn[,27])), ERNA.mn$Height_MN-ERNA.mn$Height_SE, as.vector(t(ERNA.mn[,27])), 
       ERNA.mn$Height_MN+ERNA.mn$Height_SE, angle=90, col="grey", code=3, length=0, lwd=1.6) #Error bars
lines(eff.ht.bio19$bio19[,1], eff.ht.bio19$bio19[,2],lwd=1,col="black") #Model fit
abline(v=ERNA.biovar$bio19[20], col="red")




## Phenology
Anova(ERNA.df.bv.mod)
eff.df.bio11 <- as.data.frame(predictorEffects(ERNA.df.bv.mod)[2]) #Extract values for plotting model lines

png("20250116_ERNAphenoVSclim.png", width=550, height=550, res=300, pointsize=6)
plot(as.vector(t(ERNA23.mn[,20])), ERNA23.mn$DaysToFlwr_MN, bg=ERNA23.mn$HexCode, pch=21, col="black", cex=1.2, main=NA, 
     cex.main=1.5, xlab=bv.names[2], ylab="Days to first flower", cex.lab=1, cex.axis=1, ylim=c(187,240))#,
#xlab="Mean Temp. Coldest Quarter (Deg. C)")
arrows(as.vector(t(ERNA23.mn[,20])), ERNA23.mn$DaysToFlwr_MN-ERNA23.mn$DaysToFlwr_SE, as.vector(t(ERNA23.mn[,20])), 
       ERNA23.mn$DaysToFlwr_MN+ERNA23.mn$DaysToFlwr_SE, angle=90, col="grey", code=3, length=0, lwd=0.9) #Error bars
lines(eff.df.bio11$bio11[,1], eff.df.bio11$bio11[,2],lwd=0.9,col="black") #Model fit
#abline(v=ERNA.biovar$bio11[20], col="red")
dev.off()

## Blank plot
plot.new()
par(mfrow=c(1,1))
legend("center", unique(ERNA23.mn$SdZnAbbrev[order(ERNA23.mn$SdZnOrder, decreasing=FALSE)]), col="black",
       pt.bg=unique(ERNA23.mn$HexCode[order(ERNA23.mn$SdZnOrder, decreasing=FALSE)]), cex=1.65, pch=21)
## -----------------------------------------------------------------------------------------------------
