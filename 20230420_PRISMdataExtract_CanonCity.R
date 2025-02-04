## April Goebl
## Script modified 06-20-2022
## BLM Restoration project at Denver Botanic Gardens
## Extract PRISM climate data for experimental populations 





rm(list=ls())



## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(raster)
library(rgdal)
library(dplyr)
library(stringr)
library(tidyr)
library(prism)
## ------------------------------------------------------------------------------------------------




## DOWNLOAD AND EXTRACT SEASONAL VALUES FROM PAST (E.G. 30 YEARS) --------------------------------

## FROM PRISM PACKAGE TUTORIAL 
## https://github.com/ropensci/prism

## Set download directory
prism_set_dl_dir("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate")

## Download data
#get_prism_monthlys(type=c("ppt"), year=2001:2020, mon=1:12, keepZip=FALSE)
## ----------------------------------------------------------------------------------------------



## LOAD LAT/ LONGS AND ASSIGN DESIRED TEMPORAL RANGE FOR CLIMATE DATA ---------------------------
#sites_LatLong <- read.table(file ='20221129_BOGR_LatLong.csv', sep=',', header = TRUE)  
sites <- c("PC", "EP")
longs <- as.numeric(c(-105.030438, -105.26563))
lats <- as.numeric(c(38.503163, 38.42146))
sites_LatLong <- data.frame(sites,longs,lats)
num.sites <- nrow(sites_LatLong)
#colnames(sites_LatLong)
#colnames(sites_LatLong)[1] <- "sites" 



## MAKE LIST OF DESIRED YEAR & MONTH TO EXTRACT CLIMATE DATA FROM
year.span <- as.character(c(2000:2021))
month.span <- str_pad(as.character(c(1:12)), width=2, side="left", pad="0")

yyyymm <- NULL 

for (yy in 1:length(year.span)) {
  for (mm in 1:length(month.span)) {
    yyyymm <- c(yyyymm,paste(year.span[yy], month.span[mm], sep=""))
  }
}

## Morph year & month labels into matrix  
yyyymm <- as.data.frame(matrix(yyyymm, length(month.span), length(year.span)), row.names = FALSE)
colnames(yyyymm) <- as.character(c(2000:2021))
## --------------------------------------------------------------------------------------------





## MONTHLY PRECIP ----------------------------------------------------------
setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.ppt <- prism_archive_subset("ppt", "monthly", year=2000:2021)
bils.ppt <- pd_to_file(dirs.ppt)
## ---------------------------------------------------



## Extract precip data  
ppt <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(ppt) <- c("Year", 1:12)
ppt$Year <- colnames(yyyymm)

ppt.means <- as.data.frame(matrix(NA, num.sites, (nrow(yyyymm)+1)))
colnames(ppt.means) <- c("Site", 1:12)
ppt.means$Site <- sites_LatLong$sites

#ppt.sd <- as.data.frame(matrix(NA, num.sites, (nrow(yyyymm)+1)))
#colnames(ppt.sd) <- c("Site", 1:12)
#ppt.sd$Site <- sites_LatLong$sites

for (pp in 1:num.sites) {
  
  for (yy in 1:nrow(ppt)) {
    
    for (mm in 1:12) {
    file.pos <- grep(yyyymm[mm,yy], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                   #Load specified raster
    ppt[yy,mm+1] <- raster::extract(x = raster_file, sites_LatLong[pp,2:3])     #Extract ppt for given gps coord & month
    
    ## Calculate monthly means (averaged over all years)
    ppt.means[pp,mm+1] <- mean(ppt[,mm+1])
    ## Calculate variation (sd) per month (for each month, sd over all years) 
    #ppt.sd[pp,mm+1] <- sd(ppt[,mm+1])
    }
  }
}

setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - BLM-Grassland/AGoebl/Seeds")
date <- Sys.Date()  
date <- gsub("-", "", date)
#species <- as.character("ARFR")
#saveRDS(ppt.means, file=paste(date,species,"pptMonthly",sep="_"))

## Sum ppt over all months (annual precip)
#ppt.means.ann <- rowSums(ppt.means[,2:13])
#ppt.means.ann <- as.data.frame(cbind(ppt.means$Site,as.numeric(ppt.means.ann)))
#colnames(ppt.means.ann) <- c("Site","Mean_Annual_Ppt")

## Save output
write.csv(ppt.means, file=paste(date,species,"pptAnnual.csv",sep="_"), row.names=FALSE)
## ------------------------------------------------------------------------------





## MONTHLY MIN TEMP ----------------------------------------------------------
setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.tmin <- prism_archive_subset("tmin", "monthly", year=1980:2021)
bils.tmin <- pd_to_file(dirs.tmin)


## OPTIONAL **
## Truncate lat-longs to get data from missing pop only
#sites_LatLong <- sites_LatLong[9:13,]
#num.sites <- nrow(sites_LatLong)
## ---------------------------------------------------


## Extract minimum temp data  
tmin <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(tmin) <- c("Year", 1:12)
tmin$Year <- colnames(yyyymm)

tmin.means <- as.data.frame(matrix(NA, num.sites, (nrow(yyyymm)+1)))
colnames(tmin.means) <- c("Site", 1:12)
tmin.means$Site <- sites_LatLong$sites

tmin.sd <- as.data.frame(matrix(NA, num.sites, (nrow(yyyymm)+1)))
colnames(tmin.sd) <- c("Site", 1:12)
tmin.sd$Site <- sites_LatLong$sites

for (pp in 1:num.sites) {
  
  for (yy in 1:nrow(tmin)) {
    
    for (mm in 1:12) {
      file.pos <- grep(yyyymm[mm,yy], bils.tmin)
      raster_file <- raster(bils.tmin[file.pos])                                   #Load specified raster
      tmin[yy,mm+1] <- raster::extract(x = raster_file, sites_LatLong[pp,2:3])      #Extract tmin for given gps coord & month
      
      ## Calculate monthly means (averaged over all years)
      tmin.means[pp,mm+1] <- mean(tmin[,mm+1])
      ## Calculate variation (sd) per month (for each month, sd over all years) 
      tmin.sd[pp,mm+1] <- sd(tmin[,mm+1])
      
    }
  }
}


setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - BLM-Grassland/AGoebl/Seeds")
date <- as.character(20230311)
species <- as.character("ARFR")
saveRDS(tmin.means, file=paste(date,species,"tminMonthly",sep="_"))
saveRDS(tmin.sd, file=paste(date,species,"tminSDmonthly",sep="_"))

#tmin.means <- readRDS("20221129_BOGR_tminMonthly")
#tmin.means <- readRDS("20221129_ERNA_tminMonthly")
tmin.means

## Average min temp from Dec-Feb (min winter temp)
winter.months <- c(2:3, 13)
tmin.means.winter <- rowMeans(tmin.means[,winter.months])
tmin.means.winter <- as.data.frame(cbind(tmin.means$Site,as.numeric(tmin.means.winter)))
colnames(tmin.means.winter) <- c("sites","Mean_MinWinter_Temp")

## Save output
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ARFR")
write.csv(tmin.means.winter, file=paste(date,species,"tminWinter.csv",sep="_"), row.names=FALSE)
## ------------------------------------------------------------------------------





## MONTHLY MAX TEMP ----------------------------------------------------------
setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.tmax <- prism_archive_subset("tmax", "monthly", year=1980:2021)
bils.tmax <- pd_to_file(dirs.tmax)



## Extract maximum temp data  
tmax <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(tmax) <- c("Year", 1:12)
tmax$Year <- colnames(yyyymm)

## Make variable for storing mean max temp (averaged over all years)
tmax.means <- as.data.frame(matrix(NA, num.sites, (nrow(yyyymm)+1)))
colnames(tmax.means) <- c("Site", 1:12)
tmax.means$Site <- sites_LatLong$sites

for (pp in 1:num.sites) {
  
  for (yy in 1:nrow(tmax)) {
    
    for (mm in 1:12) {
      file.pos <- grep(yyyymm[mm,yy], bils.tmax)
      raster_file <- raster(bils.tmax[file.pos])                                   #Load specified raster
      tmax[yy,mm+1] <- raster::extract(x = raster_file, sites_LatLong[pp,2:3])      #Extract tmax for given gps coord & month
      
      ## Calculate monthly means (averaged over all years)
      tmax.means[pp,mm+1] <- mean(tmax[,mm+1])
      
    }
  }
}


setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - BLM-Grassland/AGoebl/Seeds")
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ARFR")
saveRDS(tmax.means, file=paste(date,species,"tmaxMonthly",sep="_"))

#tmax.means <- readRDS("20230312_ARFR_tmaxMonthly")
tmax.means
## ------------------------------------------------------------------------------










