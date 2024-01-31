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
library(dismo)
## ------------------------------------------------------------------------------------------------




## DOWNLOAD AND EXTRACT SEASONAL VALUES FROM PAST (E.G. 30 YEARS) --------------------------------

## FROM PRISM PACKAGE TUTORIAL 
## https://github.com/ropensci/prism

## Set directory where data is downloaded
prism_set_dl_dir("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate")

## Download data
#get_prism_monthlys(type=c("ppt"), year=2001:2020, mon=1:12, keepZip=FALSE)
## ----------------------------------------------------------------------------------------------



## LOAD LAT/ LONGS AND ASSIGN DESIRED TEMPORAL RANGE FOR CLIMATE DATA ---------------------------
pops_LatLong <- read.table(file ='C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds/20220622_ERNA_LatLong.csv', sep=',', header = TRUE)  
#pops_LatLong <- read.table(file ='C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds/20221129_BOGR_LatLong.csv', sep=',', header = TRUE)  
#pops_LatLong <- read.table(file ='C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds/20220929_ARFR_LatLong.csv', sep=',', header = TRUE)  
num.pops <- nrow(pops_LatLong)
colnames(pops_LatLong)
colnames(pops_LatLong)[1] <- "Code" 



## MAKE LIST OF DESIRED YEAR & MONTH TO EXTRACT CLIMATE DATA FROM
year.span <- as.character(c(1980:2021))
month.span <- str_pad(as.character(c(1:12)), width=2, side="left", pad="0")

yyyymm <- NULL 

for (yy in 1:length(year.span)) {
  for (mm in 1:length(month.span)) {
    yyyymm <- c(yyyymm,paste(year.span[yy], month.span[mm], sep=""))
  }
}

## Morph year & month labels into matrix  
yyyymm <- as.data.frame(matrix(yyyymm, length(month.span), length(year.span)), row.names = FALSE)
colnames(yyyymm) <- as.character(c(1980:2021))
## --------------------------------------------------------------------------------------------





## MONTHLY PRECIP ----------------------------------------------------------
setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.ppt <- prism_archive_subset("ppt", "monthly", year=1980:2021)
bils.ppt <- pd_to_file(dirs.ppt)


## OPTIONAL **
## Truncate lat-longs to get data from missing pop only
#pops_LatLong <- pops_LatLong[9:13,]
#num.pops <- nrow(pops_LatLong)
## ---------------------------------------------------



## Extract precip data  
ppt <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(ppt) <- c("Year", 1:12)
ppt$Year <- colnames(yyyymm)

ppt.means <- as.data.frame(matrix(NA, num.pops, (nrow(yyyymm)+1)))
colnames(ppt.means) <- c("Pop", 1:12)
ppt.means$Pop <- pops_LatLong$Code

#ppt.sd <- as.data.frame(matrix(NA, num.pops, (nrow(yyyymm)+1)))
#colnames(ppt.sd) <- c("Pop", 1:12)
#ppt.sd$Pop <- pops_LatLong$Code

for (pp in 1:num.pops) {
  
  for (yy in 1:nrow(ppt)) {
    
    for (mm in 1:12) {
    file.pos <- grep(yyyymm[mm,yy], bils.ppt)
    raster_file <- raster(bils.ppt[file.pos])                                   #Load specified raster
    ppt[yy,mm+1] <- raster::extract(x = raster_file, pops_LatLong[pp,2:3])      #Extract ppt for given gps coord & month
    
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
species <- as.character("ARFR")
saveRDS(ppt.means, file=paste(date,species,"pptMonthly",sep="_"))

## Sum ppt over all months (annual precip)
ppt.means.ann <- rowSums(ppt.means[,2:13])
ppt.means.ann <- as.data.frame(cbind(ppt.means$Pop,as.numeric(ppt.means.ann)))
colnames(ppt.means.ann) <- c("Source","Mean_Annual_Ppt")

## Save output
write.csv(ppt.means.ann, file=paste(date,species,"pptAnnual.csv",sep="_"), row.names=FALSE)
## ------------------------------------------------------------------------------





## MONTHLY MIN TEMP ----------------------------------------------------------
setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.tmin <- prism_archive_subset("tmin", "monthly", year=1980:2021)
bils.tmin <- pd_to_file(dirs.tmin)


## OPTIONAL **
## Truncate lat-longs to get data from missing pop only
#pops_LatLong <- pops_LatLong[9:13,]
#num.pops <- nrow(pops_LatLong)
## ---------------------------------------------------


## Extract minimum temp data  
tmin <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(tmin) <- c("Year", 1:12)
tmin$Year <- colnames(yyyymm)

tmin.means <- as.data.frame(matrix(NA, num.pops, (nrow(yyyymm)+1)))
colnames(tmin.means) <- c("Pop", 1:12)
tmin.means$Pop <- pops_LatLong$Code

tmin.sd <- as.data.frame(matrix(NA, num.pops, (nrow(yyyymm)+1)))
colnames(tmin.sd) <- c("Pop", 1:12)
tmin.sd$Pop <- pops_LatLong$Code

for (pp in 1:num.pops) {
  
  for (yy in 1:nrow(tmin)) {
    
    for (mm in 1:12) {
      file.pos <- grep(yyyymm[mm,yy], bils.tmin)
      raster_file <- raster(bils.tmin[file.pos])                                   #Load specified raster
      tmin[yy,mm+1] <- raster::extract(x = raster_file, pops_LatLong[pp,2:3])      #Extract tmin for given gps coord & month
      
      ## Calculate monthly means (averaged over all years)
      tmin.means[pp,mm+1] <- mean(tmin[,mm+1])
      ## Calculate variation (sd) per month (for each month, sd over all years) 
      tmin.sd[pp,mm+1] <- sd(tmin[,mm+1])
      
    }
  }
}


setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds")
#date <- as.character(20240131)
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ERNA")
saveRDS(tmin.means, file=paste(date,species,"tminMonthly",sep="_"))
#saveRDS(tmin.sd, file=paste(date,species,"tminSDmonthly",sep="_"))

#tmin.means <- readRDS("20221129_BOGR_tminMonthly")
#tmin.means <- readRDS("20221129_ERNA_tminMonthly")
tmin.means

## Average min temp from Dec-Feb (min winter temp)
winter.months <- c(2:3, 13)
tmin.means.winter <- rowMeans(tmin.means[,winter.months])
tmin.means.winter <- as.data.frame(cbind(tmin.means$Pop,as.numeric(tmin.means.winter)))
colnames(tmin.means.winter) <- c("Code","Mean_MinWinter_Temp")

## Save output
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ARFR")
write.csv(tmin.means.winter, file=paste(date,species,"tminWinter.csv",sep="_"), row.names=FALSE)
## ------------------------------------------------------------------------------





## MONTHLY MAX TEMP ----------------------------------------------------------
#setwd("Q:/Research/All_Projects_by_Species/Astragalus SPECIES/Astragalus_microcymbus/PRISM_asmiclimate") #Directory containing climate data 
dirs.tmax <- prism_archive_subset("tmax", "monthly", year=1980:2021)
bils.tmax <- pd_to_file(dirs.tmax)


## Extract maximum temp data  
tmax <- as.data.frame(matrix(NA, ncol(yyyymm), (nrow(yyyymm)+1)))
colnames(tmax) <- c("Year", 1:12)
tmax$Year <- colnames(yyyymm)

## Make variable for storing mean max temp (averaged over all years)
tmax.means <- as.data.frame(matrix(NA, num.pops, (nrow(yyyymm)+1)))
colnames(tmax.means) <- c("Pop", 1:12)
tmax.means$Pop <- pops_LatLong$Code

for (pp in 1:num.pops) {
  
  for (yy in 1:nrow(tmax)) {
    
    for (mm in 1:12) {
      file.pos <- grep(yyyymm[mm,yy], bils.tmax)
      raster_file <- raster(bils.tmax[file.pos])                                   #Load specified raster
      tmax[yy,mm+1] <- raster::extract(x = raster_file, pops_LatLong[pp,2:3])      #Extract tmax for given gps coord & month
      
      ## Calculate monthly means (averaged over all years)
      tmax.means[pp,mm+1] <- mean(tmax[,mm+1])
      
    }
  }
}


setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds")
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ERNA")
saveRDS(tmax.means, file=paste(date,species,"tmaxMonthly",sep="_"))

#tmax.means <- readRDS("20230312_ARFR_tmaxMonthly")
tmax.means
## ------------------------------------------------------------------------------





## USE BIOVARS FUNCTION TO ESTIMATE 19 BIOVARIABLES
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds")
#ppt.means <- readRDS("20230311_ARFR_pptMonthly")
#tmin.means <- readRDS("20230311_ARFR_tminMonthly")
#tmax.means <- readRDS("20230314_ARFR_tmaxMonthly")
ppt.means <- readRDS("20230219_ERNA_pptMonthly")
tmin.means <- readRDS("20221129_ERNA_tminMonthly")
tmax.means <- readRDS("20230825_ERNA_tmaxMonthly")
#ppt.means <- readRDS("20221129_BOGR_pptMonthly")
#tmin.means <- readRDS("20221129_BOGR_tminMonthly")
#tmax.means <- readRDS("20230824_BOGR_tmaxMonthly")

ppt.means.nopops <- ppt.means[,2:13]
tmin.means.nopops <- tmin.means[,2:13]
tmax.means.nopops <- tmax.means[,2:13]
ppt.means.t <- t(ppt.means.nopops)
tmin.means.t <- t(tmin.means.nopops)
tmax.means.t <- t(tmax.means.nopops)

biovar.means <- as.data.frame(matrix(NA, num.pops, 20))
colnames(biovar.means) <- c("Pop", as.character(1:19))
biovar.means$Pop <- pops_LatLong$Code

for (pp in 1:num.pops) {
  bioVars.BOGR <- biovars(as.vector(ppt.means.t[,pp]), as.vector(tmin.means.t[,pp]),
                          as.vector(tmax.means.t[,pp]))
  
  biovar.means[pp,2:20] <- bioVars.BOGR
}


#setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/AGoebl/Seeds")
date <- Sys.Date()  
date <- gsub("-", "", date)
species <- as.character("ERNA")
saveRDS(biovar.means, file=paste(date,species,"BiovarsAvg1980_2021",sep="_"))
#write.csv(biovar.means, file=paste(date,species,"BiovarsAvg1980_2021.csv",sep="_"), row.names=FALSE)

