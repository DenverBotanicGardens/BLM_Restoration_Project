## April Goebl
## Script started 2022-03-31
## Analyze Chatfield greenhouse emergence data 
## BLM restoration species



rm(list=ls())
setwd("C:/Users/april.goebl/Denver Botanic Gardens/Conservation - Restoration/BLM-Grassland/Greenhouse")


## LOAD PACKAGES ---------------------------------------
library(lubridate)
## -----------------------------------------------------



## LOAD DATA & CHANGE STRUCTURE ------------------------
PASM <- read.csv(file="20220331_GreenhouseDataComb_PASM.csv", sep=",", header=TRUE, dec=".")
HEVI <- read.csv(file="20220401_GreenhouseDataComb_HEVI.csv", sep=",", header=TRUE, dec=".")
PEVI <- read.csv(file="20220401_GreenhouseDataComb_PEVI.csv", sep=",", header=TRUE, dec=".")
ARFR <- read.csv(file="20220402_GreenhouseDataComb_ARFR.csv", sep=",", header=TRUE, dec=".")
ERNA <- read.csv(file="20220402_GreenhouseDataComb_ERNA.csv", sep=",", header=TRUE, dec=".")
BOGR <- read.csv(file="20220404_GreenhouseDataComb_BOGR.csv", sep=",", header=TRUE, dec=".")


str(PASM)
PASM$DateFirstRecorded <- as.Date(PASM$DateFirstRecorded, "%m/%d/%Y")
PASM$doy <- yday(PASM$DateFirstRecorded) #change date to day of year

str(HEVI)
HEVI$DateFirstRecorded <- as.Date(HEVI$DateFirstRecorded, "%m/%d/%Y")
HEVI$doy <- yday(HEVI$DateFirstRecorded) #change date to day of year

str(PEVI)
PEVI$DateFirstRecorded <- as.Date(PEVI$DateFirstRecorded, "%m/%d/%Y")
PEVI$doy <- yday(PEVI$DateFirstRecorded) #change date to day of year

str(ARFR)
ARFR$DateFirstRecorded <- as.Date(ARFR$DateFirstRecorded, "%m/%d/%Y")
ARFR$doy <- yday(ARFR$DateFirstRecorded) #change date to day of year

str(ERNA)
ERNA$DateFirstRecorded <- as.Date(ERNA$DateFirstRecorded, "%m/%d/%Y")
ERNA$doy <- yday(ERNA$DateFirstRecorded) #change date to day of year

str(BOGR)
BOGR$DateFirstRecorded <- as.Date(BOGR$DateFirstRecorded, "%m/%d/%Y")
BOGR$doy <- yday(BOGR$DateFirstRecorded) #change date to day of year
## -----------------------------------------------------



## Calculate & plot num of days to emerge 
par(mfrow=c(3,2))


## ERNA
dayOutOfStrat_ERNA <- yday(as.Date("1/3/2022", "%m/%d/%Y"))
ERNA$daysToEmg <- ERNA$doy - dayOutOfStrat_ERNA 

pops_ERNA <- unique(ERNA$Code)

boxplot(ERNA$daysToEmg ~ ERNA$Code, xlab=NA, ylab="Days to emerge", cex.lab=1.5,
        main="Ericameria nauseosa", cex.main=1.5, col="palegoldenrod", las=2, cex.axis=0.8)

      



## BOGR
dayOutOfStrat_BOGR <- yday(as.Date("1/14/2022", "%m/%d/%Y"))
BOGR$daysToEmg <- BOGR$doy - dayOutOfStrat_BOGR 

pops_BOGR <- unique(BOGR$Code)


boxplot(BOGR$daysToEmg ~ BOGR$Code,
        xlab=NA, ylab="Days to emerge", cex.lab=1.5, las=2,
        cex.axis=0.8,
        main="Bouteloua gracilis", cex.main=1.5, col="lightcyan")
#lightskyblue1




## PEVI
dayOutOfStrat_PEVI <- yday(as.Date("1/17/2022", "%m/%d/%Y"))
PEVI$daysToEmg <- PEVI$doy - dayOutOfStrat_PEVI 

pops_PEVI <- unique(PEVI$Code)

boxplot(PEVI$daysToEmg ~ PEVI$Code,
        xlab=NA, ylab="Days to emerge", cex.lab=1.5,
        main="Penstemon virens", cex.main=1.5, col="thistle")
#plum1


## ARFR
dayOutOfStrat_ARFR <- yday(as.Date("1/3/2022", "%m/%d/%Y"))
ARFR$daysToEmg <- ARFR$doy - dayOutOfStrat_ARFR 

pops_ARFR <- unique(ARFR$Code)

boxplot(ARFR$daysToEmg ~ ARFR$Code, las=2,
        xlab=NA, ylab="Days to emerge", cex.lab=1.5,
        cex.axis=0.9, main="Artemisia frigida", cex.main=1.5, col="darkseagreen3")


## HEVI
dayOutOfStrat_HEVI <- yday(as.Date("1/3/2022", "%m/%d/%Y"))
HEVI$daysToEmg <- HEVI$doy - dayOutOfStrat_HEVI 

pops_HEVI <- unique(HEVI$Code)

boxplot(HEVI$daysToEmg ~ HEVI$Code,
        xlab="Population", ylab="Days to emerge", cex.lab=1.5,
        main="Heterotheca villosa", cex.main=1.5, col="lightyellow")
#lightyellow


## PASM
dayOutOfStrat_PASM <- yday(as.Date("1/3/2022", "%m/%d/%Y"))
PASM$daysToEmg <- PASM$doy - dayOutOfStrat_PASM 

pops_PASM <- unique(PASM$Code)

boxplot(PASM$daysToEmg ~ PASM$Code,
        xlab="Population", ylab="Days to emerge", cex.lab=1.5,
        main="Pascopyrum smithii", cex.main=1.5, col="palegreen1")

