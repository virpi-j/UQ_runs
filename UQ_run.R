rm(list=ls())
sampleID <- 4
#r_no = regions = 1 ### forest center ID (metakeskus) 1:15
#nSetRuns = 10 #number of set runs
#harvestscenarios="Base"		##management scenarios it can be  ### c("Low","MaxSust","NoHarv","Base")
rcpfile="CurrClim"
#regSets<-"maakunta"
source("localSettings.r")

UQanalysis <- "True"

#devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/sampleRun.r")
##### From GitHub
devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/virpi-j/UQ_runs/main/UQ_functions.R")
#           https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
#source("Rsrc/virpiSbatch/functionsUQ.r")

# Give new set of outputs ------------------------------------------------
varOuts <- c("NEP sp","GPPspecies", "npp")
varSel <- match(varOuts,varNames)
funX <- rep("sum",length(varSel))
funX[match(varNames[c(7,11:12)],varNames[varSel])] <- "baWmean"
#----------------------------------------------------------------------------

source("Rsrc/virpiSbatch/localSettings.r")

set.seed(10)
#----------------------------------------------------------------------------

if(UQanalysis){ 
  sampleIDs <- 1:nSamples
  print(paste0("Sample size ",nSitesRun," segments"))
  opsInd <- matrix(0, nSitesRun, nSamples) 
  for(ij in 1:nSamples) opsInd[,ij] <- sample(1:nrow(data.all), nSitesRun)
  save(opsInd,file=paste0("Rsrc/virpiSbatch/results/opsInd_reg",r_no,".rdata")) 
} else {
  setX=1
  nSamples <- ceiling(dim(data.all)[1]/nSitesRun)
  sampleIDs <- split(1:nSamples,             # Applying split() function
                     cut(seq_along(1:nSamples),
                         nSetRuns,
                         labels = FALSE))[[setX]]
  ops <- split(data.all, sample(1:nSamples, nrow(data.all), replace=T))
}

#library(ff)

toMem <- ls()
startRun <- Sys.time() 
sampleX <- runModelOrig(sampleID,sampleRun=F)
endRun <- Sys.time()
timeRun <- endRun - startRun
print(paste0("Run time for sample size ", nSitesRun," = ",timeRun))
#mclapply(sampleIDs, function(jx) {
#  runModelOrig(jx)  ## Do nothing for 10 seconds
#}, mc.cores = nCores,mc.silent=FALSE)      ## Split this job across 10 cores

setwd("Rsrc/virpiSbatch/")

#source("postprocessResults.R")