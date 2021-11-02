rm(list=ls())
sampleID <- 4
r_no = regions = 1 ### forest center ID (metakeskus) 1:15
#nSetRuns = 10 #number of set runs
#harvestscenarios="Base"		##management scenarios it can be  ### c("Low","MaxSust","NoHarv","Base")
rcpfile="CurrClim"
#regSets<-"maakunta"
source("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/localSettings.r")

UQanalysis <- "True"

#devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/sampleRun.r")
##### From GitHub
devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/virpi-j/UQ_runs/main/functions.R")
#source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
#source("Rsrc/virpiSbatch/functionsUQ.r")
source("Rsrc/virpiSbatch/localSettings.r")

# Give new set of outputs ------------------------------------------------
varOuts <- c("NEP","GPPTot/1000", "npp","V")
varSel <- match(varOuts,varNames)
funX <- rep("sum",length(varSel))
funX[match(varNames[c(7,11:12)],varNames[varSel])] <- "baWmean"
#----------------------------------------------------------------------------

#sampleRun <- FALSE
set.seed(10)
uncRun <- TRUE
#nSitesRun <- 100

#----------------------------------------------------------------------------
if(uncRun){ 
  sampleIDs <- 1:nSamples
  area_total <- sum(data.all$area)
  areas <- data.all$area/area_total
  
  print(paste0("Sample size ",nSitesRun," segments"))
  opsInd <- matrix(0, nSitesRun, nSamples) 
  for(ij in 1:nSamples) opsInd[,ij] <- sample(1:nrow(data.all), nSitesRun, replace = TRUE, prob = areas)
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
#if(sampleRun){
#  toMem <- ls()
#  startRun <- Sys.time() 
#  sampleX <- runModelOrig(sampleID,sampleRun=F, uncRun = uncRun)
#  endRun <- Sys.time()
#  timeRun <- endRun - startRun
#  print(paste0("Run time for sample size ", nSitesRun," = ",timeRun))
#} else {
rcpfile = rcps
# for(rcpfile in rcps) { ## ---------------------------------------------
#  print(rcpfile)
if(rcpfile=="CurrClim"){
  load(paste(climatepath, rcpfile,".rdata", sep=""))  
  #####process data considering only current climate###
  # dat <- dat[rday %in% 1:10958] #uncomment to select some years (10958 needs to be modified)
  maxRday <- max(dat$rday)
  xday <- c(dat$rday,(dat$rday+maxRday),(dat$rday+maxRday*2))
  dat = rbind(dat,dat,dat)
  dat[,rday:=xday]
  
} else{
  load(paste(climatepath, rcpfile,".rdata", sep=""))  
}
##
niter <- ceiling(nSamples/5)

sampleOutput <- list()

for(nii in 1:niter){
  toMem <- ls()
  print("Start running iter ",nii,"/",niter,"...")
  startRun <- Sys.time() 
  #sampleXs <- lapply(sampleIDs[1:2], function(jx) {
  #  runModelUQ(jx,  ## Do nothing for 10 seconds
  #  uncRun = TRUE)})      
  sampleXs <- mclapply(sampleIDs[(1+(nii-1)*5):(nii*5)], function(jx) {
    runModelUQ(jx,  ## Do nothing for 10 seconds
    uncRun = uncRun)}, 
    mc.cores = nCores,mc.silent=FALSE)      ## Split this job across 10 cores
  timeRun <- Sys.time() - startRun
  print(paste0("Run time for ",nSamples," samples of size ", nSitesRun," = ",timeRun))
  print("End running...")

  save(sampleXs,file=paste0("Rsrc/virpiSbatch/results/samplex_",r_no,".rdata")) 
  setwd("Rsrc/virpiSbatch/")

  m <- nrow(sampleXs[[1]])
  n <- length(sampleXs)
  varNams <-  sampleXs[[1]][,"vari"]
  cS <- c(-44/12/10^12, 16^2, 16^2, 0.16^2)
  # g /m2 /year -> -44/12*16^2/10^12
  for(j in 1:m){
    x <- data.frame()
    for(k in 1:n){
      x <- rbind(x, sampleXs[[k]][j,])
    }
    x[,3:5] <- x[,3:5]*cS[j]
    #assign(varNams[j,1], x)
    sampleOutputx[[j]][((nii-1)*5):(nii*5),] <- x
  }
}
  #source("postprocessResults.R")

print("make histograms...")
pdf(paste0("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/figures/results_regionID_",r_no,".pdf"))
m <- nrow(sampleXs[[1]])
n <- length(sampleXs)
varNams <-  sampleXs[[1]][,"vari"]
cS <- c(-44/12/10^12, 16^2, 16^2, 0.16^2)
  # g /m2 /year -> -44/12*16^2/10^12
# g/m2/a -> area_tot *(g/m2/a * area_sample)/area_sample

par(mfrow=c(m,3))
#par(mfrow=c(1,1))
for(j in 1:m){
  x <- data.frame()
  for(k in 1:n){
    x <- rbind(x, sampleXs[[k]][j,])
  }
  x[,3:5] <- x[,3:5]*cS[j]
  #assign(varNams[j,1], x)
  xlims <- c(min(x[,3:5]),max(x[,3:5]))
  for(per in 1:3){
    hist(x[,2+per], main = paste0("period",per), xlab = varNams[j], xlim = xlims)  
  }
}

dev.off()
print("histograms made")