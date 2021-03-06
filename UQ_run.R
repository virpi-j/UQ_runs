rm(list=ls())
sampleID <- 4
#r_no = regions = 1 ### forest center ID (metakeskus) 1:15
#nSetRuns = 10 #number of set runs
#harvestscenarios="Base"		##management scenarios it can be  ### c("Low","MaxSust","NoHarv","Base")
rcpfile="CurrClim"
ststDeadW<-FALSE
#regSets<-"maakunta"
source("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/localSettings.r")
#nSitesRun <- nSitesRunr
#nSamples <- nSamplesr

UQanalysis <- "True"

#devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/sampleRun.r")
##### From GitHub

devtools::source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
source_url("https://raw.githubusercontent.com/virpi-j/UQ_runs/main/functions.R")
#source_url("https://raw.githubusercontent.com/virpi-j/UQ_runs/main/functions.R")


for(r_no in r_nos){
  regions_no <- r_no
  source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/finRuns/Rsrc/settings.r")
  #source("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/localSettings.r")
  #r_no = regions = regions_no 
  nSitesRun <- nSitesRunr
  #nSamples <- nSamplesr
  print(paste("start region",r_no,"- set size",nSitesRun,"- no of samples", nSamplesr))
  
  # Give new set of outputs ------------------------------------------------
  varOuts <- c("NEP","V","npp","VroundWood")
#  cS <- c(-16^2*44/(12*(10^12)), 0.16^2,16^2, 16^2, 0.16^2) # multipliers for tot.sums
  cS <- c(-100^2*44/12, 1, 100^2, 1) # multipliers for tot.sums
  
  varSel <- match(varOuts,varNames)
  funX <- rep("sum",length(varSel))
  funX[match(varNames[c(7,11:12)],varNames[varSel])] <- "baWmean"
  #----------------------------------------------------------------------------

  #sampleRun <- FALSE
  set.seed(10)
  uncRun <- TRUE
  #uncRun <- FALSE
  #nSitesRun <- 100

  #----------------------------------------------------------------------------
  if(uncRun){ 
    sampleIDs <- 1:nSamplesr
    area_total <- sum(data.all$area)
    areas <- data.all$area
    areas <- areas/area_total
    #hist(areas)
    print(paste0("Sample size ",nSitesRunr," segments"))
    opsInd <- list(); #matrix(0, nSitesRun, nSamples) 
    for(ij in 1:nSamplesr){ 
      #opsInd[,ij] <- sample(1:nrow(data.all), nSitesRun, replace = FALSE, prob = areas)
      #opsInd[,ij] <- sample(1:nrow(data.all), nSitesRun, replace = TRUE, prob = areas)
      opsInd[[ij]] <- sample(1:nrow(data.all), nSitesRunr, replace = TRUE, prob = areas)
    }
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
  }else{
    load(paste(climatepath, rcpfile,".rdata", sep=""))  
  }
  ##
  niter <- ceiling(nSamplesr/nParRuns)

  sampleOutput <- list()

  for(nii in 1:niter){
    toMem <- ls()
    print(paste0("Start running iter ",nii,"/",niter,"..."))
    startRun <- Sys.time() 
    #sampleXs <- lapply(sampleIDs[1:4], function(jx) {
    #  runModelUQ(jx,  ## Do nothing for 10 seconds
    #  uncRun = TRUE, ststDeadW=FALSE)})      
    sampleXs <- mclapply(sampleIDs[(1+(nii-1)*nParRuns):(nii*nParRuns)], function(jx) {
      runModelUQ(jx,  ## Do nothing for 10 seconds
      uncRun = uncRun)}, 
      mc.cores = nCores,mc.silent=FALSE)      ## Split this job across 10 cores
    timeRun <- Sys.time() - startRun
    print(paste0("Run time for ",nParRuns," samples of size ", nSitesRunr," = ",timeRun))
    print("End running...")

    save(sampleXs,file=paste0("Rsrc/virpiSbatch/results/samplex_",r_no,".rdata")) 
  
    m <- nrow(sampleXs[[1]])
    n <- length(sampleXs)
    varNams <-  sampleXs[[1]][,"vari"]
    # g /m2 /year -> -44/12*16^2/10^12
    for(j in 1:m){
      x <- data.frame()
      for(k in 1:n){
        x <- rbind(x, sampleXs[[k]][j,])
      }
      x[,3:5] <- x[,3:5]*cS[j]
      #assign(varNams[j,1], x)
      if(nii == 1){
        sampleOutput[[j]] <- x
      } else {
        sampleOutput[[j]] <- rbind(sampleOutput[[j]], x)
        if(j==1){sampleOutput[[j]][,"vari"] <- "NEE"}
      }
    }

    save(sampleOutput,file=paste0("Rsrc/virpiSbatch/results/samplexout",r_no,"samplesize",nSitesRunr,".rdata")) 
  
  }
  #source("postprocessResults.R")

  print("make histograms...")
  m <- length(sampleOutput)
  n <- nrow(sampleOutput[[1]])
  units_hist <- c(10^-12,10^-6,1,10^-6)
  
  #par(mfrow=c(m,3))
  for(j in 1:m){
    x <- sampleOutput[[j]]
    x[,3:5] <- x[,3:5]*units_hist[j]
    varNams <- x[1,"vari"]
    png(file = paste0("/scratch/project_2000994/PREBASruns/finRuns/Rsrc/virpiSbatch/figures/results_regionID",r_no,"_",nSitesRunr,"_",varNams,".png"))
    xlims <- c(min(x[,3:5]),max(x[,3:5]))
    par(mfrow=c(3,1))
    for(per in 1:3){
      hist(x[,2+per], main = paste0("period",per), xlab = varNams, xlim = xlims)  
    }
    dev.off()
    #print(colMeans(x[,3:5]))
  }
  print("histograms made")
  setwd("Rsrc/virpiSbatch/")

}