#rm(list=ls())
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
source_url("https://raw.githubusercontent.com/ForModLabUHel/IBCcarbon_runs/master/general/functions.r")
source("Rsrc/virpiSbatch/functionsUQ.r")

# Give new set of outputs ------------------------------------------------
varOuts <- c("NEP sp","GPPspecies", "npp")
varSel <- match(varOuts,varNames)
funX <- rep("sum",length(varSel))
funX[match(varNames[c(7,11:12)],varNames[varSel])] <- "baWmean"
#----------------------------------------------------------------------------

source("Rsrc/virpiSbatch/localSettings.r")

set.seed(10)
#----------------------------------------------------------------------------
sample_inputs <- 1

if(sample_inputs==1){
  area_quantile_limits <- quantile(data.all$area, probs = c(0,1))#0.25, 0.5, 0.75,1))
  pdf(paste0("Rsrc/virpiSbatch/figures/BootstrapInitVals_regionID_",r_no,".pdf"))
  for(kk in 1:(length(area_quantile_limits)-1)){
    nset <- which(data.all$area >= area_quantile_limits[kk] & data.all$area <= area_quantile_limits[kk+1])
    vals_mat <- data.all[nset,4:11]
    n <- nrow(vals_mat)
    m <- ncol(vals_mat)
    niter <- 1000
    Ymean <- colSums(vals_mat)/n 
    Ynames <- colnames(vals_mat)
    vals_mat <- as.matrix(vals_mat)
    vals_mat[which(vals_mat[,1]>500),1] <- 500
    
    nsamples <- c(1,2,3,4,5,10,20,30,40,50,75,100) # number of segments included
    rmses <- matrix(0,length(nsamples),m)
    areas <- matrix(0,length(nsamples),1)
    
    for(nsamplei in 1:length(nsamples)){
      # nsample <- 10#10^6
      nsample <- nsamples[nsamplei]
      Ysample <- matrix(0,niter,m)
      Asample <- matrix(0,niter,1)
      
      for(ii in 1:niter){
        nn <- sample(1:n, nsample, replace = TRUE)
        if(nsample > 1) {
          Ysample[ii,] <- colSums(vals_mat[nn,]*data.all$area[nset[nn]])/sum(data.all$area[nset[nn]])
          Asample[ii] <- sum(data.all$area[nset[nn]])
        } else {
          Ysample[ii,] <- vals_mat[nn,]
          Asample[ii] <- data.all$area[nset[nn]]
        } 
      }
      
      for(ij in 1:m){
        #hist(Y[,ij])
        if(ij<m){
          for(jj in (ij+1):m) plot(Ysample[,ij],Ysample[,jj], xlab = Ynames[ij], ylab = Ynames[jj], main = paste0("sample size ",nsample))
        }
        rmses[nsamplei,ij] <- sd(Ysample[,ij])/Ymean[ij]
      }
      areas[nsamplei] <- mean(Asample)
      rmses
      areas
    }
    for(ij in 1:m) plot(nsamples, rmses[,ij]*100, main = Ynames[ij], ylab = "Rmse %")
    for(ij in 1:m) plot(areas, rmses[,ij]*100, main = Ynames[ij], ylab = "Rmse %")
    dev.off()
  }
}
#----------------------------------------------------------------------------
figure_of_samplesizes <- 1 # draw histograms of the initial value distributions with different sample sizes
if(figure_of_samplesizes==1){
  nSamplesf <- nSamples
  nSamples <- 1000 # How many times the sample statistics is calculated = simulations made
  sampleIDs <- 1:nSamples
  nSitesRunf <- nSitesRun
  nSitesRun <- 40000 # How many segments are simulated in one sample
  
  pdf(paste0("Rsrc/virpiSbatch/figures/sampleInitVals_regionID_",r_no,".pdf"))
  sampleYnames <- colnames(data.all)
  
  opsInd <- matrix(0, nSitesRun, nSamples) 
  for(ij in 1:nSamples) opsInd[,ij] <- sample(1:nrow(data.all), nSitesRun, replace = TRUE)
  
  areaTot <- sum(data.all$area)
  volTot <- colSums(data.all[,8:11] * data.all$area)#/areaTot
  tot <- sum(volTot)
  volTot["total"] <- tot
  
  vols <- data.frame(matrix(0,nSamples,1))
  nsizes <- c(100,10000, 20000, 40000)
  for(ik in 1:length(nsizes)){ 
    areaSample <- matrix(0,nSamples,1)
    volSample <- matrix(0,nSamples,5)
    for(ij in 1:nSamples){
      inds <- opsInd[1:nsizes[ik],ij]
      areas <- data.all$area[inds]
      areaSample[ij] <- sum(areas)
      volSample[ij,1] <- sum(data.all$pine[inds] * areas)/sum(areas)*areaTot
      volSample[ij,2] <- sum(data.all$spruce[inds] * areas)/sum(areas)*areaTot
      volSample[ij,3] <- sum(data.all$birch[inds] * areas)/sum(areas)*areaTot
      volSample[ij,4] <- sum(data.all$decid[inds] * areas)/sum(areas)*areaTot    
      volSample[ij,5] <- sum(volSample[ij,1:4])
    }
    vols[[ik]] <- data.frame(volSample)
    colnames(vols[[ik]]) <- c("pine","spruce","birch","decid","total")
  }
  par(mfrow=c(length(nsizes),1))
  
  for(ik in 1:ncol(vols[[1]])){
    xmin <- min(vols[[1]][,ik])
    xmax <- max(vols[[1]][,ik])
    for(zk in 1:length(nsizes)){
      #      hist(vols[[zk]][,ik], main = paste0("Whole area vol [m^3/ha] = ", round(volTot[ik],2)," sample size = ",nsizes[zk]), xlab = colnames(vols[[zk]][ik]), xlim = c(xmin,xmax))
      hist(vols[[zk]][,ik], main = paste0("Whole area vol [m^3] = ", round(volTot[ik],2)," sample size = ",nsizes[zk]), xlab = colnames(vols[[zk]][ik]), xlim = c(xmin,xmax))
    }
    #    hist(vols[[2]][,ik], main = paste0("Whole area vol [m^3/ha] = ", round(volTot[ik],2)," sample size = 20000"), xlab = colnames(vols[[1]][ik]), xlim = c(xmin,xmax))
    #    hist(vols[[3]][,ik], main = paste0("Whole area vol [m^3/ha] = ", round(volTot[ik],2)," sample size = 40000"), xlab = colnames(vols[[1]][ik]), xlim = c(xmin,xmax))
    #    hist(vols[[3]][,ik], main = paste0("Whole area vol [m^3/ha] = ", round(volTot[ik],2)," sample size = 40000"), xlab = colnames(vols[[1]][ik]), xlim = c(xmin,xmax))
  }
  
  dev.off()
  
  nSamples <- nSamplesf
  nSitesRun <- nSitesRunf
}

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