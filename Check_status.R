### Check status of HST decompositions

sessionInfo()
.libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))

library(data.table)
library(foreach) # Parallel evaluation of for loops
library(doParallel) # Parallel backend for the foreach's %dopar% function 

machine = 'zeus'
if (machine == 'zeus'){
  nCores = 20
} else {
  nCores = 24
}

registerDoParallel(cores=nCores)

statusDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/Status/alpha'
plotsDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/Plots'
outputsDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/ProFuse_Outputs'

catFilename = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Catalogues/Subcats/DEVILS_HST-ACS_exposure_positions.csv'
dfCat = fread(catFilename, colClasses = c("UID"="character"))

dfOut = foreach(ii=1:nrow(dfCat), .inorder=TRUE, .combine='rbind') %dopar% {
  if (ii %% 1000 == 0){
    cat(ii,'\n')
  }
  
  sourceID = dfCat$UID[ii]
  numExps = dfCat$num_exps[ii]
  rdsFilename = paste0(outputsDir,'/',sourceID,'/',sourceID,'_output.rds')
  
  hasFit = NA
  started = NA
  
  if (file.exists(rdsFilename)){
    hasFit = TRUE
    started = TRUE
  } else {
    hasFit = FALSE
  
    plotFilename = paste0(plotsDir,'/',sourceID,'/',sourceID,'-segmentation_map.png')
    
    if (numExps == 0 | ! file.exists(plotFilename)) {  # if no exposures found or no plots created, fit was never started
      started = FALSE
    } else {
      started = TRUE
    }
    
  }
  
  
  return( c(sourceID, numExps, hasFit, started) )

}

dfOut = data.table(dfOut)
names(dfOut) = c('UID','num_exps','hasfit','started')

outFilename = paste0(statusDir,'/Alpha_fitting_status.csv')
write.csv(dfOut, file = outFilename, row.names = F)
