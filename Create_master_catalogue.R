# Loading all the required packages 
sessionInfo()
.libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))
#.libPaths(c("/group/pawsey0160/jthorne/R/magnus/",.libPaths()))
#.libPaths(c('/group/pawsey0160/rhwcook/r/3.6',.libPaths()))

library(ProFuse)
library(magicaxis)
library(data.table)
library(foreach)
library(doParallel)

library(LaplacesDemon)
library(celestial)

#library(scales)

power = function(base, exp) sign(base) * abs(base)^exp

machine = 'zeus'
if (machine == 'zeus'){
  nCores = 20
} else {
  nCores = 24
}

registerDoParallel(cores=nCores)

#outFilename = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Tests/Alpha_test_master_catalogue.csv'
outFilename = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/Catalogues/DEVILS_D10_structural_decomposition_catalogue_v1.1.csv'

#outputsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs/Alpha'
outputsDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/ProFuse_Outputs'

#catFilename = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Catalogues/Subcats/DEVILS_HST-ACS_exposure_positions.csv'
#catFilename = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/DEVILS_HST-ACS_exposure_positions.csv'
#catFilename = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Jobs/Status/Alpha_fitting_status.csv'
catFilename = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/Status/alpha/Alpha_fitting_status.csv'

logFile = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition/Results/Catalogues/Master_catalogue_log.txt'

dfCat = fread(catFilename, colClasses=c("UID"="character")) # Catalogue
dfSample = dfCat[dfCat$hasfit == TRUE] #dfCat[dfCat$masstot > 10^10]

models = c('single', 'psf-exp', 'free-exp')
modelNames = list(single='single', `psf-exp`='psfexp', `free-exp`='serexp')

paramRef = list(single=c('sersic.xcen', 'sersic.ycen', 'sersic.mag', 'sersic.re', 'sersic.nser', 'sersic.ang', 'sersic.axrat'),
                `psf-exp`=c('sersic.xcen', 'sersic.ycen', 'pointsource.mag', 'sersic.mag', 'sersic.re', 'sersic.ang', 'sersic.axrat'),
                `free-exp`=c('sersic.xcen1', 'sersic.ycen1', 'sersic.mag1', 'sersic.mag2', 'sersic.re1', 'sersic.re2', 'sersic.nser1', 'sersic.ang2', 'sersic.axrat2'))

unlogRef = list()
for (model in models){
  unlogRef[[model]] = c()
  for (ii in seq_along(paramRef[[model]])){
    if (grepl('re', paramRef[[model]][ii]) | grepl('nser', paramRef[[model]][ii]) | grepl('axrat', paramRef[[model]][ii])){
      unlogRef[[model]][ii] = TRUE
    } else {
      unlogRef[[model]][ii] = FALSE
    }
  }
}

colNames = c('UID', 'nexps', 'npixels')
for (model in models){
  sfx = modelNames[[model]]
  if (model == 'single'){
    colNames = c(colNames, paste(c('xcen', 'ycen', 'mag', 're', 'nser', 'ang', 'axrat'), sfx, sep='_'))
    colNames = c(colNames, paste(c('xcen_err', 'ycen_err', 'mag_err', 're_err', 'nser_err', 'ang_err', 'axrat_err'), sfx, sep='_'))
  } else if (model == 'psf-exp'){
    colNames = c(colNames, paste(c('xcen', 'ycen', 'mag1', 'mag2', 're2', 'ang2', 'axrat2'), sfx, sep='_'))
    colNames = c(colNames, paste(c('xcen_err', 'ycen_err', 'mag1_err', 'mag2_err', 're2_err', 'ang2_err', 'axrat2_err'), sfx, sep='_'))
  } else if (model == 'free-exp'){
    colNames = c(colNames, paste(c('xcen', 'ycen', 'mag1', 'mag2', 're1', 're2', 'nser1', 'ang2', 'axrat2'), sfx, sep='_'))
    colNames = c(colNames, paste(c('xcen_err', 'ycen_err', 'mag1_err', 'mag2_err', 're1_err', 're2_err', 'nser1_err', 'ang2_err', 'axrat2_err'), sfx, sep='_'))
  }
  
  colNames = c(colNames, paste(c('LP', 'DIC', 'chisq_perdof', 'time'), sfx, sep='_'))
}

#galList = list.files(outputsDir)

sample = foreach(ii = 1:nrow(dfSample), .combine='rbind') %dopar% {
  
  sourceName = dfSample$UID[ii]
  #sourceName = galList[ii]
  
  fileConn = file(logFile, open='a')
  writeLines(paste0(ii,': ',sourceName), fileConn )
  close(fileConn)
  
  outputFile = paste0(outputsDir,'/',sourceName,'/',sourceName,'_output.rds')
  if (file.exists(outputFile)){
    result = readRDS(outputFile)
  } else {
    return (NULL)
  }
  
  # sourceName = as.character(sapply(strsplit(file_name, "\\."), "[[", 1))
  
  nExps = result$numimgs
  nPix = result$dataList$single$N
  
  outList = list(sourceName, nExps, nPix) # initialise the list of outputs
  for (model in models){
    dof = length(result$highFit[[model]]$parm)
    
    ### Using LaplacesDemon MCMC chains
    #bestParams = as.numeric(result$highFit[[model]]$LD_last$Summary1[paramRef[[model]],'Mean'])
    #bestParams[unlogRef[[model]]] = 10^bestParams[unlogRef[[model]]]
    #bestParamErrors = result$highFit[[model]]$LD_last$Summary1[paramRef[[model]],'SD']
    
    ### Using ProFit monitored parameter values
    bestParams = as.numeric(result$highFit[[model]]$parm[paramRef[[model]]])
    bestParams[unlogRef[[model]]] = 10^bestParams[unlogRef[[model]]]
    bestParamErrors = as.numeric(result$highFit[[model]]$error)
    
    LP = result$highFit[[model]]$LD_last$Summary1['LP','Mean']
    DIC = result$highFit[[model]]$LD_last$DIC1[3]
    chisq = result$highFit[[model]]$RedChi2
    time = result$highFit[[model]]$time
    
    outList = c(outList, as.numeric(bestParams), as.numeric(bestParamErrors), LP, DIC, chisq, time)
  }
  
  
  # for (model in models){
  #   
  #   # Finding the singular set of model parameters that gives the best Likelihood
  #   posteriors = data.table( cbind(result$highFit[[model]]$LD_last$Posterior1, result$highFit[[model]]$LD_last$Monitor) )
  #   posteriors = posteriors[order(-LP_mon),] # sort posteriors in descending order by the LP_mon (= LP)
  #   
  #   dof = length(result$highFit[[model]]$parm)
  #   chisqCut = max(posteriors$LP_mon) - qchisq(0.68, df=dof)/2 # i.e., one standard deviation of a chi-squared distriubtion (with n degs of freedom) from the best LP
  #   bestPosteriors = posteriors[posteriors$LP_mon > chisqCut,]
  #   bestPosteriors[, LP_mon:=NULL] # now remove LP_mon from table as it's no longer needed
  #   
  #   bestParams = as.numeric(posteriors[1,1:dof])
  #   bestParamRanges = foreach(ii = 1:ncol(bestPosteriors), .combine = 'rbind') %do% {range(bestPosteriors[,..ii])}
  #   
  #   ## iteration information
  #   type = ifelse(result$highFit[[model]]$best == 'CMA', 1, 2)
  #   iter = result$highFit[[model]]$error
  #   
  #   outList = c(outList, bestParams, as.numeric(bestParamRanges[,1]), as.numeric(bestParamRanges[,2]), type, iter)
  #    
  # }
  
  return(unlist(outList))
  
}

# covert to data.table
sample = data.table(sample)
names(sample) = colNames

write.csv(sample, file = outFilename, row.names = F)
