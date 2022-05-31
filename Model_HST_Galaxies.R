### Model_HST_Galaxies.R ###
# A script designed to take in a catalogue of HST COSMOS observations sets (containing names of the 4+ exposures) as well as a catalogue of galaxies with their positions in each exposure listed.

# The code loops over each galaxy, loading into memory all the exposure frames that it is contained within, choosing only those with sufficient pixels to cover the galaxy.
# For each galaxy, the code creates cutouts and loads/creates in the necessary inputs (e.g. PSF, segmentation map) for running a model optimisation on the images.
# Using profuseFound2Fit, the Datalist object containing all 4+ exposures are collated and Highlander is then run to optimise the chosen model(s) types.

# Model types may be:
#   - 'single'
#   - 'psf-exp'
#   - 'deV-exp'
#   - 'psf-free'
#   - 'deV-free'
#   - 'free-exp'
#   - 'free-free'

## Author: R. H. W. Cook
## Date: 22/02/2022

## Usage:
# >> Rscript Model_HST_Galaxies.R --inputcat=[PATH] --models=[MODEL] --computer=['local'|'magnus'] --name=[NAME]

#######################################
### Get arguments from command-line ###
#######################################

library(argparse)
parser = ArgumentParser()
parser$add_argument("-i","--inputcat",required=TRUE,
                    action='store',dest='inputcat',type='character',default=NULL,
                    help="The catalogue of galaxy sources over which to run the optimisation routine.",metavar="PATH")
parser$add_argument("-m","--models",
                    action='store',dest='models',type='character',nargs='+',default=NULL,choices=c('single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'),
                    help="The model type to be optimised.",metavar="ID")
parser$add_argument("-c","--computer",
                    action='store',dest='computer',type='character',default=NULL,choices=c('local','magnus','zeus'),
                    help="The machine on which this code is running",metavar="computer")
parser$add_argument("-n","--name",
                    action='store',dest='name',type='character',default='test',
                    help="The name of this particular fitting run.",metavar="NAME")

# Parse arguments and validate
args = parser$parse_args()

#args = list(inputcat = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/Subcats/Pilot/DEVILS_HST_pilot_subcat_0.csv', models = 'single', computer = 'local', name = 'test')
#args$inputcat = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/Subcats/DEVILS_HST_pilot_sample_sorted.csv'


if (is.null(args$computer)){
  nCores = strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
} else if (args$computer == 'local'){
  nCores = strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
} else if (args$computer == 'magnus') {
  .libPaths(c('/group/pawsey0160/rhwcook/r/3.6',.libPaths()))
  nCores = 24 #strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
} else if (args$computer == 'zeus') {
  .libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))
  nCores = 20 #strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
}

cat(paste0(".libPaths(): ",.libPaths(),"\n"))
sessionInfo()

library(ProFound)
library(ProFit)
library(ProFuse)
# library(ProTools)

library(Rfits)
library(magicaxis)
library(Highlander)
library(celestial)
library(Rwcs)

library(RColorBrewer)
library(Cairo)

library(foreach) # Parallel evaluation of for loops
library(doParallel) # Parallel backend for the foreach's %dopar% function 

library(glue)
library(assertthat)

evalglobal = TRUE # Set global evaluate
options(stringsAsFactors=FALSE) # adjust string problems
options(scipen=999) # suppress scientific notation

###############################
### Define Useful Functions ###
###############################

### Prints a formatted string, much like python's print(f"")
printf = function(txt,end='\n'){
  
  for (sub in regmatches(txt, gregexpr("\\{\\K[^{}]+(?=\\})", txt, perl=TRUE))[[1]] ){
    splitstr = strsplit(sub,':')[[1]]
    val = splitstr[[1]]
    if (length(splitstr) > 1){
      fmt = splitstr[[2]]  
    } else {
      fmt = NULL
    }
    
    if (!is.null(fmt)){
      ndecimal = as.integer(substr(fmt,2,2))
      txt = gsub(sub, paste0('format(round(',val,', ',ndecimal,'), nsmall=',ndecimal,')'), txt)
    }
  }
  
  cat(glue(txt),sep=end)
}


### Covert numbers to binary vectors 
number2binary = function(number, numbits) {
  
  binary_vector = as.numeric(intToBits(number))
  
  if(missing(numbits)) {
    return(binary_vector)
  } else {
    return(binary_vector[1:(length(binary_vector) - numbits)])
  }
}


### Check whether a pixel should be masked based on the value in the corresponding Data Quality map.
mask_pixel = function(dqval, flags=c(3,5,8,10,13)){
  
  if (is.na(dqval)){
    return(TRUE)
  }
  
  binary = number2binary(dqval, numbits=15)
  for (flag in flags){
    if (binary[flag] == 1){
      return (TRUE)
    }
  }
  return (FALSE)
}


### Given a key, get the corresponding value from a FITS header object. Use `as` for converting to a specific data type after retrieving (e.g. as=as.numeric)
from_header = function(hdr, keys, as=NULL){
  if (is.string(keys)){ # convert to a list if only a single string given
    keys = c(keys)
  }
  
  vals = c()
  for (key in keys) {
    val = hdr[which(hdr == key)+1]
    if (!is.null(as)){
      vals = c(vals, as(val))
    } else {
      vals = c(vals, val)
    }
  }
  
  if (length(vals) == 1){
    return (vals[1])
  } else {
    return (vals)
  }
  
}


### For the list of pre-defined model types, set the appropriate model options.
get_model_opts = function(model){
  modelOpts = list(disk_nser = 1, # Disk Sersic indec always starts with exponential
                   bulge_nser = 4, # Bulge Sersic indec always starts with de Vaucouleurs
                   bulge_circ = TRUE) # A circular bulge (i.e., axrat = 1))
  
  if (model == 'single'){ # Single component Sersic (pure-disks and Ellipticals)
    modelOpts = c(modelOpts, 
                  list(n_comps = 1,
                       sing_nser_fit = TRUE,
                       disk_nser_fit = FALSE,
                       bulge_nser_fit = FALSE))
  } else if (model == 'psf-exp'){ # PSF bulge + Exponential disk
    modelOpts = c(modelOpts, 
                  list(n_comps = 1.5,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = FALSE,
                       bulge_nser_fit = FALSE))
  } else if (model == 'deV-exp'){  # de Vaucouleurs bulge + Exponential disk (Canonical light profiles)
    modelOpts = c(modelOpts, 
                  list(n_comps = 2,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = FALSE,
                       bulge_nser_fit = FALSE))
  } else if (model == 'psf-free'){  # PSF bulge + Free disk
    modelOpts = c(modelOpts, 
                  list(n_comps = 1.5,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = TRUE,
                       bulge_nser_fit = FALSE))
  } else if (model == 'deV-free'){  # de Vaucouleurs bulge + Free disk
    modelOpts = c(modelOpts, 
                  list(n_comps = 2,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = TRUE,
                       bulge_nser_fit = FALSE))
  } else if (model == 'free-exp'){  # Free bulge + Exponential disk (pseudobulges?)
    modelOpts = c(modelOpts, 
                  list(n_comps = 2,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = FALSE,
                       bulge_nser_fit = TRUE))
  } else if (model == 'free-free'){  # Free bulge + Free disk (Most flexibile model, least constrained)
    modelOpts = c(modelOpts, 
                  list(n_comps = 2,
                       sing_nser_fit = FALSE,
                       disk_nser_fit = TRUE,
                       bulge_nser_fit = TRUE))
  } else {
    cat(paste0("ERROR: Model name '",model,"' not recognised.\n"))
  }
  
  return(modelOpts)
  
}



##############################
### Define some parameters ###
##############################

# Define filepaths
date = Sys.Date()
if (args$computer == 'local'){
  baseDir = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA' # Base directory
  
  framesDir = glue('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Frames') # Where the frame FITS files are kept
  psfScriptPath = glue('{baseDir}/Programs/HST-structural-decomposition/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  psfsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs' # Where to load/save PSFs
  focusFilename = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/PSFs/HST-ACS_focus_lookup_RJM.csv' # The focus measurements table for HST PSFs
  
  resultsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs' # ProFuse outputs saved here
  plotDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Plots' # Plots saved here
  
  badIDFilename = glue('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Status/Failed_IDs_{args$name}_{date}.txt') # Which galaxies failed
  goodIDFilename = glue('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Status/Successful_IDs_{args$name}_{date}.txt')
  logsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Logs'
  
  source('/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Programs/HST-structural-decomposition/Model_plots.R')
  
} else if (args$computer == 'magnus' | args$computer == 'zeus'){
  baseDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition'
  
  framesDir = glue('{baseDir}/Data/Frames')
  psfsDir = glue('{baseDir}/Data/PSFs')
  psfScriptPath = glue('{baseDir}/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  focusFilename = glue('{baseDir}/Data/HST-ACS_focus_lookup_RJM.csv')
  
  resultsDir = glue('{baseDir}/Results/ProFuse_Outputs')
  plotDir = glue('{baseDir}/Results/Plots')
  
  badIDFilename = glue('{baseDir}/Results/Status/Failed_IDs_{args$name}_{date}.txt')
  goodIDFilename = glue('{baseDir}/Results/Status/Successful_IDs_{args$name}_{date}.txt')
  logsDir = glue('{baseDir}/Results/Logs/{args$name}')
  
  source('{baseDir}/Scripts/Model_plots.R')
  
}

if (!file.exists(logsDir)){
  cat("INFO: Logs directory did not exist, creating now.\n")
  dir.create(logsDir)
}


# Set up some parallel code
registerDoParallel(cores=nCores)
cat(paste0("INFO: Running on ",args$computer," machine with [",nCores,"] cores.\n"))

if (! file.exists(args$inputcat)){
  stop(paste0("ERROR: The input catalogue '",args$inputcat,"' does not exist."))
}

pids = c('PID09822', 'PID10092', 'PID12440')

pixScaleHST = 0.05 # The pixel scale of HST ACS non-rotated and raw .flt exposure frames.
ccdGain = 1

sizePSF = 5.0 # arcsec
refocusPSF = TRUE # whether the PSF should be refocused using the HST focus measurements from R. J. Massey
trimPSF = TRUE # whether to trim the TinyTim PSF to have pixel dimensions that match the size of sizePSF based on pixelScaleHST
clobberPSF = FALSE # IF TRUE, always creates a new PSF using TinyTim; else a new PSF will only be created if one does not already exist.

maxExps = 16 # The maximum number of exposures to use for a single fit

cutoutBox = c(601,601) # The size of the box used to cut out images of sources
skyBoxDims = c(150,150)

roughFit = FALSE

toPlot = FALSE #TRUE # Whether to show plots
toSave = TRUE # Whether to save plots
logsToFile = FALSE # Whether print statements should be sent to separate files (good when running parallel computing)

######################################
### Select Optimisation Parameters ###
######################################

optimOpts = list(optim_iters=5, Niters=c(200,200), NfinalMCMC=500, walltime=7.5*60) # max walltime of 7.5 hours
#optimOpts = list(optim_iters=3, Niters=c(100,100), NfinalMCMC=300, walltime=7.5*60)
optimOpts$CMA = list(control=list(maxit=optimOpts$Niters[1]))
optimOpts$LD = list(control=list(abstol=0.1), Iterations = optimOpts$Niters[2], Status=50, Algorithm = 'CHARM', Thinning = 1, Specs=list(alpha.star=0.44))

#############################
### Load Source Catalogue ###
#############################
cat("INFO: Loading source catalogue\n")
dfCat = read.csv(args$inputcat, colClasses=c("UID"="character")) # read in the catalogue

if (refocusPSF){ # load the PSF focuses file (if available)
  if (file.exists(focusFilename)){
    dfFocii = read.csv(focusFilename)
    dfFocii = dfFocii[!is.nan(dfFocii$Focus),] # clean any rows that are missing focus values
  } else {
    cat(paste0("WARNING: File does not exist '",focusFilename,"'\n"))
    refocusPSF = FALSE
  }
}

#for (idx in seq(nrow(dfCat))){
foreach(idx=1:nrow(dfCat), .inorder=FALSE) %dopar% {
  check=try({ # Catch any failed runs
    
    startTime = Sys.time()
    
    # The row index for this source
    sourceName = format(dfCat$UID[idx], scientific=FALSE)
    
    if (logsToFile == TRUE){
      cat(paste0("\n\n",sourceName," (",idx,"/",nrow(dfCat),")\n"))
      logFile = paste0(logsDir,'/',sourceName,'_log.txt')
      sink(logFile)
    } else {
      logFile = ''
    }
    
    cat("\n\n#############################################################################\n")
    cat(paste0("## INFO: Running structural decomposition on ",sourceName," (",idx,"/",nrow(dfCat),") ##\n"))
    cat("#############################################################################\n")
    
    numExps = min(c(dfCat$num_exps[idx], maxExps))
    if (numExps == 0){
      cat(paste0("INFO: No frames found for source: ",sourceName,"\n"))
      write(sourceName,file=badIDFilename,append=TRUE)
    } else {
      
      cat(paste0("INFO: ",dfCat$num_exps[idx]," exposures found (using ",numExps,"): \n\n"))
      
      # Create output plot file if it does not already exist
      if (!file.exists(glue("{plotDir}/{sourceName}"))){
        cat("INFO: Output directory did not exist, creating now.\n")
        dir.create(paste0(plotDir,'/',sourceName))
      }
      
      ### Load each exposure's image & DQ map from the observation set
      expNameList = list() # The names of each exposure
      magzpList = list() # The magnitude zeropoints
      gainList = list() # The CCD gain, typically 1.0 bu is 2.0 in some PIDs.
      
      imgCutList = list()
      hdrCutList = list()
      dqCutList = list()
      maskList = list()
      psfList = list()
      
      srcLocList = list()
      imgWarpList = list()
      maskWarpList = list()
      
      skyList = list()
      skyRMSList = list()
      
      for (jj in seq(1,numExps)){
        expName = dfCat[[paste0('name_exp',jj)]][[idx]]
        expNameList[[jj]] = expName
        obsetID = paste0( substr(toupper(expName), 1, 6), '010' )
        
        # Establish which proposal ID this exposure came from.
        if (substr(expName,1,3) == 'j8p'){
          proposalID = 'PID09822'
        } else if (substr(expName,1,3) == 'j8x'){
          proposalID = 'PID10092'
        } else if (substr(expName,1,3) == 'jco'){
          proposalID = 'PID13641'
        } else if (substr(expName,1,3) == 'jbo'){
          proposalID = 'PID12440'
        } else if (substr(expName,1,3) == 'jbh'){
          proposalID = 'PID12328'
        } else {
          cat(paste0("WARNING: No proposal ID could be found for exposure '",expName,"'\n"))
        }
        
        
        # The FITS extensions indices for the chip
        expChip = dfCat[[paste0('chip_exp',jj)]][idx]
        sciIndex = ifelse(expChip == 1, 2+3, 2) # The index in HST fits files that points to the science image HDU
        errIndex = ifelse(expChip == 1, 3+3, 3) # The index in HST fits files that points to the error HDU
        dqIndex = ifelse(expChip == 1, 4+3, 4) # The index in HST fits files that points to the data quality HDU
        
        # Load the image, data quality (dq) map and hdr objects
        imgFilename = paste0(framesDir,"/",proposalID,"/",obsetID,"/",expName,"_flt.fits")
        cat(paste0("\n\nINFO: Loading exposure #",jj," on chip ",expChip,": ",imgFilename,"\n"))
        hdulist = Rfits_read_all(imgFilename, pointer=FALSE, zap=c('LOOKUP','DP[1-2]'))
        
        img = hdulist[[sciIndex]] # this is an Rfits_image object
        dq = hdulist[[dqIndex]]
        
        # Get source position info in exposure
        srcLocList[[expName]] = Rwcs_s2p(RA=dfCat[['RAcen']][idx], Dec=dfCat[['DECcen']][idx], header = img$raw, pixcen='R')
        
        # Calculate magnitude zeropoint
        fluxlambda = from_header(hdr = img$hdr, keys = 'PHOTFLAM', as=as.numeric)
        pivotlambda = from_header(hdr = img$hdr, keys = 'PHOTPLAM', as=as.numeric)
        expTime = from_header(hdr = hdulist[[1]]$hdr, keys = 'EXPTIME', as=as.numeric) # exposure time from first extension header
        expStart = from_header(hdr = hdulist[[1]]$hdr, keys = 'EXPSTART', as=as.numeric) # exposure start date (MJD) from first extension header
        gainList[[jj]] = from_header(hdr = hdulist[[1]]$hdr, keys = 'CCDGAIN', as=as.numeric) # CCD gain from first extension header
        magzpList[[jj]] = -2.5*log10(fluxlambda/expTime) - 5*log10(pivotlambda) - 2.408 # instrumental zeropoint magnitudes in AB mags (from https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints)
        
        rm(hdulist)
        
        # Now get cutouts and WCS-adjusted headers
        cat(paste0("INFO: Source location in ",expName,": ",srcLocList[[jj]][1],", ",srcLocList[[jj]][2],"\n"))
        imgCutList[[expName]] = img[srcLocList[[jj]][1], srcLocList[[jj]][2], box=cutoutBox]
        hdrCutList[[expName]] = imgCutList[[expName]]$hdr
        
        dqCutList[[expName]] = dq[srcLocList[[jj]][1], srcLocList[[jj]][2], box=cutoutBox] # magcutout(dq, loc=srcLocList[[jj]], box=cutoutBox, plot=FALSE)$image
        dqDilated = profoundMakeSegimDilate(segim=dqCutList[[expName]]$imDat,size=3,expand=seq(4096,8192)) # Expand the cosmic rays by a pixel either side (i.e. size=3)
        maskList[[expName]] = apply(apply(dqDilated$segim, c(1,2), mask_pixel), c(1,2), as.integer) # Convert the binary dq map to a list of 1s and 0s
        #maskList[[expName]][is.na(imgCutList[[expName]])] = 1 # Also set missing values (NaNs) to be masked
        
        # Calculate the initial sky statistics of each exposure, so that they can be inverse-variance stacked later
        skyGrid = profoundProFound(image=imgCutList[[expName]]$imDat, mask=maskList[[expName]], box=skyBoxDims, roughpedestal=TRUE)
        
        # Warp images and dqs onto a common WCS.
        cat(paste0("INFO: Warping image and dq map to a common WCS\n"))
        imgWarpList[[expName]] = Rwcs_warp(image_in=imgCutList[[jj]], header_out = imgCutList[[1]]$raw)
        skyList[[expName]] = Rwcs_warp(image_in=skyGrid$sky, header_in = imgCutList[[jj]]$raw, header_out = imgCutList[[1]]$raw)$imDat
        skyRMSList[[expName]] = Rwcs_warp(image_in=skyGrid$skyRMS, header_in = imgCutList[[jj]]$raw, header_out = imgCutList[[1]]$raw)$imDat
        
        maskWarpList[[expName]] = apply(round(Rwcs_warp(image_in=maskList[[jj]], header_in = imgCutList[[jj]]$raw, header_out = imgCutList[[1]]$raw)$imDat, digits=0), c(1,2), as.integer)
        maskWarpList[[expName]][is.na(maskWarpList[[expName]])] = 1 # NAs produced by unoccupied regions in the warped mask should be converted to 1, i.e. masked
        
        # Load or Generate PSF
        psfFilename = paste0(psfsDir,'/',sourceName,'_',expName,'_psf.fits')
        if (!file.exists(psfFilename) | clobberPSF==TRUE){
          cat(paste0("INFO: Creating new PSF, saving to file '",psfFilename,"'.\n"))
          focus = 0.0
          if (refocusPSF){
            sel = which.min(abs(expStart-dfFocii$MJD))
            if (abs(expStart - dfFocii$MJD[sel]) <= 0.125){ # Focus measurement must be within a MJD diff of 0.125 = 3 hours (i.e. ~2 HST orbits)
              focus = dfFocii$Focus[sel]
              cat(paste0("INFO: Refocusing PSF by ",focus," micrometers.\n"))
            } else {
              cat(paste0("WARNING: Could not refocus PSF, closest focus measurement was ", expStart - dfFocii$MJD[sel], " days from observation.\n"))
            }
          }
          
          #>> Rscript Generate_HST_PSFs.R [dir] [name] [chip] [x] [y] [focus=0.0] [size=10.0] [filter=f814w] [spectrum=13]
          cat(paste0("INFO: Running TinyTim with:\n>> Rscript ",psfScriptPath," ",psfsDir," ",sourceName,"_",expName," ",expChip," ",round(srcLocList[[jj]][1])," ",round(srcLocList[[jj]][2])," ",focus," ",sizePSF,"\n"))
          if (logsToFile) {sink()}
          invisible(
            system(paste0("Rscript ",psfScriptPath," ",psfsDir," ",sourceName,"_",expName," ",expChip," ",round(srcLocList[[jj]][1])," ",round(srcLocList[[jj]][2])," ",focus," ",sizePSF))
          )
          if (logsToFile) {sink(logFile,append=TRUE)}
        }
        
        psf = Rfits_read_image(psfFilename,ext=1)$imDat # read in the PSF file
        
        if (trimPSF == TRUE){ # Trim the PSF?
          cat("INFO: Trimming PSF\n")
          dimsPSF = sizePSF/pixScaleHST
          peakPSF = which(psf==max(psf, na.rm=T) , arr.ind = T)
          psfList[[expName]] = psf[(peakPSF[1]-dimsPSF%/%2):(peakPSF[1]+dimsPSF%/%2), (peakPSF[2]-dimsPSF%/%2):(peakPSF[2]+dimsPSF%/%2)]
        } else {
          psfList[[expName]] = psf
        }
        
        cat(paste0("INFO: Dimensions of ",expName," PSF: ",dim(psfList[[expName]])[1]," x ",dim(psfList[[expName]])[2],"\n"))
        
      }
      
      cat("\nINFO: Exposure names:", end=' ')
      cat(names(imgCutList), end='\n')
      cat(paste0("INFO: Length of images: ",length(imgCutList),"\n"))
      cat(paste0("INFO: Length of dq: ",length(dqCutList),"\n"))
      cat(paste0("INFO: Length of mask: ",length(maskList),"\n"))
      cat(paste0("INFO: Length of PSF: ",length(psfList),"\n"))
      
      # Stack the images and masks
      cat(paste0("INFO: Creating median stack of ",length(imgCutList)," cutout images.\n"))
      imgCutStack = profoundMakeStack(image_list=lapply(imgWarpList, `[[`, 'imDat'), sky_list=skyList, skyRMS_list=skyRMSList, mask_list=maskWarpList, magzero_in = unlist(magzpList), masking='&')$image
      
      maskStack = maskList[[1]]
      if (numExps > 1){
        for (ii in seq(2,numExps)){
          maskStack = maskStack & maskList[[ii]]
        }
      }
      maskStack = apply(maskStack, c(1,2), as.integer)
      
      ### Run profoundProfound on the median stack of cutouts:
      cat("INFO: Creating segmentation map for median stacked image.\n")
      seg = profoundProFound(image=imgCutStack, mask=maskStack, sigma=1.5, size=9, skycut=0.9, SBdilate=2.5, tolerance=7, reltol=0.3, ext=9, redoskysize=55, roughpedestal = T, plot=F)
      
      if (toPlot|toSave){
        if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'-segmentation_map.png'), width=650, height=650, pointsize=20)}
        profoundSegimPlot(image=imgCutStack, segim=seg$segim, mask=maskStack)
        if (toSave){dev.off()}
      }
      
      ### Subtract sky
      cat("INFO: Measuring sky statistics for each exposure and subtracting from cutout images.\n")
      segimList = list()
      imgList = list() # The sky subtracted image cutouts
      hdrList = list()
      for (jj in seq(1,numExps)){
        expName = dfCat[[paste0('name_exp',jj)]][idx]
        
        segimList[[expName]] = apply(round(Rwcs_warp(image_in=seg$segim, header_in = imgCutList[[1]]$raw, header_out = imgCutList[[jj]]$raw, doscale = FALSE, interpolation='nearest')$imDat, digits=0), c(1,2), as.integer)
        segimList[[expName]] = profoundMakeSegimDilate(segim=segimList[[expName]],size=3)$segim # This dilation will catch any aliasing holes caused by warping integers maps onto different (typically rotated) grids
        segimList[[expName]][is.na(segimList[[expName]])] = 0
        
        if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'-',expNameList[jj],'_sky_statistics.png'), width=720, height=720, pointsize=16)}
        if (toPlot|toSave){par(mfrow=c(2,2), mar = c(1,1,1,1))}
        skyList[[expName]] = profoundProFound(image=imgCutList[[jj]]$imDat, segim=segimList[[expName]], mask=maskList[[jj]], magzero = magzpList[[jj]], gain=ccdGain, plot=FALSE,
                                              sigma=1.25, skycut=0.5, SBdilate=3, size=19, box=skyBoxDims)
        
        if (toPlot|toSave){
          plot(skyList[[expName]])
        }
        
        if(toSave){dev.off()}
        
        # Subtract sky
        imgList[[expName]] = imgCutList[[jj]]$imDat - skyList[[jj]]$sky
        hdrList[[expName]] = imgCutList[[jj]]$raw
        
      }
      
      ## Save images of the image cutouts + overlaid segmentation and mask maps
      if (toPlot|toSave){
        if (numExps <= 4){
          if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_cutouts.png'), width=900, height=900, pointsize=10)}
          par(mfrow=c(2,2), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
        } else if (numExps > 4 & numExps <= 8) {
          if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_cutouts.png'), width=1800, height=900, pointsize=10)}
          par(mfrow=c(2,4), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
        } else if (numExps > 8 & numExps <= 12){
          if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_cutouts.png'), width=1800, height=900+450, pointsize=10)}
          par(mfrow=c(3,4), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
        } else {
          if (toSave){png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_cutouts.png'), width=1800, height=1800, pointsize=10)}
          par(mfrow=c(4,4), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
        }
        
        for (jj in seq(1,numExps)){
          profoundSegimPlot(image=imgList[[jj]], segim=segimList[[jj]], mask=maskList[[jj]], bad=NA)
          text(0.1*cutoutBox[1],0.925*cutoutBox[2],expNameList[jj], adj=0, col='white', cex=1.5*(0.5 + 0.75 * numExps%/%4))
        }
        if (toSave){dev.off()}
      }
      
      ### Calculate offsets and rotations between frames
      cat("INFO: Calculating offsets/rotations between exposures based on header information.\n")
      offsets = list() # The offsets are just the sub-pixel offsets between positions as the cutouts are already centered on the pixel positions.
      xrots = list() # rotations from the x-axis
      yrots = list() # rotations from the y-axis (!= xrots + 90, because HST images have skew)
      
      cds = from_header(hdrCutList[[1]], keys = c('CD1_1','CD2_1','CD1_2','CD2_2'), as=as.numeric)
      cdsign = sign(det(matrix(cds, ncol=2))) # Get the sign of the determinant for the CD matrix
      
      xrot1 = rad2deg(atan2(-cdsign*cds[2], cds[1]))
      yrot1 = rad2deg(atan2(cdsign*cds[3], cds[4]))
      for (jj in seq(1,numExps)){
        
        cds = from_header(hdrCutList[[jj]], keys = c('CD1_1','CD2_1','CD1_2','CD2_2'), as=as.numeric)
        cdsign = sign(det(matrix(cds, ncol=2))) # Get the sign of the determinant for the CD matrix
        xrots[[jj]] = rad2deg(atan2(-cdsign*cds[2], cds[1]))
        yrots[[jj]] = rad2deg(atan2(cdsign*cds[3], cds[4]))
        
        #offsets[[jj]] = c(dfCat[[glue('x_exp{jj}')]][idx]%%1 - dfCat[['x_exp1']][idx]%%1, dfCat[[glue('y_exp{jj}')]][idx]%%1 - dfCat[['y_exp1']][idx]%%1, yrots[[jj]] - yrot1)
        offsets[[jj]] = c(srcLocList[[jj]][1]%%1 - srcLocList[[1]][1]%%1, srcLocList[[jj]][2]%%1 - srcLocList[[1]][2]%%1, yrots[[jj]] - yrot1)
        
      }
      
      cat("INFO: Estimating WCS offsets between exposures.\n")
      segStatsList = list()
      segExpList = list()
      wcsOffsets = list()
      maxDiffs = list(x=list(), y=list())
      cenDiffs = list(x=list(), y=list())
      
      # if (toSave){CairoPNG(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_offsets.png'), width=3000, height=1500*(1+(numExps-1)%/%4), pointsize=20)}
      if (toPlot){
        plotArr = matrix(seq(1,8), nrow = 2)
        if (numExps > 4){
          for (ii in seq((numExps-1)%/%4)){ # extra rows
            plotArr = rbind(plotArr, matrix(seq(8*ii+1,8*(ii+1)), nrow = 2))
          }
        }
        layout(plotArr, widths = rep(1/4,4), heights=rep(1/2*(1+(numExps-1)%/%4),2*(1+(numExps-1)%/%4)))
      }
      
      # This code below is a bit janky, needs some polishing.
      for (jj in seq(1,numExps)){
        # Run profound to get segment statistics
        imgWarpList[[jj]]$imDat[is.nan(imgWarpList[[jj]]$imDat)] = NA
        segExpList[[jj]] = profoundProFound(image=imgWarpList[[jj]]$imDat, header=imgWarpList[[jj]]$raw, segim=seg$segim, mask=maskWarpList[[1]] | maskWarpList[[jj]], plot=F,
                                            sky=0, skyRMS = skyRMSList[[jj]], magzero=magzpList[[jj]], pixscale=pixScaleHST, boundstats = T)
        segStatsList[[jj]] = segExpList[[jj]]$segstats
        tempseg = segExpList[[jj]]$segim
        
        if (jj != 1){
          xmatch = coordmatch(segStatsList[[1]][,c("RAcen","Deccen")], segStatsList[[jj]][,c("RAcen","Deccen")], rad=0.1, radunit="asec")$bestmatch # 0.1 arcsec = 2 pixels
          
          inboth = c()
          for (ii in seq_along(xmatch[,1])){ # check if the segment actually exists within each exposure
            #if (anyNA(imgWarpList[[1]]$imDat[segExpList[[1]]$segim == xmatch[ii,'refID']]) | anyNA(imgWarpList[[jj]]$imDat[segExpList[[jj]]$segim == xmatch[ii,'compareID']])){
            inboth[ii] = ifelse(!is.na(segStatsList[[1]][xmatch$refID[ii],'mag']) & !is.na(segStatsList[[jj]][xmatch$compareID[ii],'mag']), TRUE, FALSE)
          }
          
          tempseg[!tempseg %in% segStatsList[[jj]][xmatch[inboth,2],'segID']] = 0
          
          maxDiffs$x[[jj]] = segStatsList[[jj]][xmatch[inboth,2],"xmax"] - segStatsList[[1]][xmatch[inboth,1],"xmax"]
          maxDiffs$y[[jj]] = segStatsList[[jj]][xmatch[inboth,2],"ymax"] - segStatsList[[1]][xmatch[inboth,1],"ymax"]
          cenDiffs$x[[jj]] = segStatsList[[jj]][xmatch[inboth,2],"xcen"] - segStatsList[[1]][xmatch[inboth,1],"xcen"]
          cenDiffs$y[[jj]] = segStatsList[[jj]][xmatch[inboth,2],"ycen"] - segStatsList[[1]][xmatch[inboth,1],"ycen"]
          
          if (toPlot){
            magbin(maxDiffs$x[[jj]], maxDiffs$y[[jj]], xlim=c(-15,15), ylim=c(-15,15), step=1, dustlim=0, xlab='x offset [pixels]', ylab='y offset [pixels]', cex.axis=1.5, cex.lab=1.5)
            points(maxDiffs$x[[jj]], maxDiffs$y[[jj]], xlim=c(-15,15), ylim=c(-15,15), xlab='x offset [pixels]', ylab='y offset [pixels]', cex.axis=1.5, cex.lab=1.5)
            points(x=offsets[[jj]][1], y=offsets[[jj]][2], pch=4, col='red', cex=2) 
          }
          
          if (length(maxDiffs$x[[jj]][abs(maxDiffs$x[[jj]]) < 5]) >= 10 & length(maxDiffs$y[[jj]][abs(maxDiffs$y[[jj]]) < 5]) >= 10){ # Must have at least 10 points to calculate average offsets from
            wcsOffsets[[jj]] = c(median(maxDiffs$x[[jj]][abs(maxDiffs$x[[jj]]) < 5]), median(maxDiffs$y[[jj]][abs(maxDiffs$y[[jj]]) < 5]))
          } else {
            wcsOffsets[[jj]] = c(0, 0)
          }
          
          cat(paste0("WCS offsets in ",expNameList[[jj]],": x=",wcsOffsets[[jj]][1],", y=",wcsOffsets[[jj]][2],"\n"))
          
          if (toPlot){
            points(x=wcsOffsets[[jj]][1], y=wcsOffsets[[jj]][2], pch=3, col='magenta', cex=2)
            text(-15,15,as.character(dfCat[[paste0('name_exp',jj)]][idx]), adj=0, col='black', cex=1.5+0.5*(numExps-1)%/%4)
            text(-15,13, paste0('dx = ',wcsOffsets[[jj]][1],'; dy = ',wcsOffsets[[jj]][2]), adj=0, col='black', cex=1.5+0.5*(numExps-1)%/%4)
            
            profoundSegimPlot(image=imgWarpList[[jj]]$imDat, header=imgWarpList[[jj]]$raw, segim=tempseg, mask=maskWarpList[[1]] | maskWarpList[[jj]],
                              sky=0, skyRMS = skyRMSList[[jj]], magzero=magzpList[[jj]], pixscale=pixScaleHST, boundstats = T)
            points(x=segStatsList[[jj]][xmatch[inboth,2],"xcen"], y=segStatsList[[jj]][xmatch[inboth,2],"ycen"], pch=3, col='magenta',cex=2)
            points(x=segStatsList[[1]][xmatch[inboth,1],"xcen"], y=segStatsList[[jj]][xmatch[inboth,1],"ycen"], pch=4, col='yellow',cex=2)
          }
          
        } else {
          if (toPlot){
            plot(x=offsets[[jj]][1], y=offsets[[jj]][2], pch=4, col='red', cex=2, xlim=c(-15,15), ylim=c(-15,15), xlab='x offset [pixels]', ylab='y offset [pixels]', cex.axis=1.75, cex.lab=1.75)
            text(-15,15,as.character(dfCat[[paste0('name_exp',jj)]][idx]), adj=0, col='black', cex=1.5+0.5*(numExps-1)%/%4)
            
            profoundSegimPlot(image=imgWarpList[[jj]]$imDat, header=imgWarpList[[jj]]$raw, segim=tempseg, mask=maskWarpList[[1]] | maskWarpList[[jj]],
                              sky=0, skyRMS = skyRMSList[[jj]], magzero=magzpList[[jj]], pixscale=pixScaleHST, boundstats = T)
          }
          
          maxDiffs$x[[jj]] = NA
          maxDiffs$y[[jj]] = NA
          cenDiffs$x[[jj]] = NA
          cenDiffs$y[[jj]] = NA
          wcsOffsets[[jj]] = c(0, 0)
        }
        
        
      }
      
      # if (toSave){dev.off()}
      
      # ## Correct for minor WCS offsets if they exist
      # for (jj in seq(1,numExps)){
      #   offsets[[jj]][1] = offsets[[jj]][1] + wcsOffsets[[jj]][1]
      #   offsets[[jj]][2] = offsets[[jj]][2] + wcsOffsets[[jj]][2]
      # }
      
      ### Set up Found2Fit objects for each exposure
      outputs = list(args=args,
                     idx=idx,
                     row=dfCat[idx,],
                     models=args$models,
                     optimOptions=optimOpts,
                     cutoutBox=cutoutBox,
                     skyBoxDims=skyBoxDims,
                     magzpList=magzpList,
                     sourceName=sourceName,
                     numimgs=length(expNameList),
                     expNames=expNameList,
                     imgList=imgList,
                     hdrList=hdrList,
                     imgStack=imgCutStack,
                     dqList=dqCutList,
                     maskList=maskList,
                     skyList=skyList,
                     psfList=psfList,
                     segimList=segimList,
                     segStatsList=segStatsList,
                     offsets=offsets,
                     wcsOffsets=wcsOffsets,
                     wcsDiffs=list(max=maxDiffs, cen=cenDiffs)
      )
      
      for (model in args$models){ #'single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'
        modelOpts = get_model_opts(model)
        cat(paste0("INFO: Optimising model: '",model,"', with options:\n"))
        print(modelOpts)
        
        cat("INFO: Generating found2Fits objects for each exposure cutout.\n")
        dataList = profuseMultiImageFound2Fit(image_list = imgList, psf_list = psfList, segim_list = segimList, mask_list = maskList, offset_list = offsets,
                                              # SBdilate=2.5, reltol=1.5, ext=5,
                                              magzero=unlist(magzpList), gain = gainList,
                                              Ncomp=modelOpts$n_comps,
                                              sing_nser_fit=modelOpts$sing_nser_fit,
                                              disk_nser_fit=modelOpts$disk_nser_fit,  disk_nser=modelOpts$disk_nser,
                                              bulge_nser_fit=modelOpts$bulge_nser_fit, bulge_nser=modelOpts$bulge_nser, bulge_circ=modelOpts$bulge_circ,
                                              pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
                                              rough=roughFit, tightcrop = TRUE, deblend_extra=FALSE, fit_extra=FALSE, plot=toPlot)
        
        
        # Rename F2F objects according to their exposure name
        for (jj in seq(1, numExps)){
          names(dataList)[names(dataList) == paste0("image",jj)] = expNameList[jj]
        }
        
        cat(paste0("INFO: Running Highlander optimisation for model '",model,"'.\n"))
        cat(paste0("Optim iters: ",optimOpts$optim_iters,"\nNiters: ",optimOpts$Niters,"\nNfinalMCMC: ",optimOpts$NfinalMCMC,"\nwalltime: ",optimOpts$walltime,"\n\n"))
        highFit = profuseMultiImageDoFit(image_list = imgList, dataList, ablim=1,
                                         optim_iters=optimOpts$optim_iters, Niters=optimOpts$Niters, NfinalMCMC=optimOpts$NfinalMCMC, walltime=optimOpts$walltime, 
                                         CMAargs = optimOpts$CMA, LDargs = optimOpts$LD)
        
        if (toPlot|toSave){
          if (toSave) {CairoPNG(filename=paste0(plotDir,'/',sourceName,'/',sourceName,'_',model,'_multifit.png'), width=(numExps+1)*325+60, height=1000, pointsize=20)}
          multiImageMakePlots(datalist=dataList, imagestack=imgCutStack, segim=seg$segim, parm=highFit$parm, magzps = magzpList,
                              plottext=paste0("UID:\n",sourceName,"\n\n # Exposures: ",dataList$Nim,"\n log(Mstar): ",format(round(log10(dfCat[['StellarMass']][idx]), 2), nsmall = 2),"\n z: ",format(round(dfCat[['zBest']][ii], 3), nsmall = 2)))
          # modelnum = 2
          # if (toSave) {CairoPNG(filename=paste0('/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Tests/',result$sourceName,'_',result$models[modelnum],'_multifit.png'), width=(length(result$expNames)+1)*325+50, height=1000, pointsize=20)}
          # multiImageMakePlots(datalist=result$dataList[[modelnum]], imagestack=result$imgStack, segim=result$segimList[[1]], parm=result$highFit[[modelnum]]$parm, magzps = result$magzpList,
          #                     plottext=paste0("UID:\n",result$sourceName,"\n\n # Exposures: ",result$numimgs,"\n log(Mstar): ",format(round(log10(result$row[['StellarMass']][1]), 2), nsmall = 2),"\n z: ",format(round(result$row[['zBest']][1], 3), nsmall = 2)))
          
          if (toSave){dev.off()}
        }
        
        
        remakeModel = profitRemakeModellist(parm = highFit$parm, Data = dataList[[1]])
        optimModellist = remakeModel$modellist
        
        ### Save plots of resulting model ###
        if(toSave & modelOpts$n_comps>=2){ # fix to add fake bulge if single component chosen.
          png(filename = paste0(plotDir,'/',sourceName,'/',sourceName,'_',model,'_radial_profile.png'), width=720, height=480, pointsize=16)
          profitEllipsePlot(Data=dataList[[1]], modellist=optimModellist, pixscale=pixScaleHST, SBlim=25, fwhm=0.1) # 0.1 arcsec is the typical HST ACS PSF width at F814W from the COSMOS Survey 
          dev.off()
        }
        
        # Add model optimisation results to outputs list
        outputs$dataList[[model]]=dataList
        outputs$highFit[[model]]=highFit
        outputs$modelOpts[[model]]=modelOpts
        outputs$optimModellist[[model]]=optimModellist
        
      } # END model loop
      
      endTime = Sys.time()
      elapsedTime = (endTime - startTime)[[1]]
      outputs$elapsedTime = elapsedTime
      
      ### Save outputs to .RDS file ###
      dir.create(file.path(resultsDir, sourceName), showWarnings = FALSE) # Create output directory if one does not already exist
      outputFilename = paste0(resultsDir,'/',sourceName,'/',sourceName,'_output.rds')
      cat(paste0("INFO: Saving outputs to '",outputFilename,"'\n"))
      saveRDS(object=outputs, file=outputFilename)
      
    } # END numexps condition
    
  }) # END trycatch
  if (class(check)=='try-error'){
    endTime = Sys.time()
    elapsedTime = (endTime - startTime)[[1]]
    
    cat(paste0("ERROR: source ",sourceName," (",idx,") failed to be fit after ",format(round(elapsedTime, 2), nsmall = 2)," minutes.\n"))
    write(sourceName,file=badIDFilename,append=TRUE)
    dev.off()
    if (logsToFile) {sink()}
  } else {
    cat(paste0("INFO: source ",sourceName," (",idx,") successfully fit after ",format(round(elapsedTime, 2), nsmall = 2)," minutes.\n"))
    if (logsToFile) {sink()}
    write(sourceName,file=goodIDFilename,append=TRUE)
  }
  
  
} # END main galaxy loop


### Run in RStudio afterwards ###
#result = readRDS('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs/101505979732574752/101505979732574752_deV-exp_output.rds')

##### Model HST Galaxies Overview #####
### Find all sources that exist with this exposure ###

### Loop over all sources ###

### Load data for each of the N exposures ###
# Load image
# Load dqmask
# Load PSF

### Create sky map ###
# Run ProFound() w/ fairly relaxed parameters
# Subtract sky from image

### Create segmentation map ###
# Run ProFound() with conservative detection
# Dilate and Stitch segmentation map

### Set up ProFuse fit ###
# Set up model choice (e.g. Ncomp, bulge_circ, disk_nser, etc.)
# Get 4x Found-to-Fits (FtoF) by running profuseFound2Fit( ) or profuseDoFit( )
# Set FtoF_{2-4}$Data$offset c(hdr{2-4}['CRPIX1']-hdr1['CRPIX1'], hdr{2-4}['CRPIX2']-hdr1['CRPIX2'])
# Also set rotations
# Combine datalists using: Datalists = c(list(FtoF_1$Data), list(FtoF_2$Data), list(FtoF_3$Data), list(FtoF_4$Data))

### Run Highlander Optimisation ###
# initiate Highlander Fit with 2 iterations each of CMA and LaplacesDemon.
# Complete MCMC optimisation with a longer MCMC run to get sufficient posteriors for parameters.

### OLD CODE ###
## The old way of running Highlander
# highFit = Highlander(found2Fits[[1]]$Data$init, Data=dataList, likefunc=profitLikeModel, ablim=1, 
#                      optim_iters=optimOpts$optim_iters, Niters=optimOpts$Niters, NfinalMCMC=optimOpts$NfinalMCMC,
#                      CMAargs = optimOpts$CMA, LDargs = optimOpts$LD)

## The old way of initiating the found2Fits (without using profuseMultiImageFound2Fits)
# found2Fits0 = list()
# for (jj in seq(1,numExps)){
#   
#   expName = dfCat[[paste0('name_exp',jj)]][idx]
#   srcLoc = c(dfCat[[glue('x_exp{jj}')]][idx], dfCat[[glue('y_exp{jj}')]][idx]) # need to move this out of the loop and just use the first exposure's position as the location.
#   
#   
#   
#   found2Fits0[[jj]] = profuseFound2Fit(imgCutList[[expName]], psf=psfList[[expName]], segim = seg$segim, mask=maskList[[expName]],
#                                        SBdilate=2.5,reltol=1.5,ext=5,
#                                        magzero=magZP, gain = ccdGain, 
#                                        Ncomp=modelOpts$n_comps,
#                                        sing_nser_fit=modelOpts$sing_nser_fit,
#                                        disk_nser_fit=modelOpts$disk_nser_fit,  disk_nser=modelOpts$disk_nser,
#                                        bulge_nser_fit=modelOpts$bulge_nser_fit, bulge_nser=modelOpts$bulge_nser, bulge_circ=modelOpts$bulge_circ,
#                                        pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
#                                        rough=roughFit, tightcrop = FALSE, deblend_extra=FALSE, fit_extra=FALSE, plot=toPlot|toSave)
#   
#   if (jj > 1){# Offsets calculated as the difference between the x/y values listed in the data frame.
#     found2Fits0[[jj]]$Data$offset = c(dfCat[[glue('x_exp{jj}')]][idx]%%1 - dfCat[['x_exp1']][idx]%%1, # just the sub-pixel offsets need to be set as cutouts have been made centered on the pixel positions.
#                                       dfCat[[glue('y_exp{jj}')]][idx]%%1 - dfCat[['y_exp1']][idx]%%1)
#   }
#   
#   # Need to extend model beyond calc region
#   found2Fits0[[jj]]$Data$usecalcregion = FALSE
#   
#   if(toPlot|toSave){
#     text(0.1*cutoutBox[1],0.925*cutoutBox[2],as.character(dfCat[[paste0('name_exp',jj)]][idx]), adj=0, col='white', cex=1.75)
#   }
#   
# }
# 
# if(toSave){dev.off()}
# 
# dataList0 = lapply(found2Fits0, `[[`, 'Data') # For each `found2Fits` instance, pull out the `Data` object
# 
# # Set additional parameters to the Data list for fitting multiple images.
# dataList0$mon.names = found2Fits0[[1]]$Data$mon.names
# dataList0$parm.names = found2Fits0[[1]]$Data$parm.names
# dataList0$N = found2Fits0[[1]]$Data$N

