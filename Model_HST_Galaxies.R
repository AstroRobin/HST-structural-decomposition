### Model_HST_Galaxies.R ###
# A script designed to take in a catalogue of HST COSMOS observations sets (containing names of the 4+ exposures) as well as a catalogue of galaxies with their positions in each exposure listed.

# The code loops over each observation set, loading into memory all the exposure frames that comprise the set and locates all sources that are contained within the bounds of that observation set.
# For each found source, the code creates cutouts and loads/creates in the necessary inputs (e.g. PSF, segmentation map) for running a model optimisation on the images.
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
## Date: 12/01/2022

## Usage:
# >> Rscript Model_HST_Galaxies.R --obset=[ID] --inputcat=[PATH] --model=[MODEL] --computer=['local'|'magnus']

#######################################
### Get arguments from command-line ###
#######################################

library(argparse)
parser = ArgumentParser()
parser$add_argument("-i","--inputcat",required=TRUE,
                    action='store',dest='inputcat',type='character',default=NULL,
                    help="The catalogue of galaxy sources over which to run the optimisation routine.",metavar="PATH")
parser$add_argument("-o","--obset",
                    action='store',dest='obset',type='character',default=NULL,
                    help="The observation set ID in the format JXXXXX010",metavar="ID")
parser$add_argument("-m","--model",
                    action='store',dest='model',type='character',default=NULL,choices=c('single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'),
                    help="The model type to be optimised.",metavar="ID")
parser$add_argument("-c","--computer",
                    action='store',dest='computer',type='character',default=NULL,choices=c('local','magnus'),
                    help="The machine on which this code is running",metavar="computer")

# Parse arguments and validate
#args = parser$parse_args()

args = list()
args$inputcat = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/DEVILS_Structural_Decomposition/Catalogues/Subcats/Test_Magnus/DEVILS_HST_subcat_J8PU01010-0.csv'
args$model = 'single'
args$computer = 'local'

if (is.null(args$computer)){
  nCores = 24#strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
} else if (args$computer == 'local'){
  nCores = strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1))
} else if (args$computer == 'magnus') {
  #.libPaths(c('/home/sbellstedt/R/x86_64-suse-linux-gnu-library/3.6',.libPaths()))
  .libPaths(c('/group/pawsey0160/rhwcook/r/3.6',.libPaths()))
  nCores = strtoi(Sys.getenv('SLURM_CPUS_PER_TASK', unset=1)) #nCores = 24
}


library(ProFound)
library(ProFit)
library(ProFuse)
# library(ProTools)

library(Rfits)
library(magicaxis)
library(Highlander)
library(celestial)

library(foreach) # Parallel evaluation of for loops
library(doParallel) # Parallel backend for the foreach's %dopar% function 

library(glue)

evalglobal = TRUE # Set global evaluate


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
  
  cat(glue(txt),end)
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
  
  binary = number2binary(dqval, numbits=15)
  for (flag in flags){
    if (binary[flag] == 1){
      return (TRUE)
    }
  }
  return (FALSE)
}


### Stack a set of images and take the median. If dqmaps is given, use this to remove cosmic rays and bad pixels.
median_stack = function(images, dqmaps=NULL, plot=FALSE){
  
  # Convert pixels with DQ flags into NA values
  if (!is.null(dqmaps)){
    for (ii in seq_along(images)){
      dqmask = apply(dqmaps[[ii]], c(1,2), mask_pixel)
      images[[ii]][dqmask == TRUE] = NA
    }
  }
  
  # Coerce list of images into a 3D matrix
  imgMat = do.call(cbind, images)
  imgMat = array(imgMat, dim=c(dim(images[[1]]), length(images)))
  
  # Get the median along the z-axis
  imgMedian = apply(imgMat, c(1, 2), median, na.rm = TRUE)
  
  if (plot){
    magplot(imgMedian)
  }
  
  return (imgMedian)
}


### Given a key, get the corresponding value from a FITS header object. Use `to` for converting to a specific data type
from_header = function(hdr, keys, to=as.numeric){
  if (is.string(keys)){ # convert to a list if only a single string given
    keys = c(keys)
  }
  
  vals = c()
  for (key in keys) {
    val = hdr[which(hdr == key)+1]
    if (!is.null(to)){
      vals = c(vals, to(val))
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



##############################
### Define some parameters ###
##############################

if (args$computer == 'local'){
  baseDir = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA'
  
  framesDir = glue('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Frames/PID09822')
  psfScriptPath = glue('{baseDir}/Programs/HST-structural-decomposition/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  psfsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs'
  
  resultsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Fitting_Outputs'
  plotDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Plots'
  
  badIDFilename = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Logs/Bad_IDs_single_test.txt'
  
} else if (args$computer == 'magnus'){
  baseDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition'
  
  framesDir = glue('{baseDir}/Data/Frames/PID09822')
  psfsDir = glue('{baseDir}/Data/PSFs')
  psfScriptPath = glue('{baseDir}/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  
  resultsDir = glue('{baseDir}/Results/ProFuse_Outputs')
  plotDir = glue('{baseDir}/Results/Plots')
  
  badIDFilename = '{baseDir}/Results/Logs/Bad_IDs_single_test.txt'
}

# Set up some parallel code
registerDoParallel(cores=nCores)
printf("INFO: Running on {args$computer} machine with [{nCores}] cores.")


if (! file.exists(args$inputcat)){
  stop(glue("ERROR: The input catalogue '{args$inputcat}' does not exist."))
}

if (!is.null(args$obset)){ # Check if observation id is valid
  if (nchar(args$obset) != 9 | substr(args$obset,7,9) != '010'){
    stop(glue("ERROR: The observation set ID '{args$obset}' does not have the correct pattern (expected nine characters ending in '010')."))
  }
}


# Open bad IDs file
write('failed_IDs',file=badIDFilename,append=FALSE)

sciIndex = 2 # The index in HST fits files that points to the science image HDU
errIndex = 3 # The index in HST fits files that points to the error  HDU
dqIndex = 4 # The index in HST fits files that points to the data quality HDU
pixScaleHST = 0.05 # The pixel scale of HST ACS non-rotated and raw .flt exposure frames.

toPlot = FALSE#TRUE # Whether to show plots
toSave = TRUE # Whether to save plots

cutoutBox = c(400,400) # The size of the box used to cut out images of sources
skyBoxDims = c(75,75)

roughFit = FALSE


#################################################
### Select Light Profile Model and Parameters ###
#################################################

modelName = args$model #'single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'

modelOpts = list(disk_nser = 1, # Disk Sersic indec always starts with exponential
                   bulge_nser = 4, # Bulge Sersic indec always starts with de Vaucouleurs
                   bulge_circ = TRUE) # A circular bulge (i.e., axrat = 1))

if (modelName == 'single'){ # Single component Sersic (pure-disks and Ellipticals)
  modelOpts = c(modelOpts, 
                list(n_comps = 1,
                     sing_nser_fit = TRUE,
                     disk_nser_fit = FALSE,
                     bulge_nser_fit = FALSE))
} else if (modelName == 'psf-exp'){ # PSF bulge + Exponential disk
  modelOpts = c(modelOpts, 
                list(n_comps = 1.5,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = FALSE,
                    bulge_nser_fit = FALSE))
} else if (modelName == 'deV-exp'){  # de Vaucouleurs bulge + Exponential disk (Canonical light profiles)
  modelOpts = c(modelOpts, 
                list(n_comps = 2,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = FALSE,
                    bulge_nser_fit = FALSE))
} else if (modelName == 'psf-free'){  # PSF bulge + Free disk
  modelOpts = c(modelOpts, 
                list(n_comps = 1.5,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = TRUE,
                    bulge_nser_fit = FALSE))
} else if (modelName == 'deV-free'){  # de Vaucouleurs bulge + Free disk
  modelOpts = c(modelOpts, 
                list(n_comps = 2,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = TRUE,
                    bulge_nser_fit = FALSE))
} else if (modelName == 'free-exp'){  # Free bulge + Exponential disk (pseudobulges?)
  modelOpts = c(modelOpts, 
                list(n_comps = 2,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = FALSE,
                    bulge_nser_fit = TRUE))
} else if (modelName == 'free-free'){  # Free bulge + Free disk (Most flexibile model, least constrained)
  modelOpts = c(modelOpts, 
                list(n_comps = 2,
                    sing_nser_fit = FALSE,
                    disk_nser_fit = TRUE,
                    bulge_nser_fit = TRUE))
} else {
  printf("ERROR: Model name '{modelName}' not recognised.")
}

optimOpts = list(Niters=c(100,500),NfinalMCMC=1000,
                 CMA = list(control=list(maxit=100)),
                 LD = list(control= list(abstol=0.1), Iterations=1e3, Status=1e2, Algorithm = 'CHARM', Thinning = 1, Specs=list(alpha.star=0.44)))


###########################################
### Load Source Catalogue and Exposures ###
###########################################
dfCat = read.csv(args$inputcat) # read in the catalogue

if (is.null(args$obset)){ # If no observation set ID given, pull this directly from the catalogue
  args$obset = dfCat$obset_id[1]
}

cat(glue("Obsevation Set ID: {args$obset}\n"))

### Load each exposure's image & DQ map from the observation set
obsetDir = glue('{framesDir}/{args$obset}')
expFiles = list.files(path=obsetDir)

imgList = list() # The science images
dqList = list() # The data quality maps
hdrList = list() # The headers for the fits files
for ( jj in seq_along(expFiles) ){
  hdulist = Rfits_read_all(glue( "{obsetDir}/{expFiles[jj]}" ), pointer=FALSE, header=TRUE)
  
  expName = substr(expFiles[jj],1,9)
  imgList[[expName]] = list(hdulist[[sciIndex]]$imDat, hdulist[[sciIndex+3]]$imDat)
  dqList[[expName]] = list( profoundMakeSegimDilate(segim=hdulist[[dqIndex]]$imDat,size=3,expand=c(4096,8192))$segim, # Expand the cosmic rays by a pixel either side (i.e. size=3)
                            profoundMakeSegimDilate(segim=hdulist[[dqIndex+3]]$imDat,size=3,expand=c(4096,8192))$segim )
  hdrList[[expName]] = hdulist[[sciIndex]]$hdr
}

# Get all DEVILS sources that are present in this exposure
foundSources = which(toupper(substr(dfCat$name_exp1,1,6)) == substr(args$obset,1,6))

#Example foreach from sbellstedt:  foreach(ii=1:length(IDs$CATAID), .inorder=FALSE)%dopar%{ # .packages='packages needed for execution', .export/.noexport, .combine='c','+','*','rbind','cbind'
foreach(ii=1:length(foundSources))%dopar%{ #for (idx in foundSources){
  #check=try({ # Catch any failed runs
  
  idx = foundSources[[ii]] # The row index for this source
  cat(idx,'\n')
  
  sourceName = format(dfCat$UID[idx],scientific=FALSE)
  print(sourceName)
  # Create output plot file if it does not already exist
  if (!file.exists(glue("{plotDir}/{sourceName}"))){
    dir.create(glue("{plotDir}/{sourceName}"))
  }
  
  numExps = min(c(dfCat$num_exps[idx],4))
  cat(glue("INFO: UID: {sourceName}\nINFO: {dfCat$num_exps[idx]} exposures found (using {numExps}): \n\n"))
  
  imgCutList = list()
  dqCutList = list()
  maskList = list()
  psfList = list()
  for (jj in seq(1,numExps)){
    
    expName = dfCat[[glue('name_exp{jj}')]][idx]
    expExt = dfCat[[glue('chip_exp{jj}')]][idx] - (-1)^dfCat[[glue('chip_exp{jj}')]][idx]
    srcLoc = c(dfCat[[glue('x_exp{jj}')]][idx], dfCat[[glue('y_exp{jj}')]][idx])
      
    cat('    > ',jj,': ',expName,' [',expExt,']\n', sep='')
    
    imgCutList[[expName]] = magcutout(imgList[[expName]][[expExt]], loc=srcLoc, box=cutoutBox, plot=FALSE)$image
    dqCutList[[expName]] = magcutout(dqList[[expName]][[expExt]], loc=srcLoc, box=cutoutBox, plot=FALSE)$image
    maskList[[expName]] = apply(apply(dqCutList[[expName]], c(1,2), mask_pixel), c(1,2), as.numeric )
    
    psfFilename = glue('{psfsDir}/{sourceName}_{dfCat$name_exp1[idx]}_psf.fits')
    if (!file.exists(psfFilename)){
      system(glue("Rscript {psfScriptPath} {psfsDir} {sourceName}_{dfCat$name_exp1[idx]} {dfCat$chip_exp1[idx]} {dfCat$x_exp1[idx]} {dfCat$y_exp1[idx]} 0.0 13.0"))
    }
    psfList[[expName]] = Rfits_read_image(psfFilename,ext=1)$imDat
    
  }
  
  printf("INFO: Creating median stack of {length(imgCutList)} cutout images.")
  imgCutMedian = median_stack(images=imgCutList, dqmaps=dqCutList, plot=FALSE)
  
  ### Run profoundProfound on the median stack of cutouts:
  printf("INFO: Creating segmentation map for median stacked image.")
  seg = profoundProFound(image=imgCutMedian, sigma=1.5, dilete=5, skycut=1.2, SBdilate=2.5, tol=7, reltol=2, ext=7, plot=toPlot)
  
  ### Set up Found2Fit objects for each exposure
  printf("INFO: Generating found2Fits objects for each exposure cutout.")
  if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}_{modelName}_frame_cutouts.png"), width=720, height=720, pointsize=16)}
  if (toPlot|toSave){par(mfrow=c(2,2), mar = c(1,1,1,1))}
  found2Fits = list()
  
  for (jj in seq(1,numExps)){
    ### Read required info from header.
    magZP = 35.796 #from_header(hdrList[[jj]], 'PHOTZPT', to=as.numeric) 
    ccdGain = 1 #from_header(hdrList[[jj]], 'CCDGAIN', to=as.numeric) # This is in the previous HDU, but it's always 1 for HST ACS
    
    expName = dfCat[[glue('name_exp{jj}')]][idx]
    expExt = dfCat[[glue('chip_exp{jj}')]][idx] - (-1)^dfCat[[glue('chip_exp{jj}')]][idx]
    srcLoc = c(dfCat[[glue('x_exp{jj}')]][idx],dfCat[[glue('y_exp{jj}')]][idx]) # need to move this out of the loop and just use the first exposure's position as the location.
    
    
    found2Fits[[jj]] = profuseFound2Fit(imgList[[expName]][[expExt]], psf=psfList[[expName]], segim = seg$segim,
                                        mask=maskList[[expName]],SBdilate=2.5,reltol=1.5,ext=5,
                                        loc=srcLoc, cutbox=cutoutBox, loc_use=FALSE, loc_fit=TRUE,
                                        magzero=magZP, gain = ccdGain, 
                                        Ncomp=modelOpts$n_comps,
                                        sing_nser_fit=modelOpts$sing_nser_fit,
                                        disk_nser_fit=modelOpts$disk_nser_fit,  disk_nser=modelOpts$disk_nser,
                                        bulge_nser_fit=modelOpts$bulge_nser_fit, bulge_nser=modelOpts$bulge_nser, bulge_circ=modelOpts$bulge_circ,
                                        pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
                                        tightcrop = FALSE, deblend_extra=FALSE, fit_extra=FALSE, rough=roughFit, plot=toPlot|toSave)

    found2Fits[[jj]]$mask = maskList[[expName]] # Why doesn't the mask get passed through Found2Fit into the $Data structures?
    
    if (jj > 1){
      # Offsets calculated as the difference between the x/y values listed in the data frame.
      found2Fits[[jj]]$Data$offset = c(0,0)
      #found2Fits[[jj]]$Data$offset = c(dfCat[[glue('x_exp{jj}')]][idx] - dfCat[['x_exp1']][idx],
      #                                 dfCat[[glue('y_exp{jj}')]][idx] - dfCat[['y_exp1']][idx])
      
    }
    
    # Need to extend model beyond calc region
    found2Fits[[jj]]$Data$usecalcregion = FALSE
    
    if(toPlot|toSave){
      text(0.1*cutoutBox[1],0.925*cutoutBox[2],as.character(dfCat[[glue('name_exp{jj}')]][idx]), adj=0, col='white', cex=1.75)
    }
    
  }
  
  if(toSave){dev.off()}
  
  dataList = lapply(found2Fits, `[[`, 'Data') # For each `found2Fits` instance, pull out the `Data` object
  
  # Set additional parameters to the Data list for fitting multiple images.
  dataList$mon.names = found2Fits[[1]]$Data$mon.names
  dataList$parm.names = found2Fits[[1]]$Data$parm.names
  dataList$N = found2Fits[[1]]$Data$N
  
  printf("INFO: Running Highlander optimisation of model.")
  highFit = Highlander(found2Fits[[1]]$Data$init, Data=dataList, likefunc=profitLikeModel)
                       #CMAargs = optimOpts$CMA, LDargs = optimOpts$LD)
  

  if (toPlot|toSave){
    for (jj in seq(1,numExps)){
      if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-exp{jj}_{modelName}_model_residuals.png"), width=1000, height=650, pointsize=20)}
      profitLikeModel(highFit$parm, Data=dataList[[jj]], makeplots=TRUE, rough=FALSE, plotchisq=TRUE)
      if (toSave){dev.off()}
    }
  }
  
  
  remakeModel = profitRemakeModellist(parm = highFit$parm, Data = dataList[[1]])
  optimModellist = remakeModel$modellist
  
  ### Save plots of resulting model ###
  if(toSave & modelOpts$n_comps>1){ # fix to add fake bulge if single component chosen.
    png(filename = glue("{plotDir}/{sourceName}/{sourceName}_{modelName}_radial_profile.png"), width=720, height=480, pointsize=16)
    profitEllipsePlot(Data=dataList[[1]],modellist=optimModellist,pixscale=pixScaleHST,SBlim=25)
    dev.off()
  }
  
  ### Save outputs to .RDS file ###
  printf("INFO: Creating segmentation map for median stacked image.")
  outputs = list(idx=idx,
                 sourceName=sourceName,
                 modelName=modelName,
                 modelOptions=modelOpts,
                 optimOptions=optimOpts,
                 cutoutBox=cutoutBox,
                 exposureFiles=expFiles,
                 dqList=dqCutList,
                 imgList=imgCutList,
                 imgMedian=imgCutMedian,
                 hdrList=hdrList,
                 maskList=maskList,
                 psfList=psfList,
                 segMap=seg,
                 found2FitList=found2Fits,
                 dataList=dataList,
                 highFit=highFit,
                 optimModellist=optimModellist)
  
  dir.create(file.path(resultsDir, sourceName), showWarnings = FALSE) # Create output directory if one does not already exist
  outputFilename = glue('{resultsDir}/{sourceName}/{sourceName}_{modelName}_output.rds')
  saveRDS(object=outputs, file=outputFilename)
  
  #}) # end try
  #if (class(check)=='try-error'){
  #  printf("ERROR: Galaxy {ii} failed to be fit.")
  #  write(ii,file=badIDFilename,append=TRUE)
  #} else {
  #  printf("ERROR: Galaxy {idx} successfully fit.")
  #}
  
} # end foreach()


### Run in RStudio afterwards ###
#result = readRDS('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Fitting_Outputs/101505979732574752/101505979732574752_deV-exp_output.rds')

### Running dead-Simple Queue
# dSQ --job-file .../Jobs/[jobfile] --batch-file [path/batchfile-YYYY-MM-DD.sh] --job-name HSTdecomp_SS --output .../Jobs/Outputs/jobfile_%A_%a-%N.out --status-dir .../Jobs/Status/ --account=pawsey0160 --nodes=2 --time=24:00:00 --mem=64GB --mail-type=ALL --mail-user=robin.cook@uwa.edu.au

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

# The old way of initiating the found2Fits
### Produce the sigma map
#skyMap = profoundProFound(imgCutList[[jj]], mask = maskList[[jj]], segim = seg$segim, size=5,
#                          box = c(80,80), type='bicubic',
#                          magzero=magZP, gain=ccdGain, #pixscale=pixScale, # header=header,
#                          stats=TRUE, plot=FALSE)

#sigma = profoundMakeSigma(imgCutList[[jj]],sky=skyMap$sky,skyRMS=skyMap$skyRMS,objects=seg$objects,mask=maskList[[jj]],gain=ccdGain,plot=FALSE)

#found2Fits[[jj]] = profuseFound2Fit(imgCutList[[jj]], segim=seg$segim, sigma=sigma, psf=psfList[[jj]],
#                                    magzero=magZP, gain = ccdGain, 
#                                    Ncomp=n_comps,
#                                    sing_nser_fit=sing_nser_fit,bulge_circ=bulge_circ,
#                                    disk_nser_fit=disk_nser_fit,  disk_nser=disk_nser,
#                                    bulge_nser_fit=bulge_nser_fit, bulge_nser=bulge_nser, bulge_circ=bulge_circ,
#                                    pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
#                                    tightcrop = FALSE, deblend_extra=FALSE, fit_extra=FALSE, rough=roughFit, plot=toPlot)

    