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

# Author: R. H. W. Cook
# Date: 12/01/2022

library(ProFound)
library(ProFit)
library(ProFuse)
# Could replace with library(ProTools)

library(Rfits)
library(magicaxis)
library(Highlander)
library(celestial)

library(glue)
library(assertthat)

evalglobal = TRUE # Set global evaluate

###############################
### Define Useful Functions ###
###############################

### Prints a formatted string, much like python's print(f"")
printf = function(txt,end=''){ 
  
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
### Define User parameters ###
##############################

workDir = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/DEVILS_Structural_Decomposition'
framesDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PID09822/Frames'
psfsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs'

outDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Fitting_Outputs'
plotDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Plots'

sciIndex = 2 # The index in HST fits files that points to the science image HDU
errIndex = 3 # The index in HST fits files that points to the error  HDU
dqIndex = 4 # The index in HST fits files that points to the data quality HDU
pixScaleHST = 0.05 # The pixel scale of HST ACS non-rotated and raw .flt exposure frames.

toPlot = TRUE#TRUE # Whether to show plots
toSave = FALSE # Whether to save plots

cutoutBox = c(400,400) # The size of the box used to cut out images of sources
skyBoxDims = c(75,75)

roughFit = FALSE

#################################################
### Select Light Profile Model and Parameters ###
#################################################

modelName = 'single'#'single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'

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


optimOpts = list(Niters=c(500,1000),NfinalMCMC=2500,
                 LDargs=list(control= list(abstol=0.1), Iterations=1e3, Status=1e2, Algorithm = 'CHARM', Thinning = 1, Specs=list(alpha.star=0.44)))

### Load Catalogue of DEVILS sources ###
framePositionsFile = glue('{workDir}/Catalogues/DEVILS_D10_HST_ACS_Frame_Positions_PID09822.csv')
dfDEVILS = read.csv(framePositionsFile)

### Load Catalogue of Frames ### 
obsetsFile = glue('{workDir}/Catalogues/Obset_Reference_Tables/HST-ACS_COSMOS_PID09822_obset_reference_table.csv')
obsetsDF = read.csv(obsetsFile)
#View(obsetsDF)


### Loop over all observation sets ###
for (ii in seq(1,1)){ #nrow(obsetsDF)){
  
  cat(glue("Obsevation Set ID: {obsetsDF$dataset[ii]}\n"))
  
  ### Load each exposure's image & DQ map from the observation set
  obsetDir = glue('{framesDir}/{obsetsDF$dataset[ii]}')
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
  
  
  #if (toPlot){par(mfrow=c(2,4))} # For plotting, set the shape of the subplots
  
  # Get all DEVILS sources that are present in this exposure
  foundSources = which(toupper(substr(dfDEVILS$exp_name_1,1,6)) == substr(obsetsDF$dataset[ii],1,6))
  for (idx in foundSources[4:4]){
    cat(idx,'\n')
    sourceName = format(dfDEVILS$UID[idx],scientific=FALSE)
    cat(glue("INFO: UID: {sourceName}\n  {dfDEVILS$found_exps[idx]} exposures found: \n\n"))
    
    imgCutList = list()
    dqCutList = list()
    maskList = list()
    psfList = list()
    for (jj in seq(1,dfDEVILS$found_exps[idx])){
      
      expName = dfDEVILS[[glue('exp_name_{jj}')]][idx]
      expChip = dfDEVILS[[glue('exp_chip_{jj}')]][idx] - (-1)^dfDEVILS[[glue('exp_chip_{jj}')]][idx]
      srcLoc = c(dfDEVILS[[glue('exp_x_{jj}')]][idx],dfDEVILS[[glue('exp_y_{jj}')]][idx])
        
      cat('    > ',jj,': ',expName,' [',expChip,']\n', sep='')
      
      imgCutList[[expName]] = magcutout(imgList[[expName]][[expChip]], loc=srcLoc, box=cutoutBox, plot=FALSE)$image
      dqCutList[[expName]] = magcutout(dqList[[expName]][[expChip]], loc=srcLoc, box=cutoutBox, plot=FALSE)$image
      maskList[[expName]] = apply(apply(dqCutList[[expName]], c(1,2), mask_pixel), c(1,2), as.numeric )
      psfList[[expName]] = Rfits_read_image(glue('{psfsDir}/101506350912562912/101506350912562912_j8pu01w4q_psf.fits',ext=2))$imDat
      #psfList[[expName]] = Rfits_read_image(glue('{psfsDir}/{sourceName}/{sourceName}_{dfDEVILS$exp_name_1[idx]}_psf.fits',ext=2))$imDat
      
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
    for (jj in seq(1,dfDEVILS$found_exps[idx])){
      ### Read required info from header.
      magZP = 35.796 #from_header(hdrList[[jj]], 'PHOTZPT', to=as.numeric) 
      ccdGain = 1 #from_header(hdrList[[jj]], 'CCDGAIN', to=as.numeric) # This is in the previous HDU, but it's always 1 for HST ACS
      
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
      
      expChip = dfDEVILS[[glue('exp_chip_{jj}')]][idx] - (-1)^dfDEVILS[[glue('exp_chip_{jj}')]][idx]
      srcLoc = c(dfDEVILS[[glue('exp_x_{jj}')]][idx],dfDEVILS[[glue('exp_y_{jj}')]][idx]) # need to move this out of the loop and just use the first exposure's position as the location.
      
      
      found2Fits[[jj]] = profuseFound2Fit_(imgList[[jj]][[expChip]], psf=psfList[[jj]], segim = seg$segim,
                                          mask=maskList[[jj]],SBdilate=2.5,reltol=1.5,ext=5,
                                          loc=srcLoc, cutbox=cutoutBox, loc_use=FALSE, loc_fit=TRUE,
                                          magzero=magZP, gain = ccdGain, 
                                          Ncomp=modelOpts$n_comps,
                                          sing_nser_fit=modelOpts$sing_nser_fit,
                                          disk_nser_fit=modelOpts$disk_nser_fit,  disk_nser=modelOpts$disk_nser,
                                          bulge_nser_fit=modelOpts$bulge_nser_fit, bulge_nser=modelOpts$bulge_nser, bulge_circ=modelOpts$bulge_circ,
                                          pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
                                          tightcrop = FALSE, deblend_extra=FALSE, fit_extra=FALSE, rough=roughFit, plot=toPlot|toSave)

      
      if (jj > 1){
        # Offsets calculated as the difference between the x/y values listed in the data frame.
        found2Fits[[jj]]$Data$offset = c(0,0)
        #found2Fits[[jj]]$Data$offset = c(dfDEVILS[[glue('exp_x_{jj}')]][idx] - dfDEVILS[['exp_x_1']][idx],
        #                                 dfDEVILS[[glue('exp_y_{jj}')]][idx] - dfDEVILS[['exp_y_1']][idx])
        
      }
      
      # Need to extend model beyond calc region
      found2Fits[[jj]]$Data$usecalcregion = FALSE
      
      if(toPlot|toSave){
        text(0.1*cutoutBox[1],0.925*cutoutBox[2],as.character(dfDEVILS[[glue('exp_name_{jj}')]][idx]), adj=0, col='white', cex=1.75)
      }
      
    }
    
    if(toSave){dev.off()}
    
    dataList = lapply(found2Fits, `[[`, 'Data') # For each `found2Fits` instance, pull out the `Data` object
    
    # Set additional parameters to the Data list for fitting multiple images.
    dataList$mon.names = found2Fits[[1]]$Data$mon.names
    dataList$parm.names = found2Fits[[1]]$Data$parm.names
    dataList$N = found2Fits[[1]]$Data$N
    
    printf("INFO: Running Highlander optimisation of model.")
    highFit = Highlander(found2Fits[[1]]$Data$init, Data=dataList, likefunc=profitLikeModel, ablim=1)
                         #LDargs = optimOpts$LD)
    
    for (jj in seq(4)){
      dataList[[jj]]$usecalcregion = FALSE
    }
    
    if (toPlot|toSave){
      if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}_{modelName}_model_residuals.png"), width=1000, height=650, pointsize=20)}
      printf("INFO: Generating plots of optimised models.")
      profitLikeModel(highFit$parm, Data=dataList, makeplots=TRUE, rough=FALSE, plotchisq=TRUE)
      if (toSave){dev.off()}
    }
    
    if (toPlot|toSave){
      for (jj in seq(1,dfDEVILS$found_exps[idx])){
        if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-exp{jj}_{modelName}_model_residuals.png"), width=1000, height=650, pointsize=20)}
        profitLikeModel(highFit$parm, Data=dataList[[jj]], makeplots=TRUE, rough=FALSE, plotchisq=TRUE)
        if (toSave){dev.off()}
      }
    }
    
    
    remakeModel = profitRemakeModellist(parm = highFit$parm, Data = dataList[[1]])
    optimModellist = remakeModel$modellist
    
    ### Save plots of resulting model ###
    if(toSave & modelOpts$n_comps>1){
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
    
    dir.create(file.path(outDir, sourceName), showWarnings = FALSE)
    outputFilename = glue('{outDir}/{sourceName}/{sourceName}_{modelName}_output.rds')
    saveRDS(object=outputs, file=outputFilename)
    
  }
  
  
}


result = readRDS('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Fitting_Outputs/101505979732574752/101505979732574752_deV-exp_output.rds')



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

    