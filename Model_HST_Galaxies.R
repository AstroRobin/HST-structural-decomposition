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
# >> Rscript Model_HST_Galaxies.R --inputcat=[PATH] --model=[MODEL] --computer=['local'|'magnus'] --name=[NAME]

#######################################
### Get arguments from command-line ###
#######################################

library(argparse)
parser = ArgumentParser()
parser$add_argument("-i","--inputcat",required=TRUE,
                    action='store',dest='inputcat',type='character',default=NULL,
                    help="The catalogue of galaxy sources over which to run the optimisation routine.",metavar="PATH")
parser$add_argument("-m","--model",
                    action='store',dest='model',type='character',default=NULL,choices=c('single','psf-exp','deV-exp','psf-free','deV-free','free-exp','free-free'),
                    help="The model type to be optimised.",metavar="ID")
parser$add_argument("-c","--computer",
                    action='store',dest='computer',type='character',default=NULL,choices=c('local','magnus'),
                    help="The machine on which this code is running",metavar="computer")
parser$add_argument("-n","--name",
                    action='store',dest='name',type='character',default='test',
                    help="The name of this particular fitting run.",metavar="NAME")

# Parse arguments and validate
args = parser$parse_args()

#args = list(inputcat = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/Subcats/Pilot/DEVILS_HST_pilot_subcat_0.csv', model = 'single', computer = 'local', name = 'test')
args$inputcat = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/HST_Structural_Decomposition/Catalogues/Subcats/DEVILS_HST-ACS_exposure_positions_rotations.csv'

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
library(Rwcs)

library(RColorBrewer)

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


### Plot the fitting results of multiple exposures in a single figure
profitImageScale = function(z, zlim, col = heat.colors(12),
                            breaks, axis.pos=1, axis.padj=0, add.axis=TRUE, axis.tck = 0, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos, padj=axis.padj, tck=axis.tck)}
}

multiImageMakePlots = function(datalist, imagestack, segim, parm, plottext,
                               whichcomponents=list(sersic="all",moffat="all",ferrer="all",pointsource="all"),
                               cmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(100)),
                               errcmap = rev(c("#B00000",colorRampPalette(brewer.pal(9,'RdYlBu'))(100)[2:99],"#0000B0")),
                               errischisq=FALSE, dofs=0){
  
  magZP = datalist[[1]]$magzero
  contlw = 3
  maxsigma = 5
  
  ndofs = ifelse(missing(dofs), 0, length(dofs)) # degrees of freedom
  if(ndofs>0) stopifnot(length(dofs) <= 2)
  
  parmar = c(0.5,0.5,0.25,0)
  par(mar=parmar, oma=c(1,1,1,1), mgp=c(0,0,0))
  
  # Create the plotting layout
  nrows = 5
  ncols = datalist$Nim + 1
  plotArr = cbind(c(1,2,3,4+(ncols-1)*nrows,0) ,array(seq(4, 4+(ncols-1)*nrows-1), dim=c(nrows,ncols-1)))
  layout(plotArr, widths=c(rep((1-0.05)/ncols,ncols), 0.05),heights=rep(1/5,5))
  
  ## Segmentation map plot
  if (missing(segim)){
    magimage(image=imagestack)
  } else {
    profound = profoundProFound(image=imagestack, segim=segim)
    targetIdx = profound$segim[dim(imagestack)[1]%/%2,dim(imagestack)[2]%/%2]
    print(targetIdx)
    imageCut = magcutout(image=imagestack, loc=c(profound$segstats$xcen[[targetIdx]],profound$segstats$ycen[[targetIdx]]), box=rep(profound$segstats$R90[[targetIdx]]*5,2))$image
    segimCut = magcutout(image=segim, loc=c(profound$segstats$xcen[[targetIdx]],profound$segstats$ycen[[targetIdx]]), box=rep(profound$segstats$R90[[targetIdx]]*5,2))$image
    profoundSegimPlot(image=imageCut, segim=segimCut)
  }
  
  for (ii in seq(1,datalist$Nim)){ # loop over all exposures
    image = datalist[[ii]]$image
    region = datalist[[ii]]$region
    sigma = datalist[[ii]]$sigma
    
    fitModellist = profitRemakeModellist(parm = parm, Data = datalist[[ii]], offset = datalist[[ii]]$offset)$modellist
    modelimage = profitMakeModel(modellist = fitModellist,
                                 magzero = magZP, psf = datalist[[ii]]$psf, dim = dim(image),
                                 psfdim = dim(datalist[[ii]]$psf), whichcomponents = whichcomponents)$z
    
    data = list(x=1:dim(image)[1],y=1:dim(image)[2],z=image)
    residual = image - modelimage
    
    medimg = median(abs(image[region]))/2
    maximg = max(abs(image[region]))
    
    zlims = c(0,1)
    stretch="asinh"
    stretchscale = 1/medimg
    
    if (ii==1){
      ## Model image plot
      magimage(modelimage, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
      
      segcon = magimage(1-datalist[[ii]]$segim, add=T, col=NA)
      contour(segcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
      legend('topleft', legend='Model') # change region to segim
      
      ## Colorbar
      par(mar = parmar + c(0,10,0,0), mgp=c(0,0,0))
      breaks = seq(-maxsigma, maxsigma, length.out=length(cmap)+1)
      profitImageScale(zlim=c(-maxsigma,maxsigma), col=errcmap, breaks = breaks, axis.pos=2, axis.padj=1)
      
    }
    
    par(mar=parmar, mgp=c(0,0,0))
    
    ## Data image plot
    magimage(data$z, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
    
    tempcon=magimage(1-region, add=T, col=NA)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft', legend='Data', cex=1.25)
    
    ## Data - Model image plot
    magimage(residual, stretchscale=stretchscale, stretch=stretch, lo=-maximg, hi=maximg, zlim=zlims, type='num', col=cmap)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft', legend='Data-Model', cex=1.25)
    
    
    ## Chi-squared plots
    errsign = (-1)^(image > modelimage)
    if(!errischisq){
      error = residual/sigma
    } else {
      error = errsign*sqrt(abs(sigma))
    }
    
    errmap = error
    error = error[region]
    maxerr = max(abs(error))
    stretcherr = 1/median(abs(error))
    errmap[!region & (errmap>maxerr)] = maxerr
    minerr = -maxerr
    
    errmap[!region & (errmap<minerr)] = minerr
    errmap[errmap > maxsigma] = maxsigma
    errmap[errmap < -maxsigma] = -maxsigma
    
    magimage(errmap, magmap=FALSE, zlim=c(-maxsigma,maxsigma), col=errcmap)
    contour(tempcon, add=T, drawlabels = F, levels=1, col='darkgreen', lwd=contlw)
    legend('topleft',legend=bquote(chi*"=(Data-Model)"/sigma), cex=1.25)
    
    dx = 0.1
    xlims = c(-4,4)
    x = seq(xlims[1],xlims[2],dx)
    y = hist(error,breaks=c(-(maxerr+dx),x,maxerr+dx),plot=FALSE)$count[2:length(x)]/sum(region)/dx
    
    ylims = c(min(y[y>0]), 0.5)
    y[y<=0] = ylims[1]/10
    
    vardata = var(error)
    tdof=2*vardata/(vardata-1)
    tdof=LaplacesDemon::interval(tdof,0,Inf)
    
    magplot(x[1:(length(x)-1)]+dx,y, xlim=xlims, ylim=ylims, xlab="",ylab="", xaxs="i", type="s",log="y")
    lines(x, dnorm(x), col="blue", xaxs="i")
    lines(x, dt(x,tdof), col="red", xaxs="i")
    
    labs = c(expression(chi),bquote(norm(1)),bquote(Student-T(.(signif(tdof,5)))))
    cols = c("black","blue","red")
    ltys = c(1,1,1)
    legend("bottom",legend=labs,col=cols,lty=ltys, cex=1.25)    
    
    error = log10(error^2)
    xr = range(error)
    xlims = c(-3,min(2,max(xr)))
    x = seq(xlims[1],xlims[2],dx)
    dxbin = 10^x[2:length(x)]-10^x[1:(length(x)-1)]
    y = hist(error,breaks=c(xr[1]-dx,x,xr[2]+dx), plot=FALSE)$count[2:length(x)]
    y = y/sum(y)/dxbin
    ylim = c(min(y[y>0]),10)
    y[y<=0] = ylim[1]-1
    magplot(x[1:(length(x)-1)]+dx, y, xlim=xlims,ylim=ylim, xlab="",ylab="", xaxs="i", type="s",log="y")
    xp=10^x
    lines(x, dchisq(xp,1), col="blue", xaxs="i")
    labs = c(bquote(chi^2),expression(chi^2*(1)))
    cols = c("black","blue")
    
    if(ndofs > 0){
      dofcols = c("red", "darkgreen")
      for(i in 1:length(dofs)){
        dofstr = sprintf("%.3e",dofs[i])
        lines(x, dchisq(xp,dofs[i]), col=dofcols[i], xaxs="i")
        labs = c(labs,bquote(chi^2 (.(dofstr))))
        ltys = c(ltys,1)
        cols = c(cols, dofcols[i])
      }
    }
    
    abline(v=0,lty=2,col='red')
    legend("bottomleft",legend=labs,col=cols,lty=ltys, cex=1.25)
    
  }
  
  ### Plot text
  if (!missing(plottext)){
    par(mar = c(0,0,0,0), mgp = c(3,1,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.0, y = 0.9, plottext, adj=c(0,1),
         cex = 1.5, col = "black")
  }
  
}



##############################
### Define some parameters ###
##############################

if (args$computer == 'local'){
  baseDir = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA'
  
  framesDir = glue('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Frames')
  psfScriptPath = glue('{baseDir}/Programs/HST-structural-decomposition/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  psfsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs'
  
  resultsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs'
  plotDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Plots'
  
  badIDFilename = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/Logs/Bad_IDs_single_test.txt'
  
} else if (args$computer == 'magnus'){
  baseDir = '/group/pawsey0160/rhwcook/HST_Structural_Decomposition'
  
  framesDir = glue('{baseDir}/Data/Frames')
  psfsDir = glue('{baseDir}/Data/PSFs')
  psfScriptPath = glue('{baseDir}/Generate_HST_PSFs.R') # If no PSF exists, make one using this script
  
  resultsDir = glue('{baseDir}/Results/ProFuse_Outputs')
  plotDir = glue('{baseDir}/Results/Plots')
  
  badIDFilename = glue('{baseDir}/Results/Logs/Failed_IDs_{args$model}_{args$name}_{Sys.Date()}.txt')
}

printf("LibPaths: {.libPaths()}")

# Set up some parallel code
registerDoParallel(cores=nCores)
printf("INFO: Running on {args$computer} machine with [{nCores}] cores.")

if (! file.exists(args$inputcat)){
  stop(glue("ERROR: The input catalogue '{args$inputcat}' does not exist."))
}

magZP = 35.796 #from_header(hdrCutList[[jj]], 'PHOTZPT', to=as.numeric) 
ccdGain = 1 #from_header(hdrCutList[[jj]], 'CCDGAIN', to=as.numeric) # This is in the previous HDU, but it's always 1 for HST ACS
pixScaleHST = 0.05 # The pixel scale of HST ACS non-rotated and raw .flt exposure frames.

sizePSF = 5.0 # arcsec
trimPSF = TRUE # whether to trim the TinyTim PSF to have pixel dimensions that match the size of sizePSF based on pixelScaleHST
clobberPSF = TRUE # IF TRUE, always creates a new PSF using TinyTim; else a new PSF will only be created if one does not already exist.

maxExps = 4 # The maximum number of exposures to use for a single fit

toPlot = FALSE#TRUE # Whether to show plots
toSave = TRUE # Whether to save plots

cutoutBox = c(600,600) # The size of the box used to cut out images of sources
skyBoxDims = c(200,200)

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

optimOpts = list(optim_iters=5, Niters=c(200,200), NfinalMCMC=500)
optimOpts$CMA = list(control=list(maxit=optimOpts$Niters[1]))
optimOpts$LD = list(control=list(abstol=0.1), Iterations = optimOpts$Niters[2], Status=50, Algorithm = 'CHARM', Thinning = 1, Specs=list(alpha.star=0.44))


###########################################
### Load Source Catalogue and Exposures ###
###########################################
printf("INFO: Loading source catalogue")
dfCat = read.csv(args$inputcat) # read in the catalogue

#foreach(idx=1:nrow(dfCat))%do%{ #for (idx in foundSources){
for (idx in seq(1,nrow(dfCat))){
  check=try({ # Catch any failed runs
   
  # The row index for this source
  sourceName = format(dfCat$UID[idx], scientific=FALSE)
  printf("\n\n##################################################################################")
  printf("## INFO: Running structural decomposition on {sourceName} ({idx}/{nrow(dfCat)}) ##")
  printf("##################################################################################")
  
  # Create output plot file if it does not already exist
  if (!file.exists(glue("{plotDir}/{sourceName}"))){
    printf("INFO: Output directory did not exist, creating now.")
    dir.create(glue("{plotDir}/{sourceName}"))
  }
  
  numExps = min(c(dfCat$num_exps[idx], maxExps))
  cat(glue("INFO: {dfCat$num_exps[idx]} exposures found (using {numExps}): \n\n"))
  if (numExps == 0){
    printf("INFO: No frames found for source: {sourceName}")
    write(sourceName,file=badIDFilename,append=TRUE)
    next
  }
  
  ### Load each exposure's image & DQ map from the observation set
  expNameList = list() # The names of each exposure
  magzpList = list() # The magnitude zeropoints
  
  imgCutList = list()
  hdrCutList = list()
  dqCutList = list()
  maskList = list()
  psfList = list()
  
  imgWarpList = list()
  maskWarpList = list()
  
  skyList = list()
  skyRMSList = list()
  
  for (jj in seq(1,numExps)){
    expName = dfCat[[glue('name_exp{jj}')]][[idx]]
    expNameList[[jj]] = expName
    obsetID = paste0( substr(toupper(expName), 1, 6), '010' )
    
    # Establish which proposal ID this exposure came from.
    if (substr(expName,1,3) == 'j8p'){proposalID = 'PID09822'}
    else if (substr(expName,1,3) == 'j8x') {proposalID = 'PID10092'}
    else if (substr(expName,1,3) == 'jco') {proposalID = 'PID13641'}
    else if (substr(expName,1,3) == 'jbo') {proposalID = 'PID12440'}
    else if (substr(expName,1,3) == 'jbh') {proposalID = 'PID12328'}
    else {printf("WARNING: No proposal ID could be found for exposure {expName}.")}
    
    # Get source position info in exposure
    expChip = dfCat[[glue('chip_exp{jj}')]][idx]
    expx = dfCat[[glue('x_exp{jj}')]][idx]
    expy = dfCat[[glue('y_exp{jj}')]][idx]
    srcLoc = c(expx, expy)
    
    # The FITS extensions indices for the chip
    sciIndex = ifelse(expChip == 1, 2+3, 2) # The index in HST fits files that points to the science image HDU
    errIndex = ifelse(expChip == 1, 3+3, 3) # The index in HST fits files that points to the error HDU
    dqIndex = ifelse(expChip == 1, 4+3, 4) # The index in HST fits files that points to the data quality HDU
    
    # Load the image, data quality (dq) map and hdr objects
    imgFilename = glue("{framesDir}/{proposalID}/{obsetID}/{expName}_flt.fits")
    printf("INFO: Loading exposure #{jj}: {imgFilename}")
    hdulist = Rfits_read_all(imgFilename, pointer=FALSE)
    
    img = hdulist[[sciIndex]]$imDat
    dq = profoundMakeSegimDilate(segim=hdulist[[dqIndex]]$imDat,size=3,expand=c(4096,8192))$segim # Expand the cosmic rays by a pixel either side (i.e. size=3)
    hdr = hdulist[[sciIndex]]$hdr
    
    # Calculate magnitude zeropoint
    fluxlambda = from_header(hdr = hdr, keys = 'PHOTFLAM', as=as.numeric)
    pivotlambda = from_header(hdr = hdr, keys = 'PHOTPLAM', as=as.numeric)
    magzpList[[expName]] = -2.5*log10(fluxlambda) - 5*log10(pivotlambda) - 2.408 # instrumental zeropoint magnitudes in AB mags (from https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints)
    
    # Now get cutouts and WCS-adjusted headers
    ### TODO: change this to Rwcs_cutout. Does that exist?
    printf("INFO: Source location in {expName}: {srcLoc[1]}, {srcLoc[2]}")
    wcsCutout = magcutoutWCS(img, header=hdr, loc=srcLoc, loc.type='image', box=cutoutBox, plot=FALSE)
    imgCutList[[expName]] = wcsCutout$image
    hdrCutList[[expName]] = wcsCutout$header
    
    dqCutList[[expName]] = magcutout(dq, loc=srcLoc, box=cutoutBox, plot=FALSE)$image
    maskList[[expName]] = apply(apply(dqCutList[[expName]], c(1,2), mask_pixel), c(1,2), as.integer) # Convert the binary dq map to a list of 1s and 0s
    #maskList[[expName]][is.na(imgCutList[[expName]])] = 1 # Also set missing values (NaNs) to be masked
    
    # Calculate the initial sky statistics of each exposure, so that they can be inverse-variance stacked later
    skyGrid = profoundMakeSkyGrid(image=imgCutList[[expName]], mask=maskList[[expName]], box=skyBoxDims)
    
    # Warp images and dqs onto a common WCS.
    printf("INFO: Warping image and dq map to a common WCS")
    imgWarpList[[expName]] = Rwcs(image_in=imgCutList[[jj]], header_in=hdrCutList[[jj]], header_out=hdrCutList[[1]])$image
    maskWarpList[[expName]] = apply(round(magwarp(image_in=maskList[[jj]], header_in=hdrCutList[[jj]], header_out=hdrCutList[[1]], doscale=FALSE)$image, digits=0), c(1,2), as.integer)
    skyList[[expName]] = magwarp(image_in=skyGrid$sky, header_in=hdrCutList[[jj]], header_out=hdrCutList[[1]])$image
    skyRMSList[[expName]] = magwarp(image_in=skyGrid$skyRMS, header_in=hdrCutList[[jj]], header_out=hdrCutList[[1]])$image
    
    
    # Load or create PSF
    psfFilename = glue('{psfsDir}/{sourceName}_{expName}_psf.fits')
    if (!file.exists(psfFilename) | clobberPSF==TRUE){
      printf("INFO: Creating new PSF, saving to file '{psfFilename}'.")
      
      #>> Rscript Generate_HST_PSFs.R [dir] [name] [chip] [x] [y] [focus=0.0] [size=10.0] [filter=f814w] [spectrum=13]
      print(glue("INFO: Running TinyTim with:\n>> Rscript {psfScriptPath} {psfsDir} {sourceName}_{expName} {expChip} {expx} {expy} 0.0 {sizePSF}"))
      invisible(
        system(glue("Rscript {psfScriptPath} {psfsDir} {sourceName}_{expName} {expChip} {expx} {expy} 0.0 {sizePSF}"))
      )
    }
    
    psf = Rfits_read_image(psfFilename,ext=1)$imDat
    
    if (trimPSF == TRUE){ # Trim the PSF?
      printf("INFO: Trimming PSF")
      dimsPSF = sizePSF/pixScaleHST
      peakPSF = which(psf==max(psf,na.rm=T) , arr.ind = T)
      psfList[[expName]] = psf[(peakPSF[1]-dimsPSF%/%2):(peakPSF[1]+dimsPSF%/%2), (peakPSF[2]-dimsPSF%/%2):(peakPSF[2]+dimsPSF%/%2)]
    } else {
      psfList[[expName]] = psf
    }
    
    printf("INFO: Dimensions of {expName} PSF: {dim(psfList[[expName]])[1]} x {dim(psfList[[expName]])[2]}\n")
    
  }
  
  print("INFO: Exposure names:")
  print(names(imgCutList))
  printf("INFO: Length of images: {length(imgCutList)}")
  printf("INFO: Length of dq: {length(dqCutList)}")
  printf("INFO: Length of mask: {length(maskList)}")
  printf("INFO: Length of PSF: {length(psfList)}")
  
  
  # Stack the images and masks
  printf("INFO: Creating median stack of {length(imgCutList)} cutout images.")
  imgCutStack = profoundMakeStack(image_list=imgWarpList, sky_list=skyList, skyRMS_list=skyRMSList, mask_list=maskWarpList, magzero_in = unlist(magzpList), masking='&')$image
  #imgCutStack = profoundMakeStack(image_list=imgCutList, sky_list=skyList, skyRMS_list=skyRMSList, mask_list=maskList, magzero_in = unlist(magzpList))$image # unwarped images
  
  maskStack = maskList[[1]]
  if (numExps > 1){
    for (ii in seq(2,numExps)){
      maskStack = maskStack & maskList[[ii]]
    }
  }
  maskStack = apply(maskStack, c(1,2), as.integer)
  
  ### Run profoundProfound on the median stack of cutouts:
  printf("INFO: Creating segmentation map for median stacked image.")
  seg = profoundProFound(image=imgCutStack, mask=maskStack, sigma=1.5, dilete=5, skycut=1.0, SBdilate=2.5, tol=5, reltol=0.25, ext=5, plot=T)#toPlot)
  
  if (toPlot|toSave){
    if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-segmentation_map.png"), width=650, height=650, pointsize=20)}
    profoundSegimPlot(image=imgCutStack, segim=seg$segim, mask=maskStack)
    if (toSave){dev.off()}
  }
  
  ### Subtract sky
  printf("INFO: Measuring sky statistics for each exposure and subtracting from cutout images.")
  segimList = list()
  for (jj in seq(1,numExps)){
    expName = dfCat[[glue('name_exp{jj}')]][idx]
    
    segimList[[expName]] = apply(magwarp(image_in=seg$segim, header_in=hdrCutList[[1]], header_out=hdrCutList[[jj]], doscale=FALSE, interpolation='nearest')$image, c(1,2), as.integer)
    
    if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-{expNameList[jj]}_sky_statistics.png"), width=720, height=720, pointsize=16)}
    if (toPlot|toSave){par(mfrow=c(2,2), mar = c(1,1,1,1))}
    skyList[[expName]] = profoundProFound(image=imgCutList[[jj]], segim=segimList[[expName]], mask=maskList[[jj]],magzero = magZP, gain=ccdGain, plot=T,
                                          sigma=1.5, skycut=0.75, SBdilate=3, size=31, box=c(150,150))#skyBoxDims)
    
    if (toPlot|toSave){
      plot(skyList[[expName]])
    }
    
    imgCutList[[jj]] = imgCutList[[jj]] - skyList[[jj]]$sky
    if(toSave){dev.off()}
  }
  
  ## Save images of the image cutouts + overlaid segmentation and mask maps
  if (toPlot|toSave){
    if (numExps <= 4){
      if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-cutouts.png"), width=900, height=900, pointsize=10)}
      par(mfrow=c(2,2), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
    } else {
      if (toSave){png(filename = glue("{plotDir}/{sourceName}/{sourceName}-cutouts.png"), width=1800, height=900, pointsize=10)}
      par(mfrow=c(2,4), mar = c(1.5,1.5,1.5,1.5), oma = c(0.25,0.25,0.25,0.25))
    }
    for (jj in seq(1,numExps)){
      profoundSegimPlot(image=imgCutList[[jj]], segim=segimList[[jj]], mask=maskList[[jj]], bad=NA)
      text(0.1*cutoutBox[1],0.925*cutoutBox[2],expNameList[jj], adj=0, col='red', cex=1.6)
    }
    if (toSave){dev.off()}
  }
  
  ### Set up Found2Fit objects for each exposure
  printf("INFO: Generating found2Fits objects for each exposure cutout.")
  
  offsets = list() # The offsets are just the sub-pixel offsets between positions as the cutouts are already centered on the pixel positions.
  rotations = list()
  
  cds = from_header(hdrCutList[[1]], keys = c('CD1_1','CD1_2','CD2_1','CD2_2'), as=as.numeric)
  xrot1 = rad2deg(atan2(cds[3], cds[1]))
  for (jj in seq(1,numExps)){ 
    cds = from_header(hdrCutList[[jj]], keys = c('CD1_1','CD1_2','CD2_1','CD2_2'), as=as.numeric)
    rotations[[jj]] = rad2deg(atan2(cds[3], cds[1]))
    
    offsets[[jj]] = c(dfCat[[glue('x_exp{jj}')]][idx]%%1 - dfCat[['x_exp1']][idx]%%1, dfCat[[glue('y_exp{jj}')]][idx]%%1 - dfCat[['y_exp1']][idx]%%1, xrot1 - rotations[[jj]])
  }
  
  dataList = profuseMultiImageFound2Fit(image_list = imgCutList, psf_list = psfList, segim_list = segimList, mask_list = maskList, offset_list = offsets,
                                          SBdilate=2.5, reltol=1.5, ext=5,
                                          magzero=magZP, gain = rep(ccdGain,numExps),
                                          Ncomp=modelOpts$n_comps,
                                          sing_nser_fit=modelOpts$sing_nser_fit,
                                          disk_nser_fit=modelOpts$disk_nser_fit,  disk_nser=modelOpts$disk_nser,
                                          bulge_nser_fit=modelOpts$bulge_nser_fit, bulge_nser=modelOpts$bulge_nser, bulge_circ=modelOpts$bulge_circ,
                                          pos_delta = 10, # does this pos_delta value work for the biggest galaxies?
                                          rough=roughFit, tightcrop = TRUE, deblend_extra=FALSE, fit_extra=FALSE, plot=toPlot)
  
  
  # Can I actually do this? - yes!
  for (jj in seq(1, numExps)){
    names(dataList)[names(dataList) == glue("image{jj}")] = expNameList[jj]
  }
  
  printf("INFO: Running Highlander optimisation of model.")
  
  printf("Optim iters: {optimOpts$optim_iters}\nNiters: {optimOpts$Niters}\nNfinalMCMC: {optimOpts$NfinalMCMC}\n")
  # Change to profuseMultiImageDoFit(, ...)
  highFit = profuseMultiImageDoFit(image_list = imgCutList, dataList, ablim=1,
                                   optim_iters=optimOpts$optim_iters, Niters=optimOpts$Niters, NfinalMCMC=optimOpts$NfinalMCMC,
                                   CMAargs = optimOpts$CMA, LDargs = optimOpts$LD)
  
  if (toPlot|toSave){
    if (toSave) {png(filename=glue('{plotDir}/{sourceName}/{sourceName}_single_multifit.png'), width=(numExps+1)*300, height=1800, pointsize=20)}
    multiImageMakePlots(datalist=dataList, imagestack=imgCutStack, segim=seg$segim, parm=highFit$parm,
                        plottext=glue("UID:\n{sourceName}\n\n # Exposures: {dataList$Nim}\n log(Mstar): {format(round(log10(dfCat[['StellarMass']][idx]), 2), nsmall = 2)}\n z: {format(round(dfCat[['zBest']][ii], 3), nsmall = 2)}"))
    
    if (toSave){dev.off()}
  }
  
  profitLikeModel(highFit$parm,Data = dataList$j8pu3cijq, makeplots = T)
  
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
  outputs = list(args=args,
                 idx=idx,
                 sourceName=sourceName,
                 modelName=modelName,
                 modelOptions=modelOpts,
                 optimOptions=optimOpts,
                 cutoutBox=cutoutBox,
                 skyBoxDims=skyBoxDims,
                 dqList=dqCutList,
                 expNames=expNameList,
                 imgList=imgCutList,
                 imgStack=imgCutStack,
                 skyList=skyList,
                 hdrList=hdrCutList,
                 maskList=maskList,
                 psfList=psfList,
                 segimList=segimList,
                 offsets=offsets,
                 dataList=dataList,
                 highFit=highFit,
                 optimModellist=optimModellist)
  
  dir.create(file.path(resultsDir, sourceName), showWarnings = FALSE) # Create output directory if one does not already exist
  outputFilename = glue('{resultsDir}/{sourceName}/{sourceName}_{modelName}_output.rds')
  saveRDS(object=outputs, file=outputFilename)
  
  }) # end try
  if (class(check)=='try-error'){
    printf("ERROR: source {sourceName} ({idx}) failed to be fit.")
    write(sourceName,file=badIDFilename,append=TRUE)
  } else {
    printf("INFO: source {sourceName} ({idx}) successfully fit.")
  }
  
} # end main galaxy for loop




### Run in RStudio afterwards ###
#result = readRDS('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/ProFuse_Outputs/101505979732574752/101505979732574752_deV-exp_output.rds')

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
## The old way of running Highlander
# highFit = Highlander(found2Fits[[1]]$Data$init, Data=dataList, likefunc=profitLikeModel, ablim=1, 
#                      optim_iters=optimOpts$optim_iters, Niters=optimOpts$Niters, NfinalMCMC=optimOpts$NfinalMCMC,
#                      CMAargs = optimOpts$CMA, LDargs = optimOpts$LD)

## The old way of initiating the found2Fits (without using profuseMultiImageFound2Fits)
# found2Fits0 = list()
# for (jj in seq(1,numExps)){
#   
#   expName = dfCat[[glue('name_exp{jj}')]][idx]
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
#     text(0.1*cutoutBox[1],0.925*cutoutBox[2],as.character(dfCat[[glue('name_exp{jj}')]][idx]), adj=0, col='white', cex=1.75)
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


## The oldest way of initiating the found2Fits
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

