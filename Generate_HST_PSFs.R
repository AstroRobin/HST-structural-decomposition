### Generate_HST_PSFs.R ###
# A script for generating HST ACS Point Spread Functions (PSFs) using the PSF modelling software Tiny Tim.

# Author: R. H. W. Cook
# Date: 25/01/2022

# To-Do: 
#   * Add get_focus() function to obtain focus value from historical measurements.

library(glue)

#.libPaths(c(libPath,.libPaths()))

evalGlobal = TRUE
Sys.setenv(TINYTIM='/Users/00092380/Documents/Software/tinytim') # Need to change for Pawsey machine.

TinyTimACSR = function(input_dir, name, chip, x, y, filter, spectrum=13, size=3.0, focus=0, exBmV=0, jitter=0){
  tinyTimPath = '/Users/00092380/Documents/Software/tinytim'
  setwd(input_dir)
  
  psf_runfile <- paste(input_dir,'/',name,'_run',sep='')
  fileConn <- file(psf_runfile)
  writeLines(paste(tinyTimPath,"/tiny1 ",name,"_parameter_file ", "ebmv=",exBmV," << EOF",sep=''), fileConn)
  close(fileConn)
  
  cat(15 ,'\n',file=psf_runfile,append=TRUE)      # camera: 15 -> ACS - Wide Field Channel 
  cat(chip ,'\n',file=psf_runfile,append=TRUE)    # CCDCHIP         
  cat(paste0(x,' ',y) ,'\n',file=psf_runfile,append=TRUE)       # x/y-position      
  cat(filter ,'\n',file=psf_runfile,append=TRUE)  # filter
  cat(1 ,'\n',file=psf_runfile,append=TRUE)       # form of object spectrum: 1 -> from list
  cat(spectrum ,'\n',file=psf_runfile,append=TRUE)      # 13   G8V      1.16   0.75   0.52   0.94   1.40   1.84  ???
  cat(size ,'\n',file=psf_runfile,append=TRUE)     # the PSF diameter; recomended 3.0 arcesecond
  #cat(focus ,'\n',file=psf_runfile,append=TRUE)   # the focus (secondary mirror despace)
  cat(paste(name,'_image',sep='') ,'\n',file=psf_runfile,append=TRUE)
  cat("EOF" ,'\n',file=psf_runfile,append=TRUE)
  
  
  system(paste('source ',input_dir,'/',name,'_run',sep=''))
  Sys.sleep(1)
  system(paste(tinyTimPath,'/tiny2 ',input_dir,'/',name,'_parameter_file',sep=''))
  Sys.sleep(3)
  system(paste(tinyTimPath,'/tiny3 ',input_dir,'/',name,'_parameter_file',sep=''))
  
  system(paste('rm ', input_dir, '/', name, '_image00_psf.fits',sep= '')) 
  system(paste('rm ', input_dir, '/*_parameter_*',sep= ''))
  system(paste('rm ', input_dir, '/*.tt3',sep= ''))
  system(paste('rm ', input_dir, '/*_run',sep= ''))
  system(paste('mv ', name,'_image00.fits ', name,'_psf.fits',sep= ''))
  
  #setwd(wrk_dir)
}

### Load in any command-line arguments ###
# Rscript Generate_HST_PSFs.R [dir] [name] [chip] [x] [y] [focus=0.0] [size=10.0] [filter=F814W] [spectrum=13] 
args = commandArgs(trailingOnly = TRUE) # Parse arguments (if given)
#args = c('/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs','test1','2','1000','1000','NA','NA')
#args = c('/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/DEVILS_Structural_Decomposition/PSFs/Test_PSFs/Focus','testfocus1','2','1000','1000','auto','NA')

if (length(args) < 5) { # get command line args
  stop(glue("Require, at least, the arguments: [name] [chip] [x] [y], but only {length(args)} were given."))
} else {
  
  dir = args[1] # directory for PSFs
  if (!file.exists(dir)){
    dir.create(dir)
  }
  
  name = args[2] # name of PSF outputs
  
  chip = as.integer(args[3]) # CCD chip number
  if (is.na(chip) | chip > 2 | chip < 1){ stop(glue("Chip must be able to be coerced into integer of 1 or 2, but instead received {args[3]}.")) }
  
  x = as.integer(args[4]) # x position
  y = as.integer(args[5]) # y position
  if (is.na(x) | is.na(y)){ stop(glue("x/y must be able to be coerced into integer, but instead received {args[4]}, {args[5]}.")) }
  
  if (is.na(args[6]) | args[6]=='NA') {  # focus
    focus = 0.0
  } else if (args[6] == 'auto'){
    focus = 0.0#get_focus()
  } else {
    focus = as.numeric(args[6])
    if (is.na(focus)) {stop(glue("focus must either be 'auto' or able to be coerced into float, but instead received {args[6]}."))}
  }
  
  size = ifelse(is.na(args[7]) | args[7]=='NA', 10, as.numeric(args[7])) # Output PSF size
  if (is.na(size)) {stop(glue("size must be able to be coerced into float, but instead received {args[7]}."))}
  
  filter = ifelse(is.na(args[8]) | args[8]=='NA', 'f814w', args[8]) # Filter
  
  spectrum = ifelse(is.na(args[9]) | args[9]=='NA', 13, as.integer(args[9])) # Spectrum choice
  if (is.na(spectrum) | spectrum > 17){ stop(glue("spectrum must be able to be coerced into integer and not greater than 17, but instead received {args[9]}.")) }
  
}


TinyTimACSR(input_dir=dir, name=name, chip=chip, x=x, y=y, filter=filter, spectrum=spectrum, size=size, focus=focus)

