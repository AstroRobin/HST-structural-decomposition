### Generate_HST_PSFs.R ###
# A script for generating HST ACS Point Spread Functions (PSFs) using the PSF modelling software Tiny Tim.

# Author: R. H. W. Cook
# Date: 25/01/2022

# To-Do: 
#   * Add [dir] argument

library(glue)

if (libPath != "") {.libPaths(c(libPath,.libPaths()))}
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
# Rscript Generate_HST_PSFs.R [name] [chip] [x] [y] [focus=0.0] [size=10.0] [filter=F814W] [spectrum=13] 
#args = commandArgs(trailingOnly = TRUE) # Parse arguments (if given)
args = c('test1','2','1000','1000','NA','NA')

if (length(args) < 4) { # get command line args
  stop(glue("Require, at least, the arguments: [name] [chip] [x] [y], but only {length(args)} were given."))
} else {
  name = args[1] # name of PSF outputs
  
  chip = as.integer(args[2]) # CCD chip number
  if (is.na(chip) | chip > 2 | chip < 1){ stop(glue("Chip must be able to be coerced into integer of 1 or 2, but instead received {args[2]}.")) }
  
  x = as.integer(args[3]) # x position
  y = as.integer(args[4]) # y position
  if (is.na(x) | is.na(y)){ stop(glue("x/y must be able to be coerced into integer, but instead received {args[3]}, {args[4]}.")) }
  
  if (is.na(args[5]) | args[5]=='NA') {  # focus
    focus = 0.0
  } else if (args[5] == 'auto'){
    focus = get_focus()
  } else {
    focus = as.numeric(args[5])
    if (is.na(focus)) {stop(glue("focus must either be 'auto' or able to be coerced into float, but instead received {args[5]}."))}
  }
  
  size = ifelse(is.na(args[6]) | args[6]=='NA', 10, as.numeric(args[6])) # Output PSF size
  if (is.na(size)) {stop(glue("size must be able to be coerced into float, but instead received {args[6]}."))}
  
  filter = ifelse(is.na(args[7]) | args[7]=='NA', 'f814w', args[7]) # Filter
  
  spectrum = ifelse(is.na(args[8]) | args[8]=='NA', 13, as.integer(args[8])) # Spectrum choice
  if (is.na(spectrum) | spectrum > 20){ stop(glue("spectrum must be able to be coerced into integer, but instead received {args[8]}.")) }
  
}


TinyTimACSR(input_dir=dir, name=name, chip=chip, x=x, y=y, filter=filter, spectrum=spectrum, size=size, focus=focus)

