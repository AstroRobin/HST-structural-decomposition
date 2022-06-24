### Generate_HST_PSFs.R ###
# A script for generating HST ACS Point Spread Functions (PSFs) using the PSF modelling software Tiny Tim.

# Author: R. H. W. Cook
# Date: 25/01/2022

.libPaths(c("/group/pawsey0160/software/sles12sp3/apps/sandybridge/gcc/4.8.5/r/3.6.3/lib64/R/library",.libPaths()))

evalGlobal = TRUE

TinyTimACSR = function(input_dir, name, chip, x, y, filter, spectrum=13, size=5.0, focus=0, exBmV=0, jitter=0, machine='local'){
  
  setwd(input_dir)
  
  if (machine == 'local'){
    tinyTimPath = '/Users/00092380/Documents/Software/tinytim'
  } else if (machine %in% c('magnus', 'zeus', 'setonix')) {
    tinyTimPath = '/group/pawsey0160/rhwcook/Programs/TinyTim/'
  }
  
  psf_runfile = paste(input_dir,'/',name,'_run',sep='')
  fileConn = file(psf_runfile)
  writeLines(paste(tinyTimPath,"/tiny1 ",name,".param ", "ebmv=",exBmV," << EOF",sep=''), fileConn)
  #writeLines(glue('{tinyTimPath}/tiny1 {name}.param ebmv={exBmV} << EOF'), fileConn)
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
  
  invisible(capture.output(system(paste('source ',input_dir,'/',name,'_run',sep=''))))
  Sys.sleep(1)
  
  if (focus != 0.0){ ## Adjust focus if specified
    fileConn = file(paste0(input_dir,'/',name,'.param'), open='r+')
    parLines = readLines(fileConn)
    focusLine = grep("*# Z4 = Focus for center of ACS/WFC field", parLines)
    newFocus = as.numeric(strsplit(parLines[[focusLine]], '#')[[1]][1]) + focus
    parLines[focusLine] = paste0('  ',format(newFocus, digits=4), '   # Z4 = Focus for center of ACS/WFC field')
    writeLines(parLines, fileConn)
    close(fileConn)

    #print(parLines)
    Sys.sleep(1)
  }

  
  invisible(capture.output(system(paste(tinyTimPath,'/tiny2 ',input_dir,'/',name,'.param',sep=''))))
  Sys.sleep(3)
  invisible(capture.output(system(paste(tinyTimPath,'/tiny3 ',input_dir,'/',name,'.param',sep=''))))
  
  Sys.sleep(1)
  system(paste('rm ', input_dir, '/', name, '_image00_psf.fits',sep= '')) 
  system(paste('rm ', input_dir, '/', name, '.param',sep= ''))
  system(paste('rm ', input_dir, '/', name, '_image.tt3',sep= ''))
  system(paste('rm ', input_dir, '/', name, '_run',sep= ''))
  
  Sys.sleep(1)
  print(name)
  system(paste('mv ', name,'_image00.fits ', name,'_psf.fits',sep= ''))
  
  #setwd(wrk_dir)
}



### Load in any command-line arguments ###
# Rscript Generate_HST_PSFs.R --dir=PATH --name=NAME --machine=['local'|'magnus'|'zeus'|'setonix'] --chip=[1|2] --x --y --focus=0.0 --size=5.0 --filter='F814W' --spectrum=13

library(argparse)
parser = ArgumentParser()
parser$add_argument("-d", "--dir", required=TRUE,
                    action='store', dest='dir', type='character', default=NULL,
                    help="The directory in which to run TinyTim.", metavar="PATH")
parser$add_argument("-n", "--name", required=TRUE,
                    action='store', dest='name', type='character', default='test',
                    help="The base name for all TinyTim outputs.", metavar="NAME")
parser$add_argument("-m", "--machine",
                    action='store', dest='machine', type='character', default=NULL, choices=c('local', 'magnus', 'zeus', 'setonix'),
                    help="The machine on which this code is running", metavar="MACHINE")
parser$add_argument("-c", "--chip",
                    action='store', dest='chip', type='integer', default=1, choices=c(1,2),
                    help="HST ACS chip number. Note this is exactly opposite to the ordering of FITS extensions (chip=2 --> 1st ext, chip=1 --> 2nd ext)", metavar="CHIP")
parser$add_argument("-x", required=TRUE,
                    action='store', dest='x', type='integer', default=NULL,
                    help="The x position in the chip to generate the PSF.", metavar="X")
parser$add_argument("-y", required=TRUE,
                    action='store', dest='y', type='integer', default=NULL,
                    help="The y position in the chip to generate the PSF.", metavar="Y")
parser$add_argument("-f", "--focus",
                    action='store', dest='focus', type='double', default=0.0,
                    help="The focus value (in microns) to apply to HST instrument.", metavar="FOCUS")
parser$add_argument("-s", "--size",
                    action='store', dest='size', type='double', default=5.0,
                    help="The size of the resulting PSF image in arcseconds.", metavar="ARCSEC")
parser$add_argument("-F", "--filter",
                    action='store', dest='filter', type='character', default='F814W',
                    help="The HST filter to select, given as a string (e.g. 'F814W').", metavar="FILTER")
parser$add_argument("-S", "--spectrum",
                    action='store', dest='spectrum', type='integer', default=13,
                    help="The selected spectrum number from the list given in tiny1.", metavar="NUMBER")

# Parse arguments and validate
args = parser$parse_args()

# Input working directory for TinyTim
if (!file.exists(args$dir)){
  dir.create(args$dir)
}

if (args$machine == 'local'){
  Sys.setenv(TINYTIM='/Users/00092380/Documents/Software/tinytim')
} else if (args$machine %in% c('magnus', 'zeus', 'setonix')){
  Sys.setenv(TINYTIM='/group/pawsey0160/rhwcook/Programs/TinyTim/') # Pawsey machines.
}

# HST chip x and y positions
if (is.na(args$x) | is.na(args$y)){
  stop(paste0("x/y must be able to be coerced into integer, but instead received ",args$x,", ",args$y,"."))
}

# focus given in micron
if (is.na(args$focus) | args$focus == 'NA') {
  focus = 0.0
} else {
  focus = as.numeric(args$focus) * 0.011 # multiplied by 0.011 to convert from micron to number of wavelengths @ 547 nm
  if (is.na(focus)) {stop(paste0("focus must be able to be coerced into float, but instead received ", args$focus,"."))}
}

# PSF size
if (is.na(args$size)) {stop(past0("size must be able to be coerced into float, but instead received ", args$size,"."))}

print(args)

TinyTimACSR(input_dir=args$dir, name=args$name, chip=args$chip, x=args$x, y=args$y, filter=args$filter, spectrum=args$spectrum, size=args$size, focus=focus, machine=args$machine)

