### Generate_HST_PSF.R ###
# A script for generating HST ACS Point Spread Functions (PSFs) using the PSF modelling software Tiny Tim.

# Author: R. H. W. Cook
# Date: 25/01/2022

library(glue)

### Load Catalogue of DEVILS sources ###
workDir = '/Users/00092380/Documents/GoogleDrive/PostDoc-UWA/Tasks/DEVILS_Structural_Decomposition'

framePositionsFile = glue('{workDir}/Catalogues/DEVILS_D10_HST_ACS_Frame_Positions_PID09822.csv')
dfDEVILS = read.csv(framePositionsFile)

# Location for output PSFs
psfsDir = '/Users/00092380/Documents/Storage/PostDoc-UWA/HST_COSMOS/PSFs'



### Define the DEVILS sample of interest

dfSubset = subset(dfDEVILS, subset = (log10(StellarMass) > 9.5 & zBest < 1.0 & substring(field, 3, 5) == 'D10'))
nrow(dfSubset)

for (ii in seq(10)){#seq(nrow(dfDEVILS))){
  
  print(ii)
  
}
