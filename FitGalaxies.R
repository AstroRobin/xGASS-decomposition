###
# Script to cycle through galaxy fitting for a list of specified galaxies in ugriz frequency bands for 1 and 2 components.
#
# Author: Robin Cook
# Date: 05/05/17
###

### TO-DO: 
# Test SkyMap measurement

# Get home directory
# "HOME" specifies the machine being used:
# if (HOME == "/home/robincook"): HP-Laptop
# if (HOME == "/home/rcook"): ICRAR ~ Munro Cluster Computer
# if (HOME == "/Users/robincook"): Macbook
HOME = paste("/",unlist(strsplit(getwd(), '/'))[2],"/",unlist(strsplit(getwd(), '/'))[3],sep="")

# If running on Munro (ICRAR): must specify library path:
if (HOME == "/home/rcook") {.libPaths(c("/home/arobotham/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))}

library(ProFit) # Bayesian galaxy fitting tool
library(ProFound) # Source Extraction and Image Segmentation tool
library(EBImage) # Image processing package
library(magicaxis) # "Magically Pretty Plots" ~ ASGR
library(FITSio) # .FITS file input/output
library(LaplacesDemon) # MCMC optimisation package
if (HOME != '/home/rcook'){
  library(rgl) # 3D real-time rendering system for R
  library(rpanel) # Set of cuntions to build simple GUI controls for R functions
}

###################################################################
######################### DEFINE CONSTANTS ########################
###################################################################

### The base directory for the galaxy data
# Directory Structure:
# > GALS_DIR
#   > [galName]
#     > [band]
#       > [galName]_[band].fits          \
#       > [galName]_[band]_PSF.fits      / *inputs*
#     > Fitting
#       > [form]
#         > [galName]-[form]_[band]_[nComps]comp_{OutputType}  > *outputs*

#GALS_DIR = paste(HOME,"/Documents/PhD/GASS/Galaxies",sep = "")
GALS_DIR = paste(HOME,"/home/robincook/Google Drive/PhD/Miscellaneous/Lange2016_Galaxies/Galaxies",sep = "")
PIXSCALE = 0.396 # The pixel scale of the image (here: SDSS)

###################################################################

###################################################################
######################### DEFINE FUNCTIONS ########################
###################################################################

plotSegIm = function(image, segim, mask, sky = 0, ...) # A modified version of the profoundPlotSegim() function to adjust colours and linewidths
{
  # <param: image [array (float, 2)]> - The image matrix
  # <param: segim [array (int, 2)]> - The segmentation map matrix create with profoundProfound()
  # <param: mask [array (int, 2)]> - The mask matrix.
  # <param: sky [float]> - The value of the sky.
  
  # <return: NULL>
  
  image = image - sky
  temp = magimage(image, ...)
  if (min(segim, na.rm = TRUE) != 0) {
    segim = segim - min(segim, na.rm = TRUE)
  }
  segvec = which(tabulate(segim) > 0)
  for (i in segvec) {
    z = segim == i
    z = z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x, temp$y, z, add = T, col = rainbow(1000,start=0.49,end=0.51)[sample(1000, 1)], zlim = c(0, 1), lwd=2.5, drawlabels = FALSE, nlevels = 1)
  }
  if (!missing(mask)) {
    magimage(mask, lo = 0, hi = 1, col = c(NA, hsv(alpha = 0.3)), 
             add = T)
  }
}


divide_magnitude = function(magTot,frac=0.5) # Function to divide a magnitude by some fraction
{
  # <param: magTot [float]> - The inital total magnitude of the source
  # <param: frac [float]> - The fraction by which to divide the magnitude.
  # <return: mags [array (float, 2)]> - The divided magnitudes for the bulge and disk component.
  
  fluxTot = 10^(-0.4*magTot)
  mag1 = -2.5*log10(frac*fluxTot)
  mag2 = -2.5*log10((1.0-frac)*fluxTot)
  
  return(c(mag1,mag2))
}

find_main = function(sourceList,dims) # function to find the main (central) source ID in the image
{
  # <param: sourceList [list]> - The list of sources found in profitMakeSegIm().
  # <param: dims [array (float, 2)]> - The x and y dimensions of the image.
  # <return: mainID [int]> - The ID for the main (centremost) source.
  
  nSources = length(sourceList$segID) # get number of sources.
  
  x0 = dims[1]/2; y0 = dims[2]/2 # image centre position.
  sepArr = array(0,dim=nSources) # empty separation array
  
  # Calculate separation of source centres.
  for (ii in seq(1,nSources)) {sepArr[ii] = sqrt( (sourceList$xcen[ii] - x0)^2 + (sourceList$ycen[ii] - y0)^2 )}
  mainID = which.min(sepArr) # The main source is the one with the smallest separation from the centre
  return(mainID)
}

write_output = function(file,name,nComps,init,optim){ # Write optimisation result to file
  # <param: file [str]> - The file to append the results to.
  # <param: name [str]> - The name of the galaxy.
  # <param: nComps [int (1|2)]> - The number of components for the model (1 or 2).
  # <param: init [list]> - A list of initial values for the fit.
  # <param: optim [list]> - A list of optimised values.
  
  # <return: NULL>
  
  # Determine the number of columns in the table.
  nCols = 2*nComps*8 + 1 # 2 (inital/optimised) * nComps (1|2) * 8 (num. parameters in Sersic Profile) + 1 (Galaxy ID)
  
  # Write initial and optimised parameters ro file
  if (nComps == 1){
    write(c("GASSID","x_in","y_in","mag_in","re_in","nser_in","ang_in","axrat_in","box_in","x_out","y_out","mag_out","re_out","nser_out","ang_out","axrat_out","box_out"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name,
            init$xcen[1],init$ycen[1],init$mag[1],init$re[1],init$nser[1],init$ang[1],init$axrat[1],init$box[1],
            optim$xcen[1],optim$ycen[1],optim$mag[1],optim$re[1],optim$nser[1],optim$ang[1],optim$axrat[1],optim$box[1]),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  } else if (nComps == 2){
    write(c("GASSID","x1_in","x2_in","y1_in","y2_in","mag1_in","mag2_in","re1_in","re2_in","nser1_in","nser2_in","ang1_in","ang2_in","axrat1_in","axrat2_in","box1_in","box2_in","x1_out","x2_out","y1_out","y2_out","mag1_out","mag2_out","re1_out","re2_out","nser1_out","nser2_out","ang1_out","ang2_out,axrat1_out,axrat2_out,box1_out,box2_out"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name, 
            init$xcen[1],init$xcen[2], init$ycen[1],init$ycen[2], init$mag[1],init$mag[2], init$re[1],init$re[2], init$nser[1],init$nser[2], init$ang[1],init$ang[2], init$axrat[1],init$axrat[2], init$box[1],init$box[2], 
            optim$xcen[1],optim$xcen[2], optim$ycen[1],optim$ycen[2], optim$mag[1],optim$mag[2], optim$re[1],optim$re[2], optim$nser[1],optim$nser[2], optim$ang[1],optim$ang[2], optim$axrat[1],optim$axrat[2], optim$box[1],optim$box[2]),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  }
}


add_pseudo_bulge = function(model) # Add a zero-point magnitude bulge to the model
{
  # <param: model [list]> - The modellist from profitMakeModel()
  # <return: pseudo [list]> - A (pseudo) modellist containing a second component; i.e., a duplicate with "zero magnitude".
  
  pseudo = model # copy model list
  for (key in names(model$sersic)){pseudo$sersic[[key]][2] = model$sersic[[key]][1]} # duplicate a second component
  pseudo$sersic$mag[1] = ZERO_POINT # set magnitude of bulge (loc=1) to ZERO_POINT
  return(pseudo)
}


###################################################################

###################################################################
################### DEFINE FITTING PARAMETERS #####################
###################################################################

### Fitting Mode ("optim":BFGS , "LA":LaplacesApproximation , "LD":Full-MCMC)
mode="LD"

### Specify the number of components to be fit
#compList = c(1,2)
compList = c(2)

### Specify frequency bands
#bandList = c('u','g','r','i','z')
bandList = c('r')

### Specifying the Optimisation run metedata:
# prefix: This specifies the 'form' of this particular optimisation run and will be used in the filenames of output files (e.g. Test, NoInits, Alpha, etc.)
# flavour: The purpose of this optimisation run:
#       - testing: for testing new features
#       - example: for running an example optimisation
#       - fullsample: running the full sample of galaxies
# description: A verbose description of the optimisation run.
prefix = "Alpha"
flavour = "fullsample"
description = "Alpha: First full-sample optimisation run of xGASS galaxies.\n This uses:\n - simple first-pass segmentating routine\n - matrix sky-subtraction\n - MCMC (CHARM) optimisation w/ 1e4 samples\n - initial condition improvement via isophotal fitting/Laplaces Approximation optimisation\n\nThis particular run contains only galaxies for which there are no complicatons with segmentation maps, the galaxy is too difficult to model, or the galaxy is confirmed throgh visual inspection to be fit well with a single component."

### Specify which galaxies to fit. (Requires image and PSF files.)
args = commandArgs(trailingOnly = TRUE) # Parse arguments
n = as.integer(args[1]) # n is the line number within galFile for which to get the galaxy list

galFile = paste(HOME,"/Documents/PhD/Fitting/Samples/Alpha/Alpha_ServerList.txt",sep="") # The path to the file containing the lines of galaxy lists for each core on the
lines = readLines(galFile)
galList = strsplit(lines[n],'[,]')[[1]] # Get the list of galaxies

# Specify whether to output images or not
output = TRUE

outputInputs = TRUE # Inputs: image, segmentation, sigma map, PSF
outputIsophotes = TRUE # If isophotes are fit; displays the isophotes and 1D fit
outputSegStats = TRUE # Table of segmentation objects with their parameters
outputSkyStats = TRUE # Plot of the Sky noise statistics
outputInitial = TRUE # Likelihood + Ellipse
outputOptimised = TRUE # Likelihood + Ellipse
outputCorner = TRUE # MCMC posteriors
outputModel = TRUE # Initial+optimised model images
outputResults = TRUE # A single line result of the initial/optimised values
outputSummary = TRUE # Summary of the entire fit
outputMCMCResults = TRUE # Output results from MCMC; includes values + statistics.
outputMCMCSummary = TRUE # Summary of the MCMC fit.
outputWorkspace = TRUE # R workspace; only turn off if trying to conserve space.

# Specify verbosity
verb = TRUE

##################################################################

###################################################################
########################## OPTIMISATION ###########################
###################################################################

count = 1 # A running count of galaxies
for (galName in galList){ # loop through galaxies
  for (band in bandList){ # loop through bands
    for (nComps in compList){ # loop through number of components.
      if(verb){cat(paste("\n* ",galName," * [band = ",band,"; comps = ",nComps,"]"," (",count,"/",length(galList),")\n",sep=""))}
      ### INPUTS ### -> otherwise looped
      # galName = "GASS111029"
      # band = "r"
      # nComps = 2
      
      ### Get image file ###
      if(verb){cat("INFO: Retrieving data.\n")}
      imgFilename = paste(galName,"_",band,".fits",sep="")
      imgFile = paste(GALS_DIR,galName,band,imgFilename,sep='/')
      image0 = readFITS(imgFile)$imDat # image0 is the non sky-subtracted image
      header = readFITS(imgFile)$hdr
      dims = dim(image0)
      # Check for NaN padding:
      padded = is.element(NaN,image0)
      
      ### Get PSF file ###
      psfFilename = paste(galName,"_",band,"_PSF.fits",sep="")
      psfFile = paste(GALS_DIR,galName,band,psfFilename,sep='/')
      psf = readFITS(psfFile)$imDat
      
      ### Create outputs folder ###
      if(verb){cat("INFO: Creating output directories.\n")}
      dir.create(paste(GALS_DIR,galName,"Fitting",sep='/'), showWarnings = FALSE) # Suppress warning if directory already exists.
      # Check Prefix validity:
      if (prefix == '' || is.null(prefix)){print("WARNING: Setting prefix to 'Default'"); prefix = 'Default'}
      outputDir = paste(GALS_DIR,galName,"Fitting",prefix,sep='/')
      dir.create(outputDir, showWarnings = FALSE)  # Suppress warning if directory already exists.
      baseFilename = paste(galName,"-",prefix,"_",band,"_",nComps,"comp",sep="")
      	

      ### Get information from FITS header ###
      # Referencing keywords in header
      # <VALUE> = as.numeric(header[which(header=="<KEYWORD>")+1])
      SOFT_BIAS = as.numeric(header[which(header=="SOFTBIAS")+1])
      ZERO_POINT = as.numeric(header[which(header=="ZP")+1])
      GAIN = as.numeric(header[which(header=="GAIN")+1])
      
      ### Subtract SOFT_BIAS from image and PSF ###
      image0 = image0 - SOFT_BIAS
      psf = psf - SOFT_BIAS
      
      
      ############################################################
      ###### Measure sky statistics with profoundProfound() ######
      ############################################################
      if(verb){cat("INFO: Creating Segmentation image.\n")}
      
      # Run profund to get sky statistics
      segmentation = profoundProFound(image0, sigma=2.5, skycut=1.5, tolerance=3, size=15,
                                    magzero=ZERO_POINT, gain=GAIN, header=header,
                                    stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      
      # Extract sky measurements from segmentation output
      sky = segmentation$sky
      skyRMS = segmentation$skyRMS
      skyStats = capture.output(maghist(sky,plot=FALSE)) # A string containing the outputs from calculating statistics on the sky values
      
      ### Plot input images ###
      if (output && outputInputs){
        statsFilename = paste(galName,"_",band,"_SkyStats.png",sep='')
        png(paste(outputDir,statsFilename,sep='/'),width=900,height=300,pointsize=16)
        par(mfrow=c(1,3), mar=c(4,1,2,1))
        
        magimage(sky,axes=FALSE); text(0.1*dims[1],0.925*dims[2],"Sky",adj=0,col='white',cex=1.75)
        magimage(skyRMS,axes=FALSE); text(0.1*dims[1],0.925*dims[2],"Sky RMS",adj=0,col='white',cex=1.75)
        maghist(sky,grid=TRUE,density=25,col='red',xlab="Sky counts",verbose=FALSE)
        
        dev.off()
      }
      
      # Subtract the background sky
      image = image0 - sky
      
      
      ###########################################################
      #####  Make Segmentation map with ProFit (/ProFound)  #####
      ###########################################################
      if(verb){cat("INFO: Creating Segmentation image.\n")}
      
      # Extract sources
      segmentation = profoundProFound(image, sigma=2.0, skycut=1.0, tolerance=2, size=15,
                                      magzero=ZERO_POINT, gain=GAIN, header=header,
                                      stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      
      # Find the main (central) source
      if(verb){cat("INFO: Finding central source.\n")}
      mainID = find_main(segmentation$segstats,dims) # The main source is the one with the smallest separation from the centre
      
      # Expand Segmentation image
      segmentationDilated = profitMakeSegimDilate(image, segmentation$segim, size=25, expand=mainID,
                                                  magzero=ZERO_POINT, gain=GAIN, header=header,
                                                  stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      
      # @Hosein Hashemi:
      #segmentation = profitProFound(image, sigma=4, skycut=2, tolerance=5, size=11, pixcut = 5,
      #                              magzero=ZERO_POINT, gain=GAIN, header=header,
      #                              stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
    
      
      #Create segmentation image from only the central source
      segMap = segmentationDilated$segim
      segMap[segMap!=mainID]=0 # only use the central source
      
      # [visualisation] Make a segmentation map image using pixels from 'image'
      segMapIm = image
      segMapIm[segMap!=mainID]=0
      
      # Save to file
      if(output && outputSegStats){
        segStatsFilename = paste(galName,"_",band,"_SegmentationStats.csv",sep='')
        write.csv(segmentation$segstats,file=paste(outputDir,segStatsFilename,sep='/'),quote=FALSE,row.names=FALSE)
      }
      
      
      ##########################################
      #####   Make Sigma map with ProFit   #####
      ##########################################
      if(verb){cat("INFO: Making sigma map.\n")}
      sigma = profitMakeSigma(image,sky=0.0,skyRMS=skyRMS,gain=GAIN) # sky level is defined as 0.0 as sky has already been subtracted.
      
      
      ### Plot input images ###
      if (output && outputInputs){
        inputsFilename = paste(galName,"_",band,"_Inputs.png",sep='')
        png(paste(outputDir,inputsFilename,sep='/'),width=500,height=500,pointsize=16)
        par(mfrow=c(2,2), mar=c(0.4,0.4,1,1))
        magimage(image,axes=F,bad=0); text(0.2*dim(image)[1],0.925*dim(image)[2],"Image",col='white',cex= 1.75)
        plotSegIm(image,segim=segMap,axes=F,lwd=3,bad=0); text(0.4*dim(image)[1],0.925*dim(image)[2],"Segmentation",col='white',cex= 1.75) ## Test without foreach ***
        magimage(sigma,axes=F,bad=0); text(0.2*dim(image)[1],0.925*dim(image)[2],"Sigma",col='white',cex= 1.75)
        magimage(psf,axes=F,bad=0); text(0.2*dim(psf)[1],0.925*dim(psf)[2],"PSF",col='white',cex= 1.75)
        dev.off()
      }
      
      
      ####################################
      #####  Get Initial Conditions  #####
      ####################################
      if(verb){cat("INFO: Getting initial conditions.\n")}
      improveInits = TRUE # Whether to try isophotal fitting or quick LaplacesApproximation() to improve initial conditions
      
      # Rough Initial conditions from segmentation objects
      inits = segmentation$segstats
      if (nComps == 2){
        magInits = divide_magnitude(inits$mag[mainID],frac=0.4) # Arbitrary 40%/60% division of flux to Bulge/Disk
        reInits = c(inits$semimaj[mainID]*1.0,inits$semimaj[mainID]*3)
        nSerInits = c(4,1)
      } else {
        magInits = c(inits$mag[mainID])
        reInits = c(inits$semimaj[mainID]*1.0)
        nSerInits = c(4)
      }
        
      ## Attempt to improve initial conditions via Isophote fitting:
      # Only run if running 2 components and the image does not contain NaN padding (i.e. galaxies on frame edges).
      if (improveInits == TRUE && nComps == 2 && padded == FALSE){
        if(verb){cat("INFO: Attempting to improve initial conditions.\n")}
        if (output && outputIsophotes){
          isophotesFilename = paste(galName,"_",band,"_Isophotes.png",sep='')
          png(paste(outputDir,isophotesFilename,sep='/'),width=1050,height=500,pointsize=16)
          par(mfrow=c(1,2), mar=c(3.5,3.5,1,2))
        }
        
        # Get the ellipse isophotes
        ellipses = profitGetEllipses(image,segim=segmentation$segim,segID=mainID,levels=20,pixscale=PIXSCALE,magzero=ZERO_POINT,dobox=FALSE,plot=output)
        rMin = 1; rMax = 30; rDiff = 0.1
        rLocs=seq(rMin,rMax,by=rDiff)
        rCut = (7.5-rMin)/rDiff+1
        
        bulge=profitRadialSersic(rLocs, mag=magInits[1], re=reInits[1], nser=nSerInits[1])
        disk=profitRadialSersic(rLocs, mag=magInits[2], re=reInits[2], nser=nSerInits[2])
        
        # Define minimisation function
        sumsq1D=function(par=c(magInits[1], log10(reInits[1]), log10(nSerInits[1]), magInits[2], log10(reInits[2]), log10(nSerInits[2])),rad, SB, pixscale=PIXSCALE)
        {
          bulge=profitRadialSersic(rad, mag=par[1], re=10^par[2], nser=10^par[3])
          disk=profitRadialSersic(rad, mag=par[4], re=10^par[5], nser=10^par[6])
          total=profitFlux2SB(bulge+disk, pixscale=pixscale)
          return=sum((total-SB)^2)
        }
        
        ## Interval bounds:
        # Mag (B/D): 10 - 25
        # Re (B/D): 0.25 - 100
        # nSer (B): 1.25 - 10; nSer (D): 0.5 - 1.5
        lower=c(10,-0.6021,0.09691,10,-0.6021,-0.30103)
        upper=c(25,2,1,25,2,0.1761)
        
        # Run L-BFGS-B optimisation:
        if(verb){cat("INFO: Running BFGS 1D isophotal optimisation.\n")}
        isoFit=optim(sumsq1D, par=c(magInits[1], log10(reInits[1]), log10(nSerInits[1]), magInits[2], log10(reInits[2]), log10(nSerInits[2])),
                    rad=ellipses$ellipses$radhi[4:19], SB=ellipses$ellipses$SB[4:19], pixscale=PIXSCALE,
                    method='L-BFGS-B', lower=lower, upper=upper)$par
        
        bulge=profitRadialSersic(rLocs, mag=isoFit[1], re=10^isoFit[2], nser=10^isoFit[3])
        disk=profitRadialSersic(rLocs, mag=isoFit[4], re=10^isoFit[5], nser=10^isoFit[6])
        
        if (output && outputIsophotes) {
          magplot(ellipses$ellipses$radhi[4:19], ellipses$ellipses$SB[4:19],ylim=c(25,17), grid=TRUE, type='l',lwd=2.5,
                xlab='Radius (pixels)', ylab='Surface Brightness (mag/arcsec^2)',main="Bulge+Disk Surface Brightness Profile")
          points(ellipses$ellipses$radhi[4:19],ellipses$ellipses$SB[4:19],pch=16)
          lines(rLocs, profitFlux2SB(bulge, pixscale=PIXSCALE), col='red',lty='dashed',lwd=2)
          lines(rLocs, profitFlux2SB(disk, pixscale=PIXSCALE), col='blue',lty='dashed',lwd=2)
          lines(rLocs, profitFlux2SB(bulge+disk, pixscale=PIXSCALE), col='green',lwd=2)
          dev.off()
        }
        
        # Check for divergence
        isoConverge = TRUE # IF isoConverge=TRUE THEN use these as initial conditions. ELSE run acesApproximation()
        convParams = list("mag1"=TRUE,"re1"=TRUE,"nser1"=TRUE,"mag2"=TRUE,"re2"=TRUE,"nser2"=TRUE)
        
        if(verb){cat("INFO: Isophotal 1D fitting results:\n")}
        ii = 1
        for (param in isoFit){
          if (is.element(param,lower) || is.element(param,upper)){
            if (ii==1){convParams["mag1"]=FALSE}
            if (ii==2){convParams["re1"]=FALSE}
            if (ii==3){convParams["nser1"]=FALSE}
            if (ii==4){convParams["mag2"]=FALSE}
            if (ii==5){convParams["re2"]=FALSE}
            if (ii==6){convParams["nser2"]=FALSE}
          }
          if(verb){
            if (ii == 1 || ii == 4){
              cat(paste(names(convParams)[ii],': ',toString(param),'\n',sep=''))
            } else {
              cat(paste(names(convParams)[ii],': ',toString(10^param),'\n',sep=''))
            }
          }
          ii = ii + 1
        }
        
         
        if ((convParams['mag1']==T && convParams['mag2']==T) && (convParams['re1']==T && convParams['re2']==T)){
          isoConverge = TRUE
        } else {
          isoConverge = FALSE
        }
        
        ### Check whether bulge dominates flux in outer regions.
        for (ii in seq(rCut,(rMax-rMin)/rDiff+1,1)){
          if (bulge[ii] > disk[ii]){
             isoConverge = FALSE
             if(verb){cat(paste("WARNING: Bulge dominates beyond r > ",toString(rCut*PIXSCALE),"\" (r = ",toString((ii*rDiff+rCut)*PIXSCALE),"\").\n",sep=''))}
             break
          }
        }
        
        # If solution has converged, set the isophotal fit solution as the initial conditions
        if (isoConverge==TRUE){
          if(verb){cat("INFO: Replacing initial conditions with isophotal 1D fit solutions.\n")}
          magInits = c(isoFit[1],isoFit[4])
          reInits = c(10^isoFit[2],10^isoFit[5])
          nSerInits = c(10^isoFit[3],10^isoFit[6])
        } else {
          if(verb){cat("WARNING: isoFit solution did not converge.\n")}
        }
        
      } else { # IF nComps = 1 THEN move onto LAFit
        isoConverge = FALSE
      } # END iosphotal 1D fitting optimisation
      
      
      ##################################
      #####  DEFINE PROFIT INPUTS  #####
      ##################################
      if(verb){cat("INFO: Defining ProFit inputs.\n")}
      if (nComps == 1){
        ## SINGLE SERSIC PROFILE ##
        modellist = list(
          sersic=list(
            xcen= c(inits$xcen[mainID]),
            ycen= c(inits$ycen[mainID]),
            mag=  magInits,
            re=   reInits,
            nser= nSerInits,
            ang=  c(inits$ang[mainID]),
            axrat=c(inits$axrat[mainID]), # Bulge is initially at axrat=1
            box=c(0) # no boxiness
          )
        )
        
        # Parameters to fit
        tofit=list(
          sersic=list(
            xcen= TRUE,
            ycen= TRUE, 
            mag=  TRUE, 
            re=   TRUE, 
            nser= TRUE, #  nser(disk) _=_ free.
            ang=  TRUE, 
            axrat=TRUE, 
            box=  FALSE # do not fit for boxiness
          )
        )
        
        # Define the parameters which should be fitted in log space
        tolog=list(
          sersic=list(
            xcen= FALSE,
            ycen= FALSE,
            mag=  FALSE,
            re=   TRUE, # re is best fit in log space
            nser= TRUE, # nser is best fit in log space
            ang=  FALSE, 
            axrat=TRUE, # axrat is best fit in log space
            box=  FALSE
          )
        )
        
        # Define the sigma values for priors object
        sigmaArr = c(2,2,5,1,1,30,0.3,Inf)
        priors=function(new,init,sigmas=sigmaArr){
          return(sum(dnorm(log10(new$sersic$re),log10(init$sersic$re),sigmas[4],log=TRUE)))
        }
        
        sigmas = list(
          sersic = list(
            xcen= c(sigmaArr[1]),
            ycen= c(sigmaArr[2]),
            mag= c(sigmaArr[3]),
            re= c(sigmaArr[4]),
            nser= c(sigmaArr[5]),
            ang= c(sigmaArr[6]),
            axrat= c(sigmaArr[7]),
            box= c(sigmaArr[8])
          )
        )
        
        priors = profitMakePriors(modellist=modellist, sigmas=sigmas, tolog=tolog, tofit=tofit, allowflat = TRUE) # allowflat allows for flat priors (e.g. boxiness) where the log-likelihood will be computed as 0, rather than -Inf.
        
        # The hard intervals should also be specified in linear space.
        intervals=list(
          sersic=list(
            xcen=list(lim=c(inits$xcen[mainID]-10,inits$xcen[mainID]+10)),
            ycen=list(lim=c(inits$ycen[mainID]-10,inits$ycen[mainID]+10)),
            mag=list(lim=c(7,ZERO_POINT)),
            re=list(lim=c(0.25,100)),
            nser=list(lim=c(0.25,10)),
            ang=list(lim=c(-180,360)),
            axrat=list(lim=c(0.05,1)),
            box=list(lim=c(-0.5,0.5))
          )
        )
        
        # END Single Sersic 
      } else if (nComps == 2){
        ## DOUBLE SERSIC PROFILE ##
        modellist = list(
          sersic=list(
            xcen= c(inits$xcen[mainID],inits$xcen[mainID]),
            ycen= c(inits$ycen[mainID],inits$ycen[mainID]),
            mag = magInits,
            re=   reInits,
            nser= nSerInits,
            ang=  c(inits$ang[mainID],inits$ang[mainID]),
            axrat=c(1,inits$axrat[mainID]), # Bulge is initially at axrat=1
            box=c(0,0) # no boxiness
          )
        )
        
        # Parameters to fit
        tofit=list(
          sersic=list(
            xcen= c(TRUE,NA), # fit both sersics for xcen together *** BOUND
            ycen= c(TRUE,NA), # fit both sersics for xcen together *** BOUND
            mag=  c(TRUE,TRUE), # fit for both magnitudes
            re=   c(TRUE,TRUE), # fit for both effective radii
            nser= c(TRUE,TRUE), # fit for bulge, nser(disk) _=_ free.
            ang=  c(FALSE,TRUE), # fit for disk only
            axrat=c(FALSE,TRUE), # *********** fit for disk, bulge has axrat _=_ 1
            box=  c(FALSE,FALSE)  #
          )
        )
        
        # Define the parameters which should be fitted in log space
        tolog=list(
          sersic=list(
            xcen= c(FALSE,FALSE),
            ycen= c(FALSE,FALSE),
            mag=  c(FALSE,FALSE),
            re=   c(TRUE,TRUE), # re is best fit in log space
            nser= c(TRUE,TRUE), # nser is best fit in log space
            ang=  c(FALSE,FALSE), 
            axrat=c(TRUE,TRUE), # axrat is best fit in log space
            box=  c(FALSE,FALSE)
          )
        )
        
        
        # Define the sigmas object
        sigmaArr = c(2,2,5,1,1,30,0.3,Inf)
        priors=function(new,init,sigmas=sigmaArr){
          return(sum(2*dnorm(log10(new$sersic$re),log10(init$sersic$re),sigmas[4],log=TRUE)))
        }
        
        sigmas = list(
          sersic = list(
            xcen= c(sigmaArr[1],sigmaArr[1]),
            ycen= c(sigmaArr[2],sigmaArr[2]),
            mag= c(sigmaArr[3],sigmaArr[3]),
            re= c(sigmaArr[4],sigmaArr[4]),
            nser= c(sigmaArr[5],sigmaArr[5]),
            ang= c(sigmaArr[6],sigmaArr[6]),
            axrat= c(sigmaArr[7],sigmaArr[7]),
            box= c(sigmaArr[8],sigmaArr[8])
          )
        )
        
        priors = profitMakePriors(modellist, sigmas, tolog, allowflat = TRUE) # allowflat allows for flat priors (e.g. boxiness) where the log-likelihood will be computed as 0, rather than -Inf.
        
        # The hard intervals should also be specified in linear space.
        intervals=list(
          sersic=list(
            xcen=list(lim=c(inits$xcen[mainID]-10,inits$xcen[mainID]+10),lim=c(inits$xcen[mainID]-10,inits$xcen[mainID]+10)),
            ycen=list(lim=c(inits$ycen[mainID]-10,inits$ycen[mainID]+10),lim=c(inits$ycen[mainID]-10,inits$ycen[mainID]+10)),
            mag=list(lim=c(10,25),lim=c(10,25)),
            re=list(lim=c(0.25,100),lim=c(0.25,100)),
            nser=list(lim=c(1.25,10),lim=c(0.5,1.5)),
            ang=list(lim=c(-180,360),lim=c(-180,360)),
            axrat=list(lim=c(0.1,1),lim=c(0.1,1)),
            box=list(lim=c(-0.5,0.5),lim=c(-0.5,0.5)) 
          )
        )
      
      } # END Double Sersic
    
      
      ### Setup Data ###
      if(verb){cat("INFO: Setting up Data object.\n")}
      Data = profitSetupData(image=image,sigma=sigma,modellist=modellist,tofit=tofit,tolog=tolog,priors=priors,intervals=intervals,
                             psf=psf, magzero=ZERO_POINT,segim=segMap, algo.func='optim',like.func="t",verbose=FALSE)
      
      
      ### Plot Input Model Likelihood ###
      if (output  && outputInitial){
        initLikelihoodFilename = paste(baseFilename,"_LikelihoodInitial.png",sep='')
        png(paste(outputDir,initLikelihoodFilename,sep='/'),width=1600,height=490,pointsize=28)
        profitLikeModel(parm=Data$init,Data=Data,makeplots=TRUE,plotchisq=FALSE)
        dev.off()
      }
      
      ### Plot Input Model Ellipse ###
      if(output  && outputInitial){
        initEllipseFilename = paste(baseFilename,"_EllipseInitial.png",sep='')
        png(paste(outputDir,initEllipseFilename,sep='/'),width=1000,height=750,pointsize=20)
        if (nComps == 1){
          try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modellist),pixscale=PIXSCALE,SBlim=25))
        } else if (nComps == 2){
          try(profitEllipsePlot(Data=Data,modellist=modellist,pixscale=PIXSCALE,SBlim=25))
        }
        dev.off()
      }
      
      
      #########################################
      #####   Optimise Model (+ timing)   #####
      #########################################
      
      # IF the result from 1D isophotal fitting did not converge THEN attempt a aceApproximation() fit:
      if (improveInits==TRUE && isoConverge==FALSE) { # LaplaceApproximation LM fit
        if(verb){cat("INFO: Attempting to improve inital conditions with LaplaceApproximation()\n")}
        Data$algo.func = "LA" # Change optimising algorithm
        LAFit = LaplaceApproximation(profitLikeModel,parm = Data$init, Data = Data, Iterations=1e3,
                                     Method = 'LM', CovEst='Identity', sir = FALSE)
        
        # Remake model
        optimModel=profitRemakeModellist(LAFit$Summary1[,1],Data$modellist,Data$tofit,Data$tolog)$modellist
        
        ### Check for divergence
        if(verb){cat("INFO: Testing LAFit for divergence.\n")}
        LAConverge = TRUE # IF converged=TRUE THEN use these as initial conditions. ELSE run LaplacesApproximation()
        for (n in seq(nComps)){
          if (is.element(optimModel$sersic$mag[n],Data$intervals$sersic$mag[[n]])){LAConverge=FALSE}
          if (is.element(optimModel$sersic$re[n],Data$intervals$sersic$re[[n]])){LAConverge=FALSE}
        }
        
        if (nComps == 2){ # Create 1D profiles for bulge/disk
          rMin = 1; rMax = 30; rDiff = 0.1
          rLocs=seq(rMin,rMax,by=rDiff)
          rCut = (7.5-rMin)/rDiff+1
          
          bulge=profitRadialSersic(rLocs, mag=optimModel$sersic$mag[1], re=optimModel$sersic$re[1], optimModel$sersic$nser[1])
          disk=profitRadialSersic(rLocs, mag=optimModel$sersic$mag[2], re=optimModel$sersic$re[2], optimModel$sersic$nser[2])
        
          ### Check whether bulge dominates flux in outer regions.
          for (ii in seq(rCut,(rMax-rMin)/rDiff+1,1)){
            if (bulge[ii] > disk[ii]){
              LAConverge = FALSE
              if(verb){cat(paste("WARNING: Bulge dominates beyond r > ",toString(rCut*PIXSCALE),"\" (r = ",toString((ii*rDiff+rCut)*PIXSCALE),"\").\n",sep=''))}
              #break
            }
          }
        }
        
        # If solution has converged, set the LAFit solution as the initial conditions
        if (LAConverge==TRUE){
          if(verb){cat("INFO: Replacing initial conditions with LAFit solutions.\n")}
          Data$init = LAFit$Summary1[,1]
        } else {
          if(verb){cat("INFO: LAFit solution did not converge.\n")}
        }
        
      }  # END LAFit optimisation
      
      if(verb){cat(paste("\n\n* ",galName," * [band = ",band,"; comps = ",nComps,"]"," (",count,"/",length(galList),")\n\n",sep=""))}
      
      ### Laplaces Demon (Full MCMC)
      if (mode == "LD"){
        if(verb){cat("INFO: Starting LaplacesDemon() ~ Full MCMC optimisation.\n")}
        startTime = Sys.time()
        
        Data$algo.func = "LD"
        LDFit = LaplacesDemon(profitLikeModel, Initial.Values = Data$init, Data=Data, Iterations=1e4, Algorithm='CHARM',Thinning=1,Specs=list(alpha.star=0.44), Status=2500)
        
        #bestLD=magtri(LDFit$Posterior2,samples=500,samptype='end')
        if(output && outputCorner){
          # Determine the index of the first parameter for plotting
          if (nComps == 1){
            par0 = 3 # 
          } else {
            if (is.na(tofit$sersic$xcen[2])){
              par0 = 3 # xcen2/ycen2 are fixed, therefore not in fitting parameters
            } else{
              par0 = 5 # xcen2/ycen2 are free, therefore are included in fitting parameters
            }
          }
          
          cornerFilename = paste(baseFilename,"_CornerPlot.png",sep='')
          png(paste(outputDir,cornerFilename,sep='/'),width=2400,height=2000,pointsize=28)
          capture.output( # Use capture.output to suppress LA outputs
            magtri(LDFit$Posterior1[,par0:NCOL(LDFit$Posterior1)],samples=1000,samptype='end',grid=TRUE,tick=FALSE),
            file = '/dev/null' # Redirect output to /dev/null
          )
          
          dev.off()
        }
        
        bestLD=magtri(LDFit$Posterior1,samples=1000,samptype='end')
        dev.off()
        
        optimModel = profitRemakeModellist(bestLD[,1],Data$modellist,Data$tofit,Data$tolog)$modellist
        
        ### Plot Optimised Model Likelihood ###
        if(output  && outputOptimised){
          optimLikelihoodFilename = paste(baseFilename,"_LikelihoodOptimised.png",sep='')
          png(paste(outputDir,optimLikelihoodFilename,sep='/'),width=1600,height=490,pointsize=28)
          profitLikeModel(bestLD[,1],Data,makeplots=TRUE,plotchisq=FALSE)
          dev.off()
        }
        
        endTime = Sys.time()
        elapsedTime =  difftime(endTime,startTime,units="secs")
        
        ### Save MCMC summary to file
        if (output && outputMCMCSummary){
          MCMCSummaryFilename = paste(baseFilename,"_MCMCSummary.txt",sep='')
          sink(file=paste(outputDir,MCMCSummaryFilename,sep='/'),append=FALSE)
          print(LDFit)
          sink()
        }
        
        ### Save MCMC output parameters + statistics to file
        if (output && outputMCMCResults){
          MCMCOutputFilename = paste(baseFilename,"_MCMCOutput.csv",sep='')
          write.csv(LDFit$Summary1, file=paste(outputDir,MCMCOutputFilename,sep='/'),quote=FALSE,row.names=TRUE)
        }
      }
      
      
      ### Remake intial/optimised model ###
      initImage = profitMakeModel(modellist,dim=dim(image),magzero=ZERO_POINT)$z
      optimImage = profitMakeModel(optimModel,dim=dim(image),magzero=ZERO_POINT)$z
      
      ### Plot Model Images ###
      if(output && outputModel){
        noise = matrix( rnorm(dims[1]*dims[2],mean=0.0,sd=skyRMS), dims[1], dims[2]) 
        initModelFilename = paste(baseFilename,"_ModelImage.png",sep='')
        png(paste(outputDir,initModelFilename,sep='/'),width=1200,height=400,pointsize = 20)
        par(mfrow=c(1,3), mar=c(0.4,0.4,1,1))
        magimage(image,axes=F,bad=0); text(0.15*dim(image)[1],0.95*dim(image)[2],"Image",col='white',cex= 2)
        magimage(initImage+noise,axes=F,bad=0); 
        text(0.25*dim(image)[1],0.95*dim(image)[2],"Initial Model",col='white',cex= 2)
        text(0.715*dim(image)[1],0.05*dim(image)[2],expression(paste("N: ",mu, "=0, ", sigma,"=")),col='white',cex= 1.5)
        text(0.9*dim(image)[1],0.05*dim(image)[2],format(skyRMS,digits=3),col='white',cex= 1.5)
        magimage(optimImage+noise,axes=F,bad=0);
        text(0.3*dim(image)[1],0.95*dim(image)[2],"Optimised Model",col='white',cex= 2)
        text(0.715*dim(image)[1],0.05*dim(image)[2],expression(paste("N: ",mu, "=0, ", sigma,"=")),col='white',cex= 1.5)
        text(0.9*dim(image)[1],0.05*dim(image)[2],format(skyRMS,digits=3),col='white',cex= 1.5)
        dev.off()
      }
      
      
      ### Plot Optimised Ellipse Plot ###
      if(output && outputOptimised){
        optimEllipseFilename = paste(baseFilename,"_EllipseOptimised.png",sep='')
        png(paste(outputDir,optimEllipseFilename,sep='/'),width=1000,height=750,pointsize=20)
        if (nComps == 1){
          try(profitEllipsePlot(Data=Data,add_pseudo_bulge(optimModel),pixscale=PIXSCALE,SBlim=25))
        }else if (nComps == 2){
          try(profitEllipsePlot(Data=Data,optimModel,pixscale=PIXSCALE,SBlim=25))
        }
        dev.off()
      }
      
      ###########################
      #####  WRITE OUTPUTS  #####
      ###########################
      if(verb){cat("INFO: Writing outputs to file.\n")}
      if(output && outputResults){
        outputFilename = paste(baseFilename,"_Output.csv",sep='')
        write_output(file=paste(outputDir,outputFilename,sep='/'),name=galName,nComps=nComps,init=modellist$sersic,optim=optimModel$sersic)
      }
      
      ###################################################################
      ########################### WRITE SUMMARY #########################
      ###################################################################
      if(output && outputSummary){
        summaryFilename = paste(baseFilename,"_Summary.txt",sep='')
        sink(file=paste(outputDir,summaryFilename,sep='/'),append=FALSE)
        
        cat("##############################################################\n")
        cat("####### ProFit | Galaxy Profile Optimisation - Summary #######\n")
        cat("##############################################################\n\n")
        
        cat(paste("Galaxy: ",galName,"\n",sep=""))
        cat(paste("Band: ",band,"\n",sep=""))
        cat(paste("Num. Components: ",nComps,"\n\n",sep=""))
        cat(paste("Form: ",prefix,"\n",sep=""))
        cat(paste("Flavour: ",flavour,"\n",sep=""))
        cat(paste("\nDescription:   ","\n",description,"\n",sep=""))
        cat(paste("\nDate: ",Sys.time(),"\n",sep=""))
        cat(paste("\nElapsed time: ",elapsedTime,"\n",sep=""))
        
        cat("\n\n>> Inputs:\n")
        cat(paste("Base directory: ",outputDir,"\n",sep=""))
        cat(paste("Image: ",imgFile,"\n",sep=""))
        cat(paste(" Dimensions: ",dims[1]," x ",dims[2]," (",dims[1]*PIXSCALE/60.0,"' x ",dims[1]*PIXSCALE/60.0,"')","\n",sep=""))
        cat(paste(" Padding: ",padded,"\n",sep=""))
        cat(paste("PSF: ",psfFile,"\n",sep=""))
        cat(paste("\nZero point: ",ZERO_POINT,"\n",sep=""))
        cat(paste("Gain: ",GAIN,"\n",sep=""))
        cat(paste("Soft-bias: ",SOFT_BIAS,"\n",sep=""))
        
        cat("\n\n>> Sky Statistics:\n")
        cat(gsub("\\[1\\]","\n",skyStats)) # Print the sky statistics: the output from maghist() of profoundMakeSkyGrid()
        
        cat("\n\n>> Segmentation:\n")
        cat(paste("Num. Objects: ",length(segmentation$segstats[[1]]),"\n",sep=""))
        cat(paste("MainID: ",mainID,"\n",sep=""))
        cat("\n Stats for target object:\n")
        print(segmentation$segstats[mainID,])
        
        cat("\n\n>> Initial Conditions:\n")
        print(modellist$sersic)
        
        cat("\n\n>> Fitting Parameters:\n")
        print(tofit$sersic)
        
        #if ((mode == "LA" || mode == "LD") && nComps == 2){
        #  cat("\n\n>> Constraints:\n")
        #  print(constraints)
        #}
        
        cat("\n\n>> Intervals:\n")
        print(intervals$sersic)
        
        if(improveInits==TRUE){
          cat("\n\n>> Attempted Improvements on Initial Conditions:\n")
          if(isoConverge==TRUE){
            cat("Isophotal 1D fitting successful\n")
            cat("Fitting solutions:\n")
            for (ii in seq(length(convParams))){
              if (ii == 1 || ii == 4){cat(paste(" ",names(convParams)[ii],": ",toString(isoFit[ii]),'\n',sep=''))}
              else {cat(paste(" ",names(convParams)[ii],": ",toString(10^isoFit[ii]),'\n',sep=''))}
            }
          } else {
            cat("Isophotal 1D fitting not succesful\n")
          }
          
          if(isoConverge==FALSE){
            if(LAConverge==TRUE) {cat("\nLaplaceApproximation() 2D fitting successful\n")}
            else {cat("\nLaplaceApproximation() 2D fitting not successful\n")}
            cat("Fitting solutions:\n")
            for (ii in seq(length(LAFit$Summary1[,1]))){
              if(grepl('re',names(LAFit$Summary1[,1][ii])) || grepl('nser',names(LAFit$Summary1[,1][ii])) || grepl('axrat',names(LAFit$Summary1[,1][ii]))){
                cat(paste(" ",names(LAFit$Summary1[,1])[ii],": ",toString(10^LAFit$Summary1[,1][ii]),'\n',sep=''))
              } else {
                cat(paste(" ",names(LAFit$Summary1[,1])[ii],": ",toString(LAFit$Summary1[,1][ii]),'\n',sep=''))
              }
            }
          }
            
        }
        
        cat("\n\n>> Output Model:\n")
        print(optimModel$sersic)
        
        cat("\n\n>> Flags:\n")
        # Check if optimised parameters are stuck to interval bounds
        for (key in names(optimModel$sersic)){
          for (n in seq(1,nComps)){
            if (optimModel$sersic[[key]][n] %in% intervals$sersic[[key]][[n]] && tofit$sersic[[key]][n]){
              cat(paste(" - ",key,n,' = ',optimModel$sersic[[key]][n],'\n',sep=''))
            }
          }
        }
        if(segmentation$segstats[mainID,]$edge_frac < 0.8){cat(paste('\n - Segmentation boundary = ',segmentation$segstats[mainID,]$edge_frac,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$edge_frac < 0.8){cat(paste('\n - Segmentation boundary = ',segmentation$segstats[mainID,]$edge_frac,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$edge_excess > 1.0){cat(paste('\n - Segmentation edge excess = ',segmentation$segstats[mainID,]$edge_excess,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$asymm > 0.2){cat(paste('\n - Segmentation asymmetry = ',segmentation$segstats[mainID,]$asymm,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$flag_border != 0){cat(paste('\n - Segmentation borders with: ',segmentation$segstats[mainID,]$flag_border,'\n',sep=''))}
        
        sink()
      }
      
      ##################################
      ##### Save workspace to file #####
      ##################################
      
      if(output && outputWorkspace){  
        workspaceFilename = paste(baseFilename,"_WorkSpace.RData",sep='')
        save.image(paste(outputDir,workspaceFilename,sep='/'))
      }
      
    
    } # END nComps loop
  } # END band loop
  count = count + 1
} # END Galaxy loop

