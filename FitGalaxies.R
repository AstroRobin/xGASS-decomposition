###
# A wrapper around various likelihood optimisation and image segmentation functions in ProFit/ProFound (developed by Aaron Robotham)
# to perform bulge/disk decomposition on SDSS images of galaxies within the xGASS sample. The general flow of the program is:
#
# Specify fitting parameters: galaxy list; SDSS filter; num. components; etc. --> Get data (image + PSF) --> Create output directories
# --> Measure & Subtract sky statistics --> Get segmentation map --> Define model & fitting params --> Optimise Model! --> write results to file.
#
# Author: Robin H. W. Cook
# Date: 05/05/17
###

### TO-DO: 
# Test SkyMap measurement
# Soft-code
#   - Change to config file instead
#   - read galList

cat(paste("",
"         _________   __________             ____                                             _ __  _           \n",
"   _  __/ ____/   | / ___/ ___/            / __ \\___  _________  ____ ___  ____  ____  _____(_) /_(_)___  ____ \n",
"  | |/_/ / __/ /| | \\__ \\\\__ \\   ______   / / / / _ \\/ ___/ __ \\/ __ `__ \\/ __ \\/ __ \\/ ___/ / __/ / __ \\/ __ \\\n",
" _>  </ /_/ / ___ |___/ /__/ /  /_____/  / /_/ /  __/ /__/ /_/ / / / / / / /_/ / /_/ (__  ) / /_/ / /_/ / / / /\n",
"/_/|_|\\____/_/  |_/____/____/           /_____/\\___/\\___/\\____/_/ /_/ /_/ ____/\\____/____/_/\\__/_/\\____/_/ /_/ \n",
"                                                                       /_/                                     \n",sep=""))

###################################################################
######################## LOAD CONFIG FILES ########################
###################################################################

args = commandArgs(trailingOnly = TRUE) # Parse arguments (if given)

configFound = FALSE
if (length(args) != 0) { # get .conf file from command line arguments
  if (file.exists(args[1])){ # check if config file exists
    configFile = args[1]
    configFound = TRUE
  } else {
    cat(paste("Error: config file '",args[1],"' does not exist",sep=""))
  }
} else if ("default.conf" %in% list.files(path=".")){ # Check for default config file "config.R"
  configFile = "default.conf"
  configFound = TRUE
} else { # Look for a '.conf' file in current directory
  for (file in list.files(path=".")){
    if (grepl(".conf",file)){
      configFile = file
      configFound = TRUE
    }
  }
}

if (configFound==TRUE){source(configFile)} # If a .conf file is found, load the configuration values
if (libPath != "") {.libPaths(c(libPath,.libPaths()))}

###################################################################
########################## LOAD PACKAGES ##########################
###################################################################

library(ProFit,warn.conflicts=FALSE) # Bayesian galaxy fitting tool
library(ProFound) # Source Extraction and Image Segmentation tool
library(EBImage,warn.conflicts=FALSE) # Image processing package
library(magicaxis) # "Magically Pretty Plots" ~ ASGR
library(FITSio) # .FITS file input/output
library(LaplacesDemon,warn.conflicts=FALSE) # MCMC optimisation package


###################################################################
######################### DEFINE FUNCTIONS ########################
###################################################################

inc_gamma = function(a,x){ # Calculates the incomplete gamma function \gamma(a,x) = \int_{0}^{x} t^{a-1} e^{-t} dt
  return(pgamma(x, a) * gamma(a))
}


calc_conc = function(n,alpha=1/3){ # Calculates the concentration index for a given Sersic index, n (based upon equations given in Graham & Driver 2005)
  b = qgamma(0.5,2*n)
  return( inc_gamma(2*n, b*alpha^(1/n)) / inc_gamma(2*n,b) )
}


log_seq = function(from=1, to=1, len=50) { # logarithmically spaced sequence
  exp(seq(log(from), log(to), length.out = len))
}


calc_bn = function(n){ # Calculates the b_n integration constant for the Sersic profile.
  # <param: n (float)> - The Sersic index.
  # <return: bn> - The b_n integration constant.
  
  return( qgamma(0.5,2*n) )
}


calc_fn = function(n){ # a function of the Sersic profile
  # <param: n (float)> - The Sersic index.
  
  # <return: fn (float)> - the Sersic function integration constant.
  
  bn = calc_bn(n)
  fn = (n*exp(bn))/(bn^(2.0*n)) * 2.0*inc_gamma(2.0*n,bn)
  
  return(fn)
}


calc_x = function(n,re,r){ # Calculates the x value substitution of the Sersic profile
  # <param: n (float)> - The Sersic index.
  # <param: re (float)> - The effective radius.
  # <param: r (float)> - The radius with which to integrate until.
  
  # <return: x (float)>
  
  return ( calc_bn(n) * (r/re)^(1.0/n) )
}


get_nser = function(c,nMin=0.5,nMax=15.0){ # Given a concentration index, calculate the Sersic index with which
  # <param: c [float]> - The measured concentration index for the galaxy (R50/R90)
  # <param: nMin [float]> - The minimum Sersic index
  # <param: nMax [float]> - The maximum Sersic index
  
  # <return: n [float]> - The corresponding Sersic index for a given concentration index
  
  n_arr = seq(nMin,nMax,by=0.05) # an array of Sersic indices for which to calculate the concentration over
  c_arr = calc_conc(n_arr) # The corresponding concentrations.
  
  # Finding the index at which the concentration array is closest to the provided concentration (i.e. min of absolute differences)
  c_index = which.min(abs(c_arr-c))
  
  # Get the corresponding Sersic index
  n = n_arr[c_index]
  if (n == nMin || n == nMax) {cat("\nWARNING: Calculated Sersic index is beyond interval bounds.\n")}
  
  return(n)
}


get_B2T = function(n,nMax=20.0){ # Given a Sersic index from a single component fit, calculate an estimate for the Bulge/Total ration
  # <param: n [float]> - The single component Sersic index
  # <param: nMax [float]> - The maximum Sersic index (defined such that B2T = 0.9 @ nMax).
  
  # <return: B2T [float]> - The bulge-to-total ratio.
  
  a = 0.2; b = 0.2 # The relation is a logarithmic function that allows a steep rise at low n and a slow plateau beyond n ~ 7
  B2T = 0.9 * (a + b*log(n)) / (a + b*log(nMax)) # Defined to set B2t to 90% at nMax
  return(B2T)
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


calc_B2T = function(mag1, mag2) # Given two magnitudes (bulge and disk), get the ratios of the their.
{
  # <param: mag1 [float]> - The magnitude of the first component ("bulge")
  # <param: mag2 [float]> - The magnitude of the second component (disk)
  
  # <return: B2T [float]> - The bulge-to-total ratio of fluxes.
  
  # calculate fluxes
  flux1 = 10^(-0.4*mag1)
  flux2 = 10^(-0.4*mag2)
  
  # get bulge-to-total ratio
  B2T = flux1/(flux1+flux2)
  
  return(B2T)
}


find_main_index = function(seg) # Determine the main (central) source list index in the image
{
  # <param: seg [list]> - The segmentation object.
  
  # <return: mainIndex [int]> - The index for the main (centremost) source.
  
  nSources = length(seg$segstats$segID) # get number of sources.
  
  dims = dim(seg$segim)
  x0 = dims[1]/2; y0 = dims[2]/2 # image centre position.
  sepArr = array(0,dim=nSources) # empty separation array
  
  # Calculate separation of source centres.
  for (ii in seq(1,nSources)) {sepArr[ii] = sqrt( (seg$segstats$xcen[ii] - x0)^2 + (seg$segstats$ycen[ii] - y0)^2 )}
  mainIndex = which.min(sepArr) # The main source is the one with the smallest separation from the centre
  return(mainIndex)
}


find_main_ID = function(seg) # Determine the main (central) source ID in the image
{
  # <param: seg [list]> - The segmentation object.
  
  # <return: mainID [int]> - The ID for the main (centremost) source.
  
  nSources = length(seg$segstats$segID) # get number of sources.
  
  dims = dim(seg$segim)
  x0 = dims[1]/2; y0 = dims[2]/2 # image centre position.
  sepArr = array(0,dim=nSources) # empty separation array
  
  # Calculate separation of source centres.
  for (ii in seq(1,nSources)) {sepArr[ii] = sqrt( (seg$segstats$xcen[ii] - x0)^2 + (seg$segstats$ycen[ii] - y0)^2 )}
  mainIndex = which.min(sepArr) # The main source is the one with the smallest separation from the centre
  mainID = seg$segstats$segID[mainIndex]
  
  return(mainID)
}


calc_ave_SB = function(mag,re,pixScale=1.0){ # Calculates the average surface brightness within an effective radius.
  # <param: mag [float]> - the total magnitude of the model.
  # <param: re [float]> - the effective radius of the model.
  # <param: pixScale [float]> - The pixel scale (arcsec/pixel).
  
  # <return: aveSB [float]> - the average surface brightness within an effective radius (mag/arcsec^2)
  
  aveSB = mag + 2.5*log10(2*pi*(re*pixScale)^2)
  
  return(aveSB)
}


calc_SB = function(aveSB, n){ # Calculates the surface brightness at the effective radius numerically via the Sersic profile
  # <param: aveSB [float]> - The average surface brightness within an effective radius (mag/arcsec^2).
  # <param: n [float]> - The Sersic index.
  
  # <return: SB [float]> The surface brightness at the effective radius (mag/arcsec^2).
  
  fn = calc_fn(n)
  SB = aveSB + 2.5*log10(fn)
  
  return(SB)
}


calc_SB0 = function(SB, n){ # Calculates the central surface brightness
  # <param: SB [float]> - The surface brightness at the effective radius
  # <param: n [float]> - The Sersic index.
  
  # <return: SB0 [float]> - The central surface brightness.
  
  bn = calc_bn(n)
  SB0 = SB - 2.5*bn/log(10)
  
  return(SB0)
}


truncate_mag = function(SB,re,n,rTrun,pixScale){ # Truncates the total magnitude of a galaxy to some specified radius.
  # <param: SB (float)> - The surface brightness at the effective radius (mag/arcsec^2).
  # <param: re (float)> - The effective radius (pixels).
  # <param: n (float)> - The Sersic index.
  # <param: rTrun (float)> - The truncation radius (pixels).
  # <param: pixScale (float)> - The pixel scale (arcsec/pixel).
  
  # <return: magTrun (float)> - The truncated magnitude (mag/arcsec^2)
  
  bn = calc_bn(n)
  x = calc_x(n, re, rTrun)
  
  mTrun = SB - 5.0*log10(re*pixScale) - 2.5*log10( 2.0*pi*n * (exp(bn))/(bn^(2.0*n)) * inc_gamma(2.0*n,x) )
  
  return(mTrun)
}



write_output = function(file, name, nComps, init, optim, chisq, time, stat){ # Write optimisation result to file
  # <param: file [str]> - The file to append the results to.
  # <param: name [str]> - The name of the galaxy.
  # <param: nComps [int (1|2)]> - The number of components for the model (1 or 2).
  # <param: init [list]> - A list of initial values for the fit.
  # <param: optim [list]> - A list of optimised values.
  # <param: chisq [float]> - The average chi squared value in the fitting region where \chi^2 = ( (image-model) / sigma )[region]
  # <param: time [float]> - The time elapsed for the optimisation stage only (seconds)
  # <param: stat [boolean (1|0)]> - The stationarity (lack of trend) of the optimisation according to the BMK.Diagnostic in Laplaces Demon (1 = stationary; 2 = not stationary)
  
  # <return: NULL>
  
  # Determine the number of columns in the table.
  nCols = 2*nComps*8 + 4 # 2 (inital/optimised) * nComps (1|2) * 8 (num. parameters in Sersic Profile) + 1 (Galaxy ID + chisq + elapsedTime + stationarity)
  
  # Write initial and optimised parameters ro file
  if (nComps == 1){
    write(c("GASSID","x_in","y_in","mag_in","re_in","nser_in","ang_in","axrat_in","box_in","x_out","y_out","mag_out","re_out","nser_out","ang_out","axrat_out","box_out","chisq","elapsed_time","stationarity"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name,
            init$xcen[1],init$ycen[1],init$mag[1],init$re[1],init$nser[1],init$ang[1],init$axrat[1],init$box[1],
            optim$xcen[1],optim$ycen[1],optim$mag[1],optim$re[1],optim$nser[1],optim$ang[1],optim$axrat[1],optim$box[1],
            chisq, time, stat),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  } else if (nComps == 2){
    write(c("GASSID","x1_in","x2_in","y1_in","y2_in","mag1_in","mag2_in","re1_in","re2_in","nser1_in","nser2_in","ang1_in","ang2_in","axrat1_in","axrat2_in","box1_in","box2_in","x1_out","x2_out","y1_out","y2_out","mag1_out","mag2_out","re1_out","re2_out","nser1_out","nser2_out","ang1_out","ang2_out","axrat1_out","axrat2_out","box1_out","box2_out","chisq","elapsed_time","stationarity"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name, 
            init$xcen[1],init$xcen[2], init$ycen[1],init$ycen[2], init$mag[1],init$mag[2], init$re[1],init$re[2], init$nser[1],init$nser[2], init$ang[1],init$ang[2], init$axrat[1],init$axrat[2], init$box[1],init$box[2], 
            optim$xcen[1],optim$xcen[2], optim$ycen[1],optim$ycen[2], optim$mag[1],optim$mag[2], optim$re[1],optim$re[2], optim$nser[1],optim$nser[2], optim$ang[1],optim$ang[2], optim$axrat[1],optim$axrat[2], optim$box[1],optim$box[2],
            chisq, time, stat),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  }
}

add_pseudo_bulge = function(model) # Add a zero-point magnitude bulge to the model
{
  # <param: model [list]> - The modellist from profitMakeModel()
  # <return: pseudo [list]> - A (pseudo) modellist containing a second component; i.e., a duplicate with "zero magnitude".
  
  pseudo = model # copy model list
  for (key in names(model$sersic)){pseudo$sersic[[key]][2] = model$sersic[[key]][1]} # duplicate a second component
  pseudo$sersic$mag[1] = zeroPoint # set magnitude of bulge (loc=1) to zeroPoint
  return(pseudo)
}

get_gal_list = function(galFile, lineNum=0) # Given the path to a file, extract the galaxy list at a particular line(s)
{
  # <param: galFile [string]> - The path to the file containing the galaxy list(s)
  # <param: lineNum [int]> - The line number at which to extract the list(s) (lineNum = 0 for all lines)
  
  # <return: galList [list]> - A list containing the galaxies to be optimised
  
  # Check whether galFile exists
  if (!file.exists(galFile)){
    cat(paste("ERROR: The input galFile: '",galFile,"' does not exist!\n -- ABORTING --\n",sep=""))
    return(c())
  }
  
  # Read the galaxy lists (one list per line)
  lines = readLines(galFile)
  
  # Validate line number reference
  if (max(lineNum) > length(lines)){
    cat("\nWARNING: Line number is greater than the number of lines. Setting lineNum = 1 instead.\n")
    lineNum = c(1)
  } else if (lineNum == 0){
    lineNum = seq(1,length(lines))
  }
  
  # Extract list(s)
  if (length(lineNum) > 1){ # multiple lines specified
    galList = c()
    for (n in lineNum) {galList = c(galList, strsplit(lines[n],'[,]')[[1]])} # concatenate multiple lines into a single list
  } else { # single line specified
    galList = strsplit(lines[lineNum],'[,]')[[1]] # Get the list of galaxies at single line
  }
  
  return(galList)
}

calc_chisq = function(image, modelImage, sigma, segMap) # Calculate the average chi^2 value across the fitting region
{
  # <param: image [float (matrix)]> - The data image matrix
  # <param: modelImage [float (matrix)]> - The optimised model image matrix
  # <param: sigma [float (matrix)]> - The errors in the image
  # <param: segMap [float (matrix)]> - The region where the fitting was performed
  
  # <return: chisq [list]> - The average chi squared value in the fitting region where \chi^2 = ( (image-model) / sigma )[region]
  
  residual = image - modelImage # take the difference between the data image and the optimised model image
  errorMap = residual/sigma
  
  region = (segMap != 0)
  error = errorMap[region] # the map of errors across the image dimensions
  
  chisq = sum(error^2)/sum(region) # average over the optimisation region
  
  return(chisq)
}

##################################################################


###################################################################
#####################  RETRIEVE GALAXY LIST  ######################
###################################################################

if (length(args) > 1) { # Line number has been specified in comman-line argument
  lineNum = as.integer(args[2])
}

### Specifying which galaxies to fit. (Requires image and PSF files.) ###
if (galFile == "") { # galaxy file not given in .conf file; using galList instead.
  if (length(galNames) == 0){
    cat("WARNING: No 'galFile' or 'galList' specified.\n  -- ABORTING --\n")
    galList = c()
  }
  
  galList = galNames
  lineNum = NULL
} else { # galaxy list was given in .conf file
  galList = get_gal_list(galFile, lineNum) # extract galaxy lists from file
}

# Validate whether all galaxies in galList have both image and PSF files
for (band in bandList){
  noImg = c(); noPSF = c() # lists containing galaxies without files
  
  for (galName in galList){
    # The path to the image file
    imgFilename = paste(galName,"_",band,".fits",sep="")
    imgPath = paste(galsDir,galName,band,imgFilename,sep='/')
    
    # The path to the PSF file
    psfFilename = paste(galName,"_",band,"_PSF.fits",sep="")
    psfPath = paste(galsDir,galName,band,psfFilename,sep='/')
    
    if (!file.exists(imgPath)) {noImg = c(noImg, galName)}
    if (!file.exists(psfPath)) {noPSF = c(noPSF, galName)}
    
  }

  # Remove galaxies without image files
  if (length(noImg) > 0){
    cat(paste("WARNING: Galaxies without ",band,"-band image files:\n",sep=""))
    for (galName in noImg){cat(paste(" - ",galName,"\n",sep=""))}
    galList = setdiff(galList, noImg)
  }
  
  # Remove galaxies without PSF files
  if (length(noPSF) > 0){
    cat(paste("WARNING: Galaxies without ",band,"-band PSF files:\n",sep=""))
    for (galName in noPSF){cat(paste(" - ",galName,"\n",sep=""))}
    galList = setdiff(galList, noPSF)
  }
  
  rm(imgPath,imgFilename,psfPath,psfFilename,galName)
}

###################################################################
########################## OPTIMISATION ###########################
###################################################################

count = 1 # A running count of galaxies
for (galName in galList){ # loop through galaxies
  for (band in bandList){ # loop through bands
    for (nComps in compList){ # loop through number of components.
      if(verb){cat(paste("\n* ",galName," * [band = ",band,"; comps = ",nComps,"]"," (",count,"/",length(galList),")\n",sep=""))}
      
      ### INPUTS ### *otherwise taken from .conf file
      # galName = "GASS109005"
      # band = "r"
      # nComps = 2
      
      ### Create outputs folder ###
      if(verb){cat("INFO: Creating output directories.\n")}
      dir.create(paste(galsDir,galName,"Fitting",sep='/'), showWarnings = FALSE) # Suppress warning if directory already exists.
      # Check Prefix validity:
      if (run == '' || is.null(run)){print("WARNING: Setting run name to 'Default'"); run = 'Default'}
      outputDir = paste(galsDir,galName,"Fitting",run,sep='/')
      dir.create(outputDir, showWarnings = FALSE)  # Suppress warning if directory already exists.
      baseFilename = paste(galName,"-",run,"_",band,"_",nComps,"comp",sep="")

      ### Get image file ###
      if(verb){cat("INFO: Retrieving data.\n")}
      imgFilename = paste(galName,"_",band,".fits",sep="")
      imgFile = paste(galsDir,galName,band,imgFilename,sep='/')
      image0 = readFITS(imgFile)$imDat # image0 is the non sky-subtracted image
      header = readFITS(imgFile)$hdr
      
      dims = dim(image0) # image dimensions
      padded = is.element(NaN,image0) # Check for NaN padding in images where the frame edge is present.
      
      
      ### Get information from FITS header ###
      # Referencing keywords in header
      # <VALUE> = as.numeric(header[which(header=="<KEYWORD>")+1])
      zeroPoint = as.numeric(header[which(header=="ZP")+1])
      gain = as.numeric(header[which(header=="GAIN")+1])
      
      
      ### Get PSF file ###
      psfFilename = paste(galName,"_",band,"_PSF.fits",sep="")
      psfFile = paste(galsDir,galName,band,psfFilename,sep='/')
      psf = readFITS(psfFile)$imDat
      
      
      ### IF SDSS: Subtract softBias from image and PSF ###
      if (dataSource == "SDSS"){
        softBias = as.numeric(header[which(header=="SOFTBIAS")+1])
        if (subSoftBias){
          image0 = image0 - softBias
          psf = psf - softBias
        }
      }
      
      ############################################################
      ###### Measure sky statistics with profoundProfound() ######
      ############################################################
      
      # Get sky box car filter dimensions
      if (!is.na(skyBoxNum)){ # Box size specified as num. of boxes within image
        skyBoxDims = rep(as.integer(dims[1]/skyBoxNum),2)
      } else if (!is.na(skyBoxSize)){ # Box size given as absolute
        skyBoxDims = rep(as.integer(skyBoxSize),2)
      } else { # Default box size
        skyBoxDims = rep(101,2)
      }
      
      # Use profound to get sky measurements
      if(verb){cat("INFO: Creating initial sky mask.\n")} # initial sky mask
      skyMap0 = profoundProFound(image0, skycut=1.0, tolerance=5, redosky=TRUE, redoskysize=21, pixcut=9,
                                    #box = c(dims[1]/10,dims[2]/10), grid = c(dims[1]/12,dims[2]/12),type='bicubic',
                                    magzero=zeroPoint, gain=gain, pixscale=pixScale, # header=header,
                                    stats=FALSE, rotstats=FALSE, boundstats=FALSE, plot=FALSE)
      
      if(verb){cat("INFO: Expanding sky mask.\n")} # Expanding the sky mask further with dilation
      skyMapExp = profoundMakeSegimExpand(image0, segim=skyMap0$segim, tolerance=5, sigma=2.5, skycut=-1.0, sky=skyMap0$sky, skyRMS=skyMap0$skyRMS,
                                          magzero=zeroPoint, gain=gain, pixscale=pixScale,# header=header,
                                          stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=FALSE)
      
      skyMapExp = profoundMakeSegimDilate(image0, segim=skyMapExp$segim, tolerance=5, size=15, sky=skyMap0$sky, skyRMS=skyMap0$skyRMS,
                                          magzero=zeroPoint, gain=gain, pixscale=pixScale,# header=header,
                                          stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=FALSE)
      
      if(verb){cat("INFO: Measuring sky map.\n")} # Measuring the sky map across the image
      skyMap = profoundProFound(image0, segim = skyMapExp$segim,
                                 box = skyBoxDims, grid = skyBoxDims,type='bicubic',
                                 magzero=zeroPoint, gain=gain, pixscale=pixScale, # header=header,
                                 stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      
      # Extract sky measurement statistics from image using profitSkyEst()
      if(verb){cat("INFO: Measuring sky statistics.\n")} 
      skyEst = profoundSkyEst(image0, objects = skyMap$objects_redo, plot=FALSE)
      skyVal = skyEst$sky
      skyRMS = skyEst$skyRMS
      
      # Get sky statistics from sky map produced in profoundProFound
      skyStats = capture.output(maghist(skyMap$sky,plot=FALSE))
      
      # Subtract the background sky
      if (subSky){
        if(verb){cat("INFO: Subtracting background sky from image.\n")}
        image = if (skyAsGrid) (image0 - skyMap$sky) else (image0 - skyVal)
      } else {
        image = image0
      }
      
      ### Plot input images ###
      if (output && outputSkyStats){
        skyStatsFilename = paste(baseFilename,"_SkyStats.png",sep='')
        png(paste(outputDir,skyStatsFilename,sep='/'),width=900,height=340,pointsize=16)
        par(mfrow=c(1,3), mar=c(4,2,2,1))
        
        magimage(skyMap$sky,stretch = 'asinh'); text(0.1*dims[1],0.925*dims[2],"Sky",adj=0,col='white',cex=1.75)
        magimage(skyMap$skyRMS,stretch = 'asinh'); text(0.1*dims[1],0.925*dims[2],"Sky RMS",adj=0,col='white',cex=1.75)
        profoundSkyEst(image, objects = skyMap$objects,plot=TRUE, grid=TRUE,xlab="Sky counts")
        
        dev.off()
      }
      
      ###########################################################
      #####  Make Segmentation map with ProFit (/ProFound)  #####
      ###########################################################
      
      if (loadSegMap == TRUE){ # Load the segmentation image from galaxy data directory.
        segMapFilename = paste(galName,"_",band,"_SegMap.fits",sep="") # The name of the segmentation fits file.
        segMapFile = paste(galsDir,galName,band,segMapFilename,sep='/') # The path to the segmentation fits file.
        if (!file.exists(segMapFile)){cat(paste("\nWARNING: \"",segMapFile,"\" does not exist! Making segmentation map from image instead.\n",sep=""))}
      }
      
      if (loadSegMap == TRUE && file.exists(segMapFile)){
        if(verb){cat("INFO: Loading segmentation map.\n")}
        segmentation = list() # Set up a new segmenation[$segim,$objects,$segstats] structure
        
        # Read in the segmentation image
        segmentation$segim = readFITS(segMapFile)$imDat
        
        segmentation$objects = segmentation$segim
        segmentation$objects[segmentation$objects!=0] = 1
        
        # Measure segmentation stats for the segmentation image
        segmentation$segstats =  profoundSegimStats(image, segmentation$segim, magzero=zeroPoint, gain=gain, pixscale=pixScale, rotstats=TRUE, boundstats=TRUE)
        
      } else { # Create a new segmentation image.
      
        if(verb){cat("INFO: Creating Segmentation image.\n")}
        
        # @Robin Cook: Extract sources
        segmentation0 = profoundProFound(image, sigma=segSigma, skycut=segSkyCut, tolerance=segTol, ext=segExt,
                                        magzero=zeroPoint, gain=gain, #header=header,
                                        stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
        
        # Find the main (central) source
        if(verb){cat("INFO: Finding central source.\n")}
        mainID = find_main_ID(segmentation0) # The main source is the one with the smallest separation from the centre
        
        # @Robin Cook: Expand Segmentation image
        if(verb){cat("INFO: Expanding target segment.\n")}
        segmentationExp = profoundMakeSegimExpand(image=image, segim=segmentation0$segim, expand=mainID, skycut=expSkyCut, sigma=expSigma,
                                                    sky = if (subSky) 0.0 else skyEst$sky, skyRMS=skyMap$skyRMS,
                                                    magzero=zeroPoint, gain=gain, #header=header,
                                                    stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
        
        # @Robin Cook: Dilate Segmentation image
        
        if (dilateSize != 0){ # if dilation > 0, perform dilation.
          if(verb){cat("INFO: Performing final dilation.\n")}
          segmentation = profoundMakeSegimDilate(image=image, segim=segmentationExp$segim, expand=mainID, size=dilateSize,
                                                    magzero=zeroPoint, gain=gain, #header=header,
                                                    stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
        } else {
          segmentation = segmentationExp # set final segmentation undilated
        }
        
        # @Hosein Hashemi:
        #segmentation = profitProFound(image, sigma=4, skycut=2, tolerance=5, size=11, pixcut = 5,
        #                              magzero=zeroPoint, gain=gain, header=header,
        #                              stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
        
        # Total extra pixels added to the segmenation map
        numPixExpand = length(segmentation$objects[segmentation$objects != 0]) - length(segmentation0$objects[segmentation0$objects != 0])
      }
      
      # Find the main (central) source
      if(verb){cat("INFO: Finding central source.\n")}
      mainID = find_main_ID(segmentation) # The main source is the one with the smallest separation from the centre
      mainIndex = find_main_index(segmentation)
      
      
      # Save segmentation stats to file
      if(output && outputSegStats){
        segStatsFilename = paste(baseFilename,"_SegmentationStats.csv",sep='')
        write.csv(segmentation$segstats,file=paste(outputDir,segStatsFilename,sep='/'),quote=FALSE,row.names=FALSE)
      }
      
      #Create segmentation image from only the central source
      segMap = segmentation$segim
      segMap[segMap!=mainID]=0 # only use the central source
      
      # [visualisation] Make a segmentation map image using pixels from 'image'
      segMapIm = image
      segMapIm[segMap!=mainID]=0
      
      ##########################################
      #####   Make Sigma map with ProFit   #####
      ##########################################
      
      if (loadSigma == TRUE){ # Load the sigma image from galaxy data directory.
        sigmaFilename = paste(galName,"_",band,"_Sigma.fits",sep="") # The name of the sigma fits file.
        sigmaFile = paste(galsDir,galName,band,sigmaFilename,sep='/') # The path to the sigma fits file.
        if (!file.exists(sigmaFile)){cat(paste("\nWARNING: \"",sigmaFile,"\" does not exist! Making sigma map from image instead.\n",sep=""))}
      }
      
      if (loadSigma == TRUE && file.exists(sigmaFile)){
        if(verb){cat("INFO: Loading sigma map.\n")}
        
        # Read in the sigma image
        sigma = readFITS(sigmaFile)$imDat
      } else {
        if(verb){cat("INFO: Making sigma map.\n")}
        
        # > sky level is defined as 0.0 as sky has already been subtracted.
        # > sigma map reflects the original image, however, sky pixels are set to a fixed uncertainty and object pixels have additional shot noise.
        sigma = profoundMakeSigma(image,sky=if (subSky) 0.0 else skyEst$sky,skyRMS=skyEst$skyRMS,objects=segmentation$objects,gain=gain,plot=FALSE)
      }
      
      
      ######################################
      #####   Plot & Save Input Data   #####
      ######################################
      if (output && outputInputs){
        inputsFilename = paste(baseFilename,"_Inputs.png",sep='')
        png(paste(outputDir,inputsFilename,sep='/'),width=750,height=750,pointsize=16)
        par(mfrow=c(2,2), mar=c(0.4,0.4,1,1))
        magimage(image,axes=F,bad=0); text(0.1*dim(image)[1],0.925*dim(image)[2],"Image",pos=4,col='white',cex= 1.75)
        profoundSegimPlot(image,segim=segmentation$segim,axes=F,lwd=3,bad=0); text(0.1*dim(image)[1],0.925*dim(image)[2],"Segmentation",pos=4,col='white',cex= 1.75) ## Test without foreach ***
        magimage(sigma,axes=F,bad=0); text(0.1*dim(image)[1],0.925*dim(image)[2],"Sigma",pos=4,col='white',cex= 1.75)
        magimage(psf,axes=F,bad=0); text(0.1*dim(psf)[1],0.925*dim(psf)[2],"PSF",pos=4,col='white',cex= 1.75)
        dev.off()
      }
      
      # Save inputs if save[*] is True and not already a loaded input.
      if (savePSF && !loadPSF){ # PSF
        if(verb){cat("INFO: Saving PSF.\n")}
        psfOutFilename = paste(galName,"_",band,"_PSF.fits",sep="")
        writeFITSim(psf, file = paste(galsDir,galName,band,psfOutFilename,sep='/'))
      }
      
      if (saveSegMap && !loadSegMap){ # Segmentation Map
        if(verb){cat("INFO: Saving segmentation map.\n")}
        segOutFilename = paste(galName,"_",band,"_SegMap.fits",sep="")
        writeFITSim(segmentation$segim, file = paste(galsDir,galName,band,segOutFilename,sep='/'))
      }
      
      if (saveSigma && !loadSigma){ # Sigma Map
        if(verb){cat("INFO: Saving sigma map.\n")}
        sigmaOutFilename = paste(galName,"_",band,"_Sigma.fits",sep="")
        writeFITSim(sigma, file = paste(galsDir,galName,band,sigmaOutFilename,sep='/'))
      }
      
      
      #####################################
      #######  Get Initial Guesses  #######
      #####################################
      if(verb){cat("INFO: Getting initial guesses.\n")}
      
      # Rough Initial model from segmentation objects
      inits = segmentation$segstats
      if (nComps == 2){
        xcenInits = rep(inits$xcen[mainIndex],2)
        ycenInits = rep(inits$ycen[mainIndex],2)
        magInits = divide_magnitude(inits$mag[mainIndex],frac=bulgeFrac)
        reInits = c(inits$semimaj[mainIndex]*1/3,inits$semimaj[mainIndex]*1)
        nSerInits = c(if (nFromCon) get_nser(inits$con[mainIndex]) else 4, 1) # Bulge: n = 4 (de Vaucouleurs); Disk: n = 1 (Exponential)
        angInits = rep(inits$ang[mainIndex],2)
        axratInits = c(1,inits$axrat[mainIndex]) # Bulge is initially at axrat=1
        boxInits = rep(0,2)
      } else {
        xcenInits = c(inits$xcen[mainIndex])
        ycenInits = c(inits$ycen[mainIndex])
        magInits = c(inits$mag[mainIndex])
        reInits = c(inits$semimaj[mainIndex]*1.0)
        nSerInits = if (nFromCon) c(get_nser(inits$con[mainIndex])) else c(4)
        angInits = c(inits$ang[mainIndex])
        axratInits = c(inits$axrat[mainIndex])
        boxInits = c(0)
      }
      
      # Create initial modellist
      modellist0 = list(
        sersic=list(
          xcen= xcenInits,
          ycen= ycenInits,
          mag = magInits,
          re=   reInits,
          nser= nSerInits,
          ang=  angInits,
          axrat=axratInits,
          box=  boxInits
        )
      )
      
      # If inheritFrom, replace initial parameters with results from the specified previous optimisation run.
      if (!is.null(inheritFrom)){
        # Get the inheriting RData file
        envirFile = paste(galName,"-",inheritFrom$form,"_",inheritFrom$band,"_",inheritFrom$nComps,"comp_WorkSpace.RData",sep="")
        
        if (file.exists( paste(galsDir,galName,"Fitting",inheritFrom$form,envirFile,sep="/") )){
          tempEnvir = new.env() # Create temporary envronment to place workspace of inheriting optimisation run.
          load(paste(galsDir,galName,"Fitting",inheritFrom$form,envirFile,sep="/"), envir=tempEnvir) # load inheriting RData into temporary environment
          
          inheritModellist = get('modellist',tempEnvir)
          
          if (nComps == inheritFrom$nComps){ # The inheriting parameters have the same model structure
            for (param in names(inheritModellist$sersic)) {
              if (inheritParams$sersic[[param]] == TRUE){
                modellist0$sersic[[param]] = inheritModellist$sersic[[param]]
              }
            }
          } else { # going from 1comp model -> 2comp model
            for (param in names(inheritModellist$sersic)){
              if (inheritParams$sersic[[param]] == TRUE){
                if (param == "xcen" || param == "ycen" || param == "ang" || param == "box"){
                  modellist0$sersic[[param]] = rep(inheritModellist$sersic[[param]],2)
                } else if (param == "mag") { # use bulgeFrac IF given, ELSE estimate bulgeFrac from nSer
                  modellist0$sersic[[param]] = divide_magnitude(inheritModellist$sersic[[param]], frac = if (is.null(bulgeFrac)) get_B2T(inheritModellist$sersic$nser) else bulgeFrac)
                } else if (param == "re") { # arbitrary 0.333:1 (Bulge:Disk) divisions of effective radii
                  modellist0$sersic[[param]] = c(1/3*inheritModellist$sersic[[param]],inheritModellist$sersic[[param]])
                } else if (param == "nser") { # Assume nSer is describing towards bulge nSer; assume exponential disk.
                  modellist0$sersic[[param]] = c(inheritModellist$sersic[[param]],1)
                } else if (param == "axrat") { # Assume axial ratio is describing that of the disk only; assume spherical bulge
                  modellist0$sersic[[param]] = c(1,inheritModellist$sersic[[param]])
                }
              } # END IF param is to be inheritted
            } # END loop over parameters
          } # END 1comp --> 2comp 
        } else { # RData file did not exist
          cat(paste("\nWARNING: ",galName," does not has existing .RData file for: run='",inheritFrom$form,"'; band='",inheritFrom$band,"'; nComps='",inheritFrom$nComps,"'.\n",sep=""))
        }
        
        rm(envirFile,tempEnvir) # remove the temporary environment from memory
      }
      
      
      ####################################
      #####   DEFINE PROFIT INPUTS   #####
      ####################################
      if(verb){cat("INFO: Defining ProFit inputs.\n")}
      if (nComps == 1){
        ## SINGLE SERSIC PROFILE ##
        # modellist = list(
        #   sersic=list(
        #     xcen= ycenInits,
        #     ycen= xcenInits,
        #     mag=  magInits,
        #     re=   reInits,
        #     nser= nSerInits,
        #     ang=  angInits,
        #     axrat=axratInits, # Bulge is initially at axrat=1
        #     box=  boxInits # no boxiness
        #   )
        # )
        
        modellist = modellist0
        
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
            box=  fitBoxiness
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
        
        stdevs = list(
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
        
        priors = profitMakePriors(modellist=modellist, sigmas=stdevs, tolog=tolog, tofit=tofit, allowflat = TRUE) # allowflat allows for flat priors (e.g. boxiness) where the log-likelihood will be computed as 0, rather than -Inf.
        
        # The hard intervals should also be specified in linear space.
        intervals=list(
          sersic=list(
            xcen=list(lim=c(inits$xcen[mainIndex]-10,inits$xcen[mainIndex]+10)),
            ycen=list(lim=c(inits$ycen[mainIndex]-10,inits$ycen[mainIndex]+10)),
            mag=list(lim=c(7,zeroPoint)),
            re=list(lim=c(0.25,100)),
            nser=list(lim=c(0.25,20)),
            ang=list(lim=c(-180,360)),
            axrat=list(lim=c(0.05,1)),
            box=list(lim=c(-0.5,0.5))
          )
        )
        
        # END Single Sersic 
      } else if (nComps == 2){
        ## DOUBLE SERSIC PROFILE ##
        # modellist = list(
        #   sersic=list(
        #     xcen= c(inits$xcen[mainIndex],inits$xcen[mainIndex]),
        #     ycen= c(inits$ycen[mainIndex],inits$ycen[mainIndex]),
        #     mag = magInits,
        #     re=   reInits,
        #     nser= nSerInits,
        #     ang=  c(inits$ang[mainIndex],inits$ang[mainIndex]),
        #     axrat=c(1,inits$axrat[mainIndex]), # Bulge is initially at axrat=1
        #     box=c(0,0) # no boxiness
        #   )
        # )
        
        modellist = modellist0
        
        # Parameters to fit
        tofit=list(
          sersic=list(
            xcen= if (fixedCentres) c(TRUE,NA) else c(TRUE,TRUE),
            ycen= if (fixedCentres) c(TRUE,NA) else c(TRUE,TRUE),
            mag=  rep(TRUE,2), 
            re=   rep(TRUE,2),
            nser= c(TRUE,freeDisk),
            ang=  c(freeBulge,TRUE),
            axrat=c(freeBulge,TRUE),
            box = rep(fitBoxiness,2)
          )
        )
        
        # Define the parameters which should be fitted in log space
        tolog=list(
          sersic=list(
            xcen= rep(FALSE,2),
            ycen= rep(FALSE,2),
            mag=  rep(FALSE,2),
            re=   rep(TRUE,2), # re is best fit in log space
            nser= rep(TRUE,2), # nser is best fit in log space
            ang=  rep(FALSE,2), 
            axrat=rep(TRUE,2), # axrat is best fit in log space
            box=  rep(FALSE,2)
          )
        )
        
        
        # Define the sigmas object
        sigmaArr = c(2,2,5,1,1,30,0.3,Inf)
        
        stdevs = list(
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
        
        priors = profitMakePriors(modellist=modellist, sigmas=stdevs, tolog=tolog, tofit=tofit, allowflat = TRUE) # allowflat allows for flat priors (e.g. boxiness) where the log-likelihood will be computed as 0, rather than -Inf.
        
        # The hard intervals should also be specified in linear space.
        intervals=list(
          sersic=list(
            xcen=list(lim=c(inits$xcen[mainIndex]-10,inits$xcen[mainIndex]+10),lim=c(inits$xcen[mainIndex]-10,inits$xcen[mainIndex]+10)),
            ycen=list(lim=c(inits$ycen[mainIndex]-10,inits$ycen[mainIndex]+10),lim=c(inits$ycen[mainIndex]-10,inits$ycen[mainIndex]+10)),
            mag=list(lim=c(10,25),lim=c(10,25)),
            re=list(lim=c(0.25,100),lim=c(0.25,100)),
            nser=list(lim=c(1.25,20.0),lim=c(0.5,1.5)),
            ang=list(lim=c(-180,360),lim=c(-180,360)),
            axrat=list(lim=c(0.1,1),lim=c(0.1,1)),
            box=list(lim=c(-0.5,0.5),lim=c(-0.5,0.5)) 
          )
        )
      
      } # END Double Sersic
    
      if (subSky == FALSE){ # If sky was not subtracted, add a sky level to the modellist
        if(verb){cat("INFO: Adding sky component to model.\n")}
        modellist$sky = list(bg=skyEst$sky)
        tofit$sky = list(bg=TRUE)
        tolog$sky = list(bg=FALSE)
      }
      
      
      ### Setup Data ###
      if(verb){cat("INFO: Setting up Data object.\n")}
      Data = profitSetupData(image=image, psf=psf, segim=segMap, sigma=sigma,
                             modellist=modellist, tofit=tofit, tolog=tolog, intervals=intervals, priors = if (usePriors) priors else NULL,
                             magzero=zeroPoint, algo.func=fitMode, like.func=likeFunction, verbose=FALSE)
      
      
      ### Plot Input Model Likelihood ###
      if (output  && outputInitial){
        initLikelihoodFilename = paste(baseFilename,"_LikelihoodInitial.png",sep='')
        png(paste(outputDir,initLikelihoodFilename,sep='/'),width=1600,height=1000,pointsize=28)
        profitLikeModel(parm=Data$init,Data=Data,makeplots=TRUE,plotchisq=TRUE)
        dev.off()
      
        ### Plot Input Model Ellipse ###
        initEllipseFilename = paste(baseFilename,"_EllipseInitial.png",sep='')
        png(paste(outputDir,initEllipseFilename,sep='/'),width=1000,height=750,pointsize=20)
        if (nComps == 1){
          try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modellist),pixscale=pixScale,SBlim=25))
        } else if (nComps == 2){
          try(profitEllipsePlot(Data=Data,modellist=modellist,pixscale=pixScale,SBlim=25))
        }
        dev.off()
      }
      
      
      ##################################################
      #####   Improve Initial Guesses (Isophotal)  #####
      ##################################################
      
      ## Attempt to improve initial guesses via Isophote fitting:
      # Only run if running 2 components and the image does not contain NaN padding (i.e. galaxies on frame edges).
      if (improveInits == TRUE && nComps == 2 && is.null(inheritFrom)){
        if(verb){cat("INFO: Attempting to improve initial guess.\n")}
        if (output && outputIsophotes){
          isophotesFilename = paste(galName,"_",band,"_Isophotes.png",sep='')
          png(paste(outputDir,isophotesFilename,sep='/'),width=1050,height=500,pointsize=16)
          par(mfrow=c(1,2), mar=c(3.5,3.5,1,2))
        }
        
        # Get the ellipse isophotes
        ellipses = profoundGetEllipses(image,segim=segmentation$segim,segID=mainID,levels=20,pixscale=pixScale,magzero=zeroPoint,dobox=FALSE,plot=output)
        rMin = 1; rMax = 30; rDiff = 0.1
        rLocs=seq(rMin,rMax,by=rDiff)
        rCut = (7.5-rMin)/rDiff+1
        
        bulgeInit = profitRadialSersic(rLocs, mag=magInits[1], re=reInits[1], nser=nSerInits[1])
        diskInit = profitRadialSersic(rLocs, mag=magInits[2], re=reInits[2], nser=nSerInits[2])
        
        # Define minimisation function
        sumsq1D=function(par=c(magInits[1], log10(reInits[1]), log10(nSerInits[1]), magInits[2], log10(reInits[2]), log10(nSerInits[2])),rad, SB, pixscale=pixScale)
        {
          bulge=profitRadialSersic(rad, mag=par[1], re=10^par[2], nser=10^par[3])
          disk=profitRadialSersic(rad, mag=par[4], re=10^par[5], nser=10^par[6])
          total=profitFlux2SB(bulge+disk, pixscale=pixscale)
          return=sum((total-SB)^2)
        }
        
        ## Interval bounds:
        # Mag (B/D): 10 - 25
        # Re (B/D): 0.25 - 100
        # nSer (B): 1.25 - 20; nSer (D): 0.5 - 1.5
        lower = c(10,log10(0.25),log10(1.25),10,log10(0.25),log10(0.5))
        upper = c(25,log10(100), log10(20.0),25,log10(100), log10(1.5))
        
        # Run L-BFGS-B optimisation:
        if(verb){cat("INFO: Running BFGS 1D isophotal optimisation.\n")}
        isoFit=optim(sumsq1D, par=c(magInits[1], log10(reInits[1]), log10(nSerInits[1]), magInits[2], log10(reInits[2]), log10(nSerInits[2])),
                    rad=ellipses$ellipses$radhi[4:19], SB=ellipses$ellipses$SB[4:19], pixscale=pixScale,
                    method='L-BFGS-B', lower=lower, upper=upper)$par
        
        bulgeIsofit = profitRadialSersic(rLocs, mag=isoFit[1], re=10^isoFit[2], nser=10^isoFit[3])
        diskIsofit = profitRadialSersic(rLocs, mag=isoFit[4], re=10^isoFit[5], nser=10^isoFit[6])
        
        if (output && outputIsophotes) {
          magplot(ellipses$ellipses$radhi[4:19], ellipses$ellipses$SB[4:19],ylim=c(25,17), grid=TRUE, type='l',lwd=2.5,
                xlab='Radius (pixels)', ylab='Surface Brightness (mag/arcsec^2)',main="Bulge+Disk Surface Brightness Profile")
          points(ellipses$ellipses$radhi[4:19],ellipses$ellipses$SB[4:19],pch=16)
          lines(rLocs, profitFlux2SB(bulgeIsofit, pixscale=pixScale), col='red',lty='dashed',lwd=2)
          lines(rLocs, profitFlux2SB(diskIsofit, pixscale=pixScale), col='blue',lty='dashed',lwd=2)
          lines(rLocs, profitFlux2SB(bulgeIsofit+diskIsofit, pixscale=pixScale), col='green',lwd=2)
          dev.off()
        }
        
        # Check for divergence
        isoValid = TRUE # IF isoValid=TRUE THEN use these as initial guesses. ELSE run acesApproximation()
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
        
        # If Re or Magnitude parameters are divergent then the isophotal solution is not valid
        isoValid = if ((convParams['mag1']==T && convParams['mag2']==T) && (convParams['re1']==T && convParams['re2']==T)) TRUE else FALSE
        
        ### Check whether bulge dominates flux in outer regions.
        bulgeSB = profitFlux2SB(bulgeIsofit, pixscale=pixScale)
        diskSB = profitFlux2SB(diskIsofit, pixscale=pixScale)
        for (ii in seq(rCut,(rMax-rMin)/rDiff+1,1)){
          if (bulgeSB[ii] < diskIsofit[ii]){ # is bulge < disk? As lower values of surface magnitude means brigher
             isoValid = FALSE
             if(verb){cat(paste("WARNING: Bulge dominates beyond r_cut = ",toString(rCut*pixScale),"\" (@ r = ",toString((ii*rDiff+rCut)*pixScale),"\").\n",sep=''))}
             break
          }
        }
        
        # If solution has converged, set the isophotal fit solution as the initial guesses
        if (isoValid==TRUE){
          if(verb){cat("INFO: Replacing initial model with isophotal 1D fit solution.\n")}
          Data$init['sersic.mag1'] = isoFit[1]
          Data$init['sersic.re1'] = 10^isoFit[2]
          Data$init['sersic.nser1'] = 10^isoFit[3]
          Data$init['sersic.mag2'] = isoFit[4]
          Data$init['sersic.re2'] = 10^isoFit[5]
          Data$init['sersic.nser2'] = 10^isoFit[6]
        } else {
          if(verb){cat("WARNING: isoFit solution did not converge.\n")}
        }
      } else { # IF nComps = 1 THEN move onto LAFit
        isoValid = FALSE
      } # END iosphotal 1D fitting optimisation
      
      
      ##############################################################
      #####   Improve Initial Guesses (LaplaceApproximation)   #####
      ##############################################################
      
      # IF the result from 1D isophotal fitting did not converge THEN attempt a LaplaceApproximation() fit:
      if (improveInits==TRUE && isoValid==FALSE && is.null(inheritFrom)) { # LaplaceApproximation LM fit
        if(verb){cat("INFO: Attempting to improve inital guess with LaplaceApproximation()\n")}
        Data$algo.func = "LA" # Change optimising algorithm
        LAFit = LaplaceApproximation(profitLikeModel,parm = Data$init, Data = Data, Iterations=1e3, Method = 'LM', CovEst='Identity', sir = FALSE)
        
        # Remake model
        optimModel=profitRemakeModellist(LAFit$Summary1[,1],Data$modellist,Data$tofit,Data$tolog)$modellist
        
        ### Check for divergence
        if(verb){cat("INFO: Testing LAFit for divergence.\n")}
        LAValid = TRUE # IF converged=TRUE THEN use these as initial guess. ELSE run LaplacesApproximation()
        for (n in seq(nComps)){
          if (is.element(optimModel$sersic$mag[n],Data$intervals$sersic$mag[[n]])){LAValid=FALSE}
          if (is.element(optimModel$sersic$re[n],Data$intervals$sersic$re[[n]])){LAValid=FALSE}
        }
        
        if (nComps == 2){ # Create 1D profiles for bulge/disk
          rMin = 1; rMax = 30; rDiff = 0.1
          rLocs=seq(rMin,rMax,by=rDiff)
          rCut = (7.5-rMin)/rDiff+1
          
          bulgeLAfit = profitRadialSersic(rLocs, mag=optimModel$sersic$mag[1], re=optimModel$sersic$re[1], optimModel$sersic$nser[1])
          diskLAfit = profitRadialSersic(rLocs, mag=optimModel$sersic$mag[2], re=optimModel$sersic$re[2], optimModel$sersic$nser[2])
        
          ### Check whether bulge dominates flux in outer regions.
          bulgeSB = profitFlux2SB(bulgeLAfit, pixscale=pixScale)
          diskSB = profitFlux2SB(diskLAfit, pixscale=pixScale)
          for (ii in seq(rCut,(rMax-rMin)/rDiff+1,1)){
            if (bulgeSB[ii] < diskSB[ii]){
              LAValid = FALSE
              if(verb){cat(paste("WARNING: Bulge dominates beyond r > ",toString(rCut*pixScale),"\" (r = ",toString((rCut+ii*rDiff)*pixScale),"\").\n",sep=''))}
              break
            }
          }
        }
        
        # If solution has converged, set the LAFit solution as the initial guess
        if (LAValid==TRUE){
          if(verb){cat("INFO: Replacing initial model with LAFit solution.\n")}
          Data$init = LAFit$Summary1[,1]
        } else {
          if(verb){cat("INFO: LAFit solution did not converge.\n")}
        }
        
      } else {
        LAValid = FALSE
      } # END LAFit optimisation
      
      #####################################
      #####   Optimise Model (MCMC)   #####
      #####################################
      
      ### Laplaces Demon (Full MCMC)
      if(verb){cat(paste("\n\n* ",galName," * [band = ",band,"; comps = ",nComps,"]"," (",count,"/",length(galList),")\n\n",sep=""))}
      if (fitMode == "LD"){
        if(verb){cat("INFO: Starting LaplacesDemon() ~ Full MCMC optimisation.\n")}
        startTime = Sys.time()
        
        Data$algo.func = fitMode
        LDFit = LaplacesDemon(profitLikeModel, Initial.Values=Data$init, Data=Data, Iterations=MCMCIters, Algorithm=MCMCAlgo,Thinning=1,Specs=list(alpha.star=0.44), Status=MCMCStatus)
        
        ### Plot posterior distributions
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
        optimFit = bestLD[,1] # the optimised fitting parameters
        
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
      
      ##########################################
      ###   Remake intial/optimised models   ###
      ##########################################
      
      Data$usecalcregion = FALSE # Turn usecalcregion off to produce model over entire image.
      
      initImage = profitMakeModel(modellist,dim=dim(image),magzero=zeroPoint)$z
      
      remakeModel = profitRemakeModellist(parm = optimFit, Data = Data)
      optimModellist = remakeModel$modellist
      optimParams = remakeModel$parm
      
      optimModel = profitMakeModel(modellist = optimModellist,
                                   magzero = zeroPoint, psf = psf, dim = dim(image), 
                                   psfdim = dim(Data$psf), whichcomponents = list(sersic="all"), 
                                   rough = FALSE, magmu = Data$magmu)
      
      optimImage = optimModel$z # image of the optimised model
      
      ### Plot Model Images ###
      if(output && outputModel){
        noise = matrix( rnorm(dims[1]*dims[2],mean=0.0,sd=skyEst$skyRMS), dims[1], dims[2])
        
        initModelFilename = paste(baseFilename,"_ModelImage.png",sep='')
        png(paste(outputDir,initModelFilename,sep='/'),width=1200,height=800,pointsize = 20)
        par(mfrow=c(1,2), mar=c(0.4,0.4,1,1))
        
        magimage(image,axes=F,bad=0)
        text(0.025*dim(image)[1],0.95*dim(image)[2],"Image",col='white',pos=4,cex=1.5)
        
        magimage(optimImage+noise,axes=F,bad=0)
        text(0.025*dim(image)[1],0.95*dim(image)[2],"Optimised Model",col='white',pos=4,cex=1.5)
        text(0.875*dim(image)[1],0.05*dim(image)[2],expression(paste("N: ",mu, "=0.0, ", sigma,"=")),col='white',pos=2,cex=1.0)
        text(0.825*dim(image)[1],0.05*dim(image)[2],format(skyRMS,digits=3),col='white',pos=4,cex=1.0)
        
        dev.off()
      }
      
      ### Plot Optimisation results ###
      if(output  && outputOptimised){
        
        ### Plot Optimised Model Likelihood:
        optimLikelihoodFilename = paste(baseFilename,"_LikelihoodOptimised.png",sep='')
        png(paste(outputDir,optimLikelihoodFilename,sep='/'),width=1600,height=1000,pointsize=28)
        profitLikeModel(optimFit,Data,makeplots=TRUE,plotchisq=TRUE)
        dev.off()
      
      
        ### Plot Optimised Ellipse Plot:
        optimEllipseFilename = paste(baseFilename,"_EllipseOptimised.png",sep='')
        png(paste(outputDir,optimEllipseFilename,sep='/'),width=1000,height=750,pointsize=20)
        if (nComps == 1){
          try(profitEllipsePlot(Data=Data,add_pseudo_bulge(optimModellist),pixscale=pixScale,SBlim=25,FWHM=1.4))
        }else if (nComps == 2){
          try(profitEllipsePlot(Data=Data,optimModellist,pixscale=pixScale,SBlim=25,FWHM=1.4,raw=FALSE))
        }
        dev.off()
      
      }
      
      rTrun = numTrun*optimModellist$sersic$re
      
      if (nComps == 2){ # Calculate the Bulge-to-total ratio
        B2T = calc_B2T(optimModellist$sersic$mag[1],optimModellist$sersic$mag[2])
      }
      
      # Calculate various surface brightness measures (mag/arcsec^2)
      aveSBe = calc_ave_SB(mag = optimModellist$sersic$mag, re = optimModellist$sersic$re, pixScale = pixScale) # The average surf. brightness within the effective radius
      SBe = calc_SB(aveSB = aveSBe, n = optimModellist$sersic$nser) # The surf. brightness at the effective radius
      SB0 = calc_SB0(SB = SBe, n = optimModellist$sersic$nser) # The central surf. brightness
      
      # Calculate the magnitude of the model truncated at some radius.
      magTrun = truncate_mag(SB=SBe, re=optimModellist$sersic$re, n=optimModellist$sersic$nser, rTrun = rTrun, pixScale = pixScale)
      
      ### Calculate the chi^2 statistics
      chisq = calc_chisq(image, optimImage, sigma, segMap)
      
      ### The stationarity of the optimisation:
      # is stationary = 1; not stationary = 0
      stationarity = if (any(is.na(LDFit$Summary2))) 0 else 1
      
      
      ###########################
      #####  WRITE OUTPUTS  #####
      ###########################
      if(verb){cat("INFO: Writing outputs to file.\n")}
      if(output && outputResults){
        outputFilename = paste(baseFilename,"_Output.csv",sep='')
        write_output(file=paste(outputDir,outputFilename,sep='/'),name=galName,nComps=nComps,init=modellist$sersic,optim=optimModellist$sersic,chisq=chisq,time=elapsedTime,stat=stationarity)
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
        cat(paste("Run: ",run,"\n",sep=""))
        cat(paste("Flavour: ",flavour,"\n",sep=""))
        cat(paste("\nDescription:   ","\n",description,"\n",sep=""))
        cat(paste("\nDate: ",Sys.time(),"\n",sep=""))
        
        if (galFile != ""){
          cat(paste("\nGalaxy retrieved from file:  \n  ",galFile,"\n",sep=""))
          cat(paste("Line number = ",lineNum,"\n",sep=""))
        }
        
        cat(paste("\nElapsed time: ",elapsedTime,"\n",sep=""))
        
        cat("\n\n>> Inputs:\n")
        cat(paste("Base directory: ",outputDir,"\n",sep=""))
        cat(paste("Image: ",imgFile,"\n",sep=""))
        cat(paste(" Dimensions: ",dims[1]," x ",dims[2]," (",dims[1]*pixScale/60.0,"' x ",dims[1]*pixScale/60.0,"')","\n",sep=""))
        cat(paste(" Padding: ",padded,"\n",sep=""))
        if (loadPSF && file.exists(psfFile)){ cat(paste("PSF: ",psfFile,"\n",sep=""))}
        if (loadSegMap && file.exists(segMapFile)){ cat(paste("Segmentation map: ",segMapFile,"\n",sep=""))}
        if (loadSigma && file.exists(sigmaFile)){ cat(paste("Sigma map: ",sigmaFile,"\n",sep=""))}
        cat(paste("\nZero point: ",zeroPoint,"\n",sep=""))
        cat(paste("gain: ",gain,"\n",sep=""))
        if (dataSource == "SDSS") {cat(paste("Soft-bias: ",softBias,"\n",sep=""))}
        
        cat("\n\n>> Sky Statistics:\n")
        cat(gsub("\\[1\\]","\n",skyStats)) # Print the sky statistics: the output from maghist() of profoundMakeSkyGrid()
        
        cat("\n\n>> Segmentation:\n")
        cat(paste("Num. Objects: ",length(segmentation$segstats[[1]]),"\n",sep=""))
        cat(paste("Main Segment ID: ",mainID,"\n",sep=""))
        cat(paste("Main Segment index: ",mainIndex,"\n",sep=""))
        cat("\n Stats for target object:\n")
        print(segmentation$segstats[mainIndex,])
        
        if (!loadSegMap){ cat(paste("Num. expanded pixels in segment: ",numPixExpand,"\n",sep=""))}
        
        cat("\n\n>> Model controls:\n")
        cat("Fixed component centres: ",fixedCentres,"\n")
        cat("Fit for boxiness: ",fitBoxiness,"\n")
        cat("Free disk Sersic: ",freeDisk,"\n")
        cat("Free bulge shape: ",freeBulge,"\n")
        cat("Sersic index from Concentration: ",nFromCon,"\n")
        cat("Bulge/Total Ratio: ",bulgeFrac,"\n")
        cat("Used Priors: ",usePriors,"\n")
        
        cat("\n\n>> Initial Model:\n")
        print(modellist$sersic)
        
        cat("\n\n>> Fitting Parameters:\n")
        print(tofit$sersic)
        
        cat("\n\n>> Intervals:\n")
        print(intervals$sersic)
        
        if(improveInits==TRUE && is.null(inheritFrom)){
          cat("\n\n>> Attempted Improvements on Initial Guesses:\n")
          if(isoValid==TRUE){
            cat("Isophotal 1D fitting successful\n")
            cat("Fitting solutions:\n")
            for (ii in seq(length(convParams))){
              if (ii == 1 || ii == 4){cat(paste(" ",names(convParams)[ii],": ",toString(isoFit[ii]),'\n',sep=''))}
              else {cat(paste(" ",names(convParams)[ii],": ",toString(10^isoFit[ii]),'\n',sep=''))}
            }
          } else {
            cat("Isophotal 1D fitting not succesful\n")
          }
          
          if(isoValid==FALSE){
            if(LAValid==TRUE) {cat("\nLaplaceApproximation() 2D fitting successful\n")}
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
        print(optimModellist$sersic)
        if (subSky == FALSE) {print(optimModellist$sky)}
        
        cat(paste("\n chi^2 = ",sprintf("%.3f", chisq),"\n", sep=""))
        cat(paste("\n Stationarity = ",if (stationarity == 1) "TRUE" else "FALSE", sep=""))
        
        cat("\n\n>> Flags:\n")
        # Check if optimised parameters are stuck to interval bounds
        for (key in names(optimModellist$sersic)){
          for (n in seq(1,nComps)){
            if (optimModellist$sersic[[key]][n] %in% intervals$sersic[[key]][[n]] && tofit$sersic[[key]][n]){
              cat(paste(" - ",key,n,' = ',optimModellist$sersic[[key]][n],'\n',sep=''))
            }
          }
        }
        if(segmentation$segstats[mainIndex,]$edge_frac < 0.7){cat(paste('\n - Segmentation boundary = ',segmentation$segstats[mainIndex,]$edge_frac,'\n',sep=''))}
        if(segmentation$segstats[mainIndex,]$edge_excess > 1.0){cat(paste('\n - Segmentation edge excess = ',segmentation$segstats[mainIndex,]$edge_excess,'\n',sep=''))}
        if(segmentation$segstats[mainIndex,]$asymm > 0.2){cat(paste('\n - Segmentation asymmetry = ',segmentation$segstats[mainIndex,]$asymm,'\n',sep=''))}
        if(segmentation$segstats[mainIndex,]$flag_border != 0){cat(paste('\n - Segmentation borders with: ',segmentation$segstats[mainIndex,]$flag_border,'\n',sep=''))}
        
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

