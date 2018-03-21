### Script to loop through a list of galaxies, rerun the segmentation until and save when happy.


###################################################################
########################## LOAD PACKAGES ##########################
###################################################################

library(ProFit,warn.conflicts=FALSE) # Bayesian galaxy fitting tool
library(ProFound) # Source Extraction and Image Segmentation tool
library(EBImage,warn.conflicts=FALSE) # Image processing package
library(magicaxis) # "Magically Pretty Plots" ~ ASGR
library(FITSio) # .FITS file input/output


get_gal_list = function(galFile, lineNum=0) # Given the path to a file, extract the galaxy list at a particular line(s)
{
  # <param: galFile [string]> - The path to the file containing the galaxy list(s)
  # <param: lineNum [int]> - The line number at which to extract the list(s) (lineNum = 0 for all lines)
  
  # <return: galList [list]> - A list containing the galaxies to be optimised
  
  # Check whether galFile exists
  if (!file.exists(galFile)){
    cat(paste("ERROR: The input galFile: '",galFile,"' does not exist!\n -- ABORTING --\n",sep=""))
    quit(status=1)
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

find_main = function(seg) # Determine the main (central) source ID in the image
{
  # <param: seg [list]> - The segmentation object.
  
  # <return: mainID [int]> - The ID for the main (centremost) source.
  
  nSources = length(seg$segstats$segID) # get number of sources.
  
  dims = dim(seg$segim)
  x0 = dims[1]/2; y0 = dims[2]/2 # image centre position.
  sepArr = array(0,dim=nSources) # empty separation array
  
  # Calculate separation of source centres.
  for (ii in seq(1,nSources)) {sepArr[ii] = sqrt( (seg$segstats$xcen[ii] - x0)^2 + (seg$segstats$ycen[ii] - y0)^2 )}
  mainID = which.min(sepArr) # The main source is the one with the smallest separation from the centre
  return(mainID)
}



find_main_group = function(seg, mainID = NULL){ # Determine the main (central) group num.; i.e. contains the main.
  # <param: seg [list]> - The segmentation object.
  # <param: mainID [int]> - The ID of the target segment (If not given, will be determined as the centre-most segment)
  
  # <return: mainGroupID [int]> - The ID for the main (centre-most) group.
  
  if (is.null(mainID)){
    mainID = find_main(seg)
  }
  
  nGroups = length(seg$group$groupsegID[,1])
  
  for (id in seq(1,nGroups)){
    if (mainID %in% seg$group$groupsegID[id,]$segID[[1]]) {
      mainGroupID = id
      break
    }
  }
  
  return(mainGroupID)
}

remove_stars = function(seg, list){
  # <param: seg [list]> - The segmentation object.
  # <param: list [list (int)]> - The list of segment IDs to check for point-sources
  
  # <return: newList [list (int)]> - The original list with point-sources removed.
  
  newList = list
  for (id in list){
    if (!is.na(seg$segstats[id,]$con)){
      if (seg$segstats[id,]$con < 0.55 && seg$segstats[id,]$semimaj * 2.0 < 10.0){
        newList = setdiff(newList,id)
      }
    }
  }
  
  return(newList)
  
}

combine_segments = function(seg, list, mainID=NULL, plot=FALSE){ # combines a list of segments together
  # <param: seg [list]> - The segmentation object.
  # <param: list [list (int)]> - The list of segment IDs to be combined
  # <param: mainID [int]> - The ID of the target segment (If not given, will be determined as the centre-most segment)
  # <param: plot [boolean]> - Whether to plot the resulting segmentation map.
  
  # <return: mainID [int]> - The combined segmentation object.
  
  if (is.null(mainID)){
    mainID = find_main(seg)
  }
  
  for (id in list){
    seg$segim[seg$segim == id] = mainID
  }
  
  if(plot){profoundSegimPlot(image,seg$segim)}
  
  return(seg)
}



find_embedded = function(seg0, segDiff, mainID=NULL){ # Find embedded sources; i.e. sources which are within the original combined/expanded segment but not the central core.
  # <param: seg0 [list]> - The main segmentation object.
  # <param: segDiff [list]> - The segmentation object from the difference map.
  # <param: mainID [int]> - The ID of the target segment (If not given, will be determined as the centre-most segment)
  
  # <return: embeddedIDs [int]> - The list of embedded segment IDs.
  
  if (is.null(mainID)){
    mainID = find_main(seg0)
  }
  
  embeddedIDs = c()
  for (id in segDiff$segstats[,"segID"]){
    if (any( (seg0$segim==mainID) * (segDiff$segim==id) )){ # Check for any overlap of difference segment with main segment
      embeddedIDs = c(embeddedIDs,id)
    }
  }
  
  # Calculate separation of source centres.
  dims = dim(seg0$segim)
  x0 = dims[1]/2; y0 = dims[2]/2 # image centre position.
  
  sepArr = array(0.0,dim=length(embeddedIDs))
  for (ii in seq(1,length(embeddedIDs))) {sepArr[ii] = sqrt( (segDiff$segstats$xcen[embeddedIDs[ii]] - x0)^2 + (segDiff$segstats$ycen[embeddedIDs[ii]] - y0)^2 )}
  
  # Assume that the galaxy core should be ~ within 10 pixels of the centre, if so, remove this as an 'embedded segment'
  if (length(sepArr) != 0){
    centreID = which.min(sepArr) # Get the centre-most embedded segment (assumed to be the galaxy core)
    if (sepArr[centreID] <= 10.0){embeddedIDs = setdiff(embeddedIDs,embeddedIDs[centreID])} 
  }
  
  return(embeddedIDs)
}



mask_segments = function(seg0, segMask, list, plot=FALSE){
  # <param: seg0 [list]> - The main segmentation object.
  # <param: segMask [list]> - The masking segmentation object.
  # <param: list [list (int)]> - The list of segment IDs to mask the pixels
  # <param: plot [boolean]> - Whether to plot the resulting segmentation map.
  
  # <return: mainID [int]> - The combined segmentation object.
  
  maxID = max(seg0$segstats$segID)
  
  count = 0
  for (id in list){
    count = count + 1
    seg0$segim[segMask$segim == id] = maxID + count
  }
  
  if(plot){profoundSegimPlot(image, segim=seg0$segim)}
  
  return(seg0)
}



label_segments = function(seg){ # Labels the segments on a segmentation contour plot with their respective segment IDs
  # <param: seg [list]> - The segmentation object.
  
  # <return: NULL>
  
  col = rainbow(max(seg$segim), end = 2/3)
  for (ii in seg$segstats[,"segID"]){
    text(x=seg$segstats$xcen[ii], y=seg$segstats$ycen[ii], labels=seg$segstats$segID[ii], col = col[ii], pos = 2)
  }
  
}


##########################################################################################################

configFile = "/home/robincook/Google Drive/PhD/Fitting/Runs/Beta/Beta.conf"
source(configFile)

band = "r"
nComps = 2
run = "Beta"

galsDir = "/home/robincook/Documents/PhD/GASS/Galaxies"

#galFile = "/home/robincook/Google Drive/PhD/Fitting/Runs/Beta/Beta_redo_IDs.txt"
#galList = get_gal_list(galFile)

galList = c("GASS108026")

if (!exists("current")){current = 1}

toPlot = TRUE

###########################################################################################################

#for (gal in galList){
gal = galList[current]
print(gal)
  
envirFile = paste(gal,"-",run,"_",band,"_",nComps,"comp_WorkSpace.RData",sep="")
envirPath = paste(galsDir,gal,"Fitting",run,envirFile,sep="/")

if (file.exists( envirPath )){
  tempEnvir = new.env() # Create temporary envronment to place workspace of inheriting optimisation run.
  load(envirPath, envir=tempEnvir) # load inheriting RData into temporary environment
}
galName = gal

outputDir = paste(galsDir,gal,"Fitting",run,sep="/")

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
skyMap0 = profoundProFound(image0, skycut=1.0, tolerance=5, redosky=TRUE, redoskysize=21, pixcut=9,ext = 5,
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

#magimage(image)


segOutFilename = paste(gal,"-",run,"_",band,"_",nComps,"comp_SegMap.png",sep="")
segMapFilename = paste(gal,"_",band,"_SegMap.fits",sep="")


### Initial Segmentation ###
segmentation0 = profoundProFound(image, sigma=1.25, skycut=segSkyCut, tolerance=segTol, ext=5.0,
                                 magzero=zeroPoint, gain=gain, #header=header,
                                 stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)

# Find the main (central) source
mainID = find_main(segmentation0) # The main source is the one with the smallest separation from the centre

### Expansion of all Segments ###
segmentationExp = segmentation0

segmentationExp = profoundMakeSegimExpand(image=image, segim=segmentation0$segim, expand=mainID, skycut=-1.0, sigma=2.5,
                                          sky=0.0,skyRMS=skyMap$skyRMS,
                                          magzero=zeroPoint, gain=gain, #header=header,
                                          stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
# segmentationExp = profoundMakeSegimDilate(image=image, segim=segmentationExp$segim, expand=mainID, size=15, sigma=2.5,
#                                           sky=0.0,skyRMS=skyMap$skyRMS,
#                                           magzero=zeroPoint, gain=gain, #header=header,
#                                           stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)

# Find the main (central) source
mainID = find_main(segmentationExp) # The main source is the one with the smallest separation from the centre

### Combine Grouped Segment ###
segmentationExp$group = profoundSegimGroup(segmentationExp$segim)

mainGroupID = find_main_group(segmentationExp,mainID)
mainGroupList = segmentationExp$group$groupsegID[mainGroupID,]$segID[[1]]
#writeFITSim(segmentationExp$segim, file = paste(galsDir,gal,band,"Test_Seg.fits",sep='/'))

label_segments(segmentationExp)
#combineList = setdiff(mainGroupList,c(10,30))

combineList = remove_stars(segmentationExp,mainGroupList)

segmentationCom = combine_segments(segmentationExp,combineList,mainID,plot=TRUE)

### Point-source Segmentation ###
tempImage = image
if (any(is.nan(image))){
   tempImage[is.nan(tempImage)] = 0
}

imDiff = profoundImDiff(tempImage,sigma=2)
imDiff[is.nan(image)] = NaN

segmentationDiff = profoundProFound(imDiff,sigma=segSigma, skycut=segSkyCut,tolerance=segTol,
                                    magzero=zeroPoint, gain=gain, #header=header,
                                    stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)


mainDiffID = find_main(segmentationDiff) # The main source is the one with the smallest separation from the centre

### Expand point-source Segmentation ###
segmentationDiffExp = profoundMakeSegimDilate(image=image, expand = segmentationDiff$segstats[,"segID"],
                                              segim=segmentationDiff$segim, size = 21,
                                              sky=0.0,skyRMS=skyMap$skyRMS,
                                              magzero=zeroPoint, gain=gain, #header=header,                                               
                                              stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)


# segmentationDiffExp = profoundMakeSegimExpand(image=image, expand = segmentationDiff$segstats[,"segID"],
#                                            segim=segmentationDiff$segim, skycut=0.4, sigma=expSigma,
#                                            sky=0.0,skyRMS=skyMap$skyRMS,
#                                           magzero=zeroPoint, gain=gain, #header=header,
#                                           stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)


### Merge/Mask Segmentations
### maskIDs = find_embedded(segmentationCom,segmentationDiffExp)

#label_segments(segmentationDiffExp)
#maskIDs = setdiff(maskIDs,c(13))

#segmentationMerge = list()
#segmentationMerge$segim = profoundSegimMerge(image,segmentationCom$segim,segmentationDiffExp$segim)

#writeFITSim(segmentationDiffExp$segim, file = paste(galsDir,gal,band,segMapFilename,sep='/'))
###segmentationMerge = mask_segments(segmentationCom,segmentationDiffExp,maskIDs,plot=toPlot)

# @Robin Cook: Expand Segmentation image
if (exists("segmentationMerge")){
  segmentationMergeExp = profoundMakeSegimExpand(image=image, segim=segmentationMerge$segim, expand=mainID, skycut=expSkyCut, sigma=expSigma,
                                            sky=-1.0,skyRMS=skyMap$skyRMS,
                                            magzero=zeroPoint, gain=gain, #header=header,
                                            stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=toPlot)
} else if (exists("segmentationCom")) {
  segmentationMergeExp = profoundMakeSegimExpand(image=image, segim=segmentationCom$segim, expand=mainID, skycut=expSkyCut, sigma=expSigma,
                                                 sky=-1.0,skyRMS=skyMap$skyRMS,
                                                 magzero=zeroPoint, gain=gain, #header=header,
                                                 stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
} else {
  segmentationMergeExp = profoundMakeSegimExpand(image=image, segim=segmentationExp$segim, expand=mainID, skycut=expSkyCut, sigma=expSigma,
                                                 sky=-1.0,skyRMS=skyMap$skyRMS,
                                                 magzero=zeroPoint, gain=gain, #header=header,
                                                 stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
}


# @Robin Cook: Dilate Segmentation image
if (dilateSize != 0){ # if dilation > 0, perform dilation.
  segmentation = profoundMakeSegimDilate(image=image, segim=segmentationMergeExp$segim, expand=mainID, size=35,
                                         magzero=zeroPoint, gain=gain, #header=header,
                                         stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
} else {
  segmentation = segmentationExp # set final segmentation undilated
}


if (exists("segmentationMerge")){
  png(paste(outputDir,segOutFilename,sep='/'),width=1200,height=800,pointsize=16)
  par(mfrow=c(2,3), mar=c(0.25,0.25,0.25,0.25))
  profoundSegimPlot(image, segim=segmentation0$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Initial Segmentation",pos=4,col='white',cex= 1.75)
  profoundSegimPlot(image, segim=segmentationExp$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Segmentation Expansion",pos=4,col='white',cex= 1.75)
  profoundSegimPlot(image, segim=segmentationCom$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Central Group Combining",pos=4,col='white',cex= 1.75)
  profoundSegimPlot(imDiff, segim=segmentationDiffExp$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Embedded Source Detection",pos=4,col='white',cex= 1.75)
  profoundSegimPlot(image, segim=segmentationMergeExp$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Embedded Source Masking",pos=4,col='white',cex= 1.75)
  profoundSegimPlot(image, segim=segmentation$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Final Segmentation",pos=4,col='white',cex= 1.75)
  dev.off()
} else {
  png(paste(outputDir,segOutFilename,sep='/'),width=1200,height=800,pointsize=16)
  profoundSegimPlot(image, segim=segmentation$segim,axes=FALSE,bad=0); text(0.025*dim(image)[1],0.925*dim(image)[2],"Final Segmentation",pos=4,col='white',cex= 1.75)
  dev.off()
}

#writeFITSim(oldSegMap$segim, file = paste(galsDir,gal,band,segMapFilename,sep='/'))
rm(envirFile,envirPath,zeroPoint,gain,image,skyMap,dims,outputDir,segOutFilename,segMapFilename,
   segmentation0,segmentationExp,segmentationCom,segmentationDiff,segmentationDiffExp,segmentationMerge,segmentationMergeExp,segmentation,
   mainDiffID,mainID,mainGroupID,maskIDs,mainGroupList,tempImage,imDiff)


#}
###########################################################################################################
current = current + 1

