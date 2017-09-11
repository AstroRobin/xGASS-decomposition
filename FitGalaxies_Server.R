###
# Script to cycle through galaxy fitting for a list of specified galaxies in ugriz frequency bands for 1 and 2 components.
#
# Author: Robin Cook
# Date: 05/05/17
###

.libPaths(c("/home/arobotham/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))

library(ProFit)
library(EBImage)
library(magicaxis)
library(FITSio)
library(LaplacesDemon)

# Set up parallelism
library(doParallel)
library(foreach)
registerDoParallel(cores=8)

library(rgl)
library(rpanel)

###################################################################
######################### DEFINE CONSTANTS ########################
###################################################################

# Get home directory
#HOME = "/home" OR "/Users"
HOME = paste("/",unlist(strsplit(getwd(), '/'))[2],sep="")

#GALS_DIR = paste(HOME,"/robincook/Google Drive/PhD/GASS/Galaxies",sep = "")
GALS_DIR = paste("home/rcook/Documents/PhD/GASS/Galaxies",sep = "")
PIXSCALE = 0.396 ### Specifically for SDSS

###################################################################

###################################################################
######################### DEFINE FUNCTIONS ########################
###################################################################

plotSegIm = function(image, segim, mask, sky = 0, ...) 
{
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
  #print(sepArr)
  mainID = which.min(sepArr) # The main source is the one with the smallest separation from the centre
  return(mainID)
}

write_component = function(model,comp,state) # Write a model sersic component to file
{
  # <param: model [list]> - The modellist$sersic modellist
  # <param: comp [int]> - The component to be written
  # <param: state [char] - The status of the model [initial | optimised] >
  # <return: None>
  
  cat(paste(comp,
            state,
            model$xcen[comp],
            model$ycen[comp],
            model$mag[comp],
            model$re[comp],
            model$nser[comp],
            model$ang[comp],
            model$axrat[comp],
            model$box[comp],
            sep=','),append=T)
  cat('\n',append=T)
}

write_output = function(file,name,nComps,init,optim){ # Write optimisation result to file
  # <param: file [str]> - The file to append the results to.
  # <param: name [str]> - The name of the galaxy.
  # <param: nComps [int (1|2)]> - The number of components for the model (1 or 2).
  # <param: init [list]> - A list of initial values for the fit.
  # <param: optim [list]> - A list of optimised values.
  # <return: None>
  
  # Determine the number of columns in the table.
  nCols = 2*nComps*8 + 1 # 2 (inital/optimised) * nComps (1|2) * 8 (num. parameters in Sersic Profile) + 1 (Galaxy ID)
  
  if (nComps == 1){
    write(c("GASSID","x_in","y_in","mag_in","re_in","nser_in","ang_in","axrat_in","box_in","x_out","y_out","mag_out","re_out","nser_out","ang_out","axrat_out","box_out"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name,
            init$xcen[1],
            init$ycen[1],
            init$mag[1],
            init$re[1],
            init$nser[1],
            init$ang[1],
            init$axrat[1],
            init$box[1],
            optim$xcen[1],
            optim$ycen[1],
            optim$mag[1],
            optim$re[1],
            optim$nser[1],
            optim$ang[1],
            optim$axrat[1],
            optim$box[1]),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  } else if (nComps == 2){
    write(c("GASSID","x1_in","x2_in","y1_in","y2_in","mag1_in","mag2_in","re1_in","re2_in","nser1_in","nser2_in","ang1_in","ang2_in","axrat1_in","axrat2_in","box1_in","box2_in","x1_out","x2_out","y1_out","y2_out","mag1_out","mag2_out","re1_out","re2_out","nser1_out","nser2_out","ang1_out","ang2_out,axrat1_out,axrat2_out,box1_out,box2_out"),file=file,ncolumns=nCols,sep=',',append=FALSE)
    write(c(name,
            init$xcen[1],init$xcen[2],
            init$ycen[1],init$ycen[2],
            init$mag[1],init$mag[2],
            init$re[1],init$re[2],
            init$nser[1],init$nser[2],
            init$ang[1],init$ang[2],
            init$axrat[1],init$axrat[2],
            init$box[1],init$box[2],
            optim$xcen[1],optim$xcen[2],
            optim$ycen[1],optim$ycen[2],
            optim$mag[1],optim$mag[2],
            optim$re[1],optim$re[2],
            optim$nser[1],optim$nser[2],
            optim$ang[1],optim$ang[2],
            optim$axrat[1],optim$axrat[2]
            ,optim$box[1],optim$box[2]),
          file=file,ncolumns=nCols,sep=',',append=TRUE)
  }
  
}

append_output = function(file,name,nComps,init,optim) # Append optimisation results to master file
{
  # <param: file [str]> - The file to append the results to.
  # <param: name [str]> - The name of the galaxy.
  # <param: nComps [int (1|2)]> - The number of components for the model (1 or 2).
  # <param: init [list]> - A list of initial values for the fit.
  # <param: optim [list]> - A list of optimised values.
  # <return: None>
  
  # Check if file already exists; if not, create it!
  if (!file.exists(file)){
    sink(resultFilename,type="output")
    if (nComps == 1){
      writeLines("GASSID,x_in,y_in,mag_in,re_in,nser_in,ang_in,axrat_in,box_in,x_out,y_out,mag_out,re_out,nser_out,ang_out,axrat_out,box_out")
    }else if (nComps == 2){
      writeLines("GASSID,x1_in,x2_in,y1_in,y2_in,mag1_in,mag2_in,re1_in,re2_in,nser1_in,nser2_in,ang1_in,ang2_in,axrat1_in,axrat2_in,box1_in,box2_in,x1_out,x2_out,y1_out,y2_out,mag1_out,mag2_out,re1_out,re2_out,nser1_out,nser2_out,ang1_out,ang2_out,axrat1_out,axrat2_out,box1_out,box2_out")
    }
    sink()
  }
  
  lines = read.csv(file,stringsAsFactors=FALSE)
  index = which(apply(lines, 1, function(x) any(grepl(name, x))))
  if (length(index)[1]!=0){
    if (nComps==2){
      lines[index,] = c(name,init$xcen[1],init$xcen[2],init$ycen[1],init$ycen[2],init$mag[1],init$mag[2],init$re[1],init$re[2],init$nser[1],init$nser[2],init$ang[1],init$ang[2],init$axrat[1],init$axrat[2],init$box[1],init$box[2],optim$xcen[1],optim$xcen[2],optim$ycen[1],optim$ycen[2],optim$mag[1],optim$mag[2],optim$re[1],optim$re[2],optim$nser[1],optim$nser[2],optim$ang[1],optim$ang[2],optim$axrat[1],optim$axrat[2],optim$box[1],optim$box[2])
    }else{
      lines[index,] = c(name,init$xcen[1],init$ycen[1],init$mag[1],init$re[1],init$nser[1],init$ang[1],init$axrat[1],init$box[1],optim$xcen[1],optim$ycen[1],optim$mag[1],optim$re[1],optim$nser[1],optim$ang[1],optim$axrat[1],optim$box[1])
    }
  }else{
    if (nComps==2){
      lines[nrow(lines) + 1,] = c(name,init$xcen[1],init$xcen[2],init$ycen[1],init$ycen[2],init$mag[1],init$mag[2],init$re[1],init$re[2],init$nser[1],init$nser[2],init$ang[1],init$ang[2],init$axrat[1],init$axrat[2],init$box[1],init$box[2],optim$xcen[1],optim$xcen[2],optim$ycen[1],optim$ycen[2],optim$mag[1],optim$mag[2],optim$re[1],optim$re[2],optim$nser[1],optim$nser[2],optim$ang[1],optim$ang[2],optim$axrat[1],optim$axrat[2],optim$box[1],optim$box[2])
    }else{
      lines[nrow(lines) + 1,] = c(name,init$xcen[1],init$ycen[1],init$mag[1],init$re[1],init$nser[1],init$ang[1],init$axrat[1],init$box[1],optim$xcen[1],optim$ycen[1],optim$mag[1],optim$re[1],optim$nser[1],optim$ang[1],optim$axrat[1],optim$box[1])
    }
  }

  write.csv(lines,file,quote=FALSE,row.names=FALSE)
}


add_pseudo_bulge = function(model) # Add a zero-point magnitude bulge
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
compList = c(1)

### Specify frequency bands
#bandList = c('u','g','r','i','z')
bandList = c('r')

### Specify any Prefixes and descriptions to the output filename:
prefix = "SingleMCMC"
flavour = "fullsample"
description = "Running an MCMC fit on a sample of ~67 galaxies which were previously not able to be cropped."

### Specify which galaxies to fit. (Requires image and PSF files.)
### Galaxies in the Tail end of the sample with good Cropped images (n=728)
galList = c('GASS109079','GASS20448','GASS109035','GASS109034','GASS109135','GASS109050','GASS33737','GASS109088','GASS109132','GASS26151','GASS18755','GASS14863','GASS5035','GASS110038','GASS33777','GASS14831','GASS18830','GASS110036','GASS14943','GASS26221','GASS110058','GASS18900',
            'GASS110047','GASS18875','GASS18877','GASS18872','GASS110060','GASS18887','GASS26407','GASS26406','GASS110056','GASS54233','GASS110004','GASS18862','GASS54240','GASS26311','GASS26503','GASS26336','GASS26319','GASS26436','GASS110040','GASS23029','GASS23026','GASS26535','GASS5204',
            'GASS23039','GASS23070','GASS110019','GASS110057','GASS23102','GASS54577','GASS110042','GASS23114','GASS110003','GASS23120','GASS23194','GASS23195','GASS26602','GASS26598','GASS110017','GASS23203','GASS26586','GASS23213','GASS54763','GASS26569','GASS110059','GASS23227','GASS23228','GASS110041',
            'GASS15181','GASS26629','GASS15211','GASS26640','GASS26639','GASS15155','GASS15166','GASS23315','GASS15151','GASS110049','GASS8914','GASS110054','GASS38255','GASS29594','GASS110007','GASS15257','GASS110022','GASS29510','GASS8971','GASS29505','GASS15242','GASS23381','GASS23419','GASS34723','GASS23437',
            'GASS8953','GASS8945','GASS23408','GASS110063','GASS23445','GASS9085','GASS23450','GASS23496','GASS54986','GASS23453','GASS17659','GASS110081','GASS9114','GASS17635','GASS17673','GASS57415','GASS5442','GASS17684','GASS17622','GASS29624','GASS23518','GASS34989','GASS15570','GASS23563','GASS29699','GASS111042',
            'GASS111035','GASS47677','GASS48356','GASS47825','GASS48205','GASS48160','GASS111062','GASS17840','GASS12371','GASS15612','GASS111066','GASS17824','GASS48208','GASS12318','GASS17843','GASS23539','GASS5701','GASS111101','GASS111064','GASS48521','GASS48518','GASS24496','GASS12452','GASS48544','GASS12460',
            'GASS12458','GASS29842','GASS5848','GASS111083','GASS24513','GASS23685','GASS29883','GASS23715','GASS23703','GASS48604','GASS17865','GASS29892','GASS111030','GASS29871','GASS48974','GASS6015','GASS111079','GASS57697','GASS23761','GASS111100','GASS111086','GASS111071','GASS111080','GASS23789','GASS23781',
            'GASS111029','GASS48994','GASS111047','GASS12597','GASS111036','GASS111024','GASS111090','GASS111002','GASS111057','GASS111044','GASS18087','GASS111070','GASS23815','GASS18084','GASS18004','GASS111054','GASS111075','GASS111088','GASS111091','GASS111023','GASS111027','GASS111004','GASS111049','GASS111074',
            'GASS111063','GASS111007','GASS111056','GASS111028','GASS112065','GASS112090','GASS18138','GASS18185','GASS112072','GASS112106','GASS112014','GASS49727','GASS18131','GASS112027','GASS18225','GASS18220','GASS18202','GASS112088','GASS112019','GASS112117','GASS24149','GASS30192','GASS112058','GASS30175','GASS24168',
            'GASS24736','GASS24740','GASS24183','GASS112113','GASS112100','GASS112078','GASS18421','GASS18279','GASS24094','GASS24757','GASS18422','GASS30338','GASS6375','GASS28062','GASS18482','GASS112044','GASS50404','GASS30439','GASS18581','GASS12967','GASS12970','GASS24426','GASS18576','GASS12966','GASS50406','GASS28143',
            'GASS112038','GASS30471','GASS112051','GASS112109','GASS24366','GASS24436','GASS24437','GASS12983','GASS28168','GASS112053','GASS50550','GASS112080','GASS28482','GASS13037','GASS112081','GASS112054','GASS28462','GASS28461','GASS40393','GASS6565','GASS112089','GASS13074','GASS112016','GASS50856','GASS112116','GASS50866',
            'GASS40495','GASS35497','GASS28551','GASS40494','GASS112004','GASS112003','GASS112017','GASS112006','GASS28526','GASS112060','GASS112035','GASS112084','GASS30508','GASS35475','GASS30532','GASS112011','GASS13227','GASS113040','GASS113032','GASS113060','GASS113100','GASS113051','GASS113004','GASS35437','GASS6679','GASS113078',
            'GASS113049','GASS113052','GASS113001','GASS40439','GASS13156','GASS13158','GASS25154','GASS113137','GASS6749','GASS113126','GASS25215','GASS113022','GASS40570','GASS113141','GASS40566','GASS25213','GASS25214','GASS113098','GASS25202','GASS40686','GASS26911','GASS51025','GASS51094','GASS40781','GASS113000','GASS40790',
            'GASS113110','GASS113056','GASS44354','GASS113062','GASS113124','GASS51150','GASS113024','GASS113038','GASS51161','GASS34531','GASS25448','GASS40024','GASS40014','GASS13512','GASS113118','GASS43963','GASS35659','GASS44021','GASS43929','GASS51190','GASS13551','GASS13549','GASS113005','GASS7050','GASS7025','GASS13624',
            'GASS113010','GASS113134','GASS113063','GASS113003','GASS44856','GASS44846','GASS40317','GASS113057','GASS44892','GASS13618','GASS113008','GASS113105','GASS44889','GASS113109','GASS113067','GASS13674','GASS113123','GASS113150','GASS40257','GASS113122','GASS113011','GASS40247','GASS38409','GASS114141','GASS114137',
            'GASS114001','GASS114098','GASS114113','GASS114053','GASS9301','GASS114072','GASS114061','GASS114039','GASS114047','GASS114075','GASS114078','GASS38458','GASS25575','GASS7121','GASS25572','GASS41303','GASS30675','GASS114055','GASS114123','GASS114115','GASS114122','GASS114057','GASS30746','GASS38529','GASS114110',
            'GASS38538','GASS114033','GASS38546','GASS114031','GASS114082','GASS114027','GASS7286','GASS114046','GASS9384','GASS114058','GASS38462','GASS41444','GASS38472','GASS9343','GASS7310','GASS114054','GASS114127','GASS41323','GASS30854','GASS30847','GASS7405','GASS114091','GASS30811','GASS45613','GASS9507','GASS114121',
            'GASS9483','GASS114085','GASS114042','GASS114041','GASS114125','GASS9463','GASS114025','GASS114077','GASS114117','GASS114056','GASS114076','GASS7457','GASS114038','GASS7493','GASS114065','GASS114067','GASS114136','GASS9551','GASS114012','GASS45940','GASS114036','GASS28703','GASS114005','GASS114129','GASS114142',
            'GASS114008','GASS114074','GASS9619','GASS114037','GASS114144','GASS114010','GASS114090','GASS114000','GASS114048','GASS114024','GASS7509','GASS28875','GASS9607','GASS9604','GASS38198','GASS41482','GASS7520','GASS9814','GASS9601','GASS28810','GASS9585','GASS9776','GASS9572','GASS46068','GASS31095','GASS29090',
            'GASS7581','GASS9748','GASS41621','GASS9917','GASS9702','GASS41703','GASS9938','GASS41699','GASS9695','GASS31131','GASS9942','GASS41718','GASS31478','GASS9948','GASS29205','GASS41723','GASS29225','GASS9863','GASS38751','GASS38758','GASS38752','GASS29420','GASS38748','GASS29371','GASS29377','GASS38717',
            'GASS38716','GASS38718','GASS10032','GASS38729','GASS38703','GASS38732','GASS10031','GASS38721','GASS38706','GASS10019','GASS10010','GASS10040','GASS10012','GASS10005','GASS42191','GASS10058','GASS38935','GASS41783','GASS38923','GASS10132','GASS10145','GASS41793','GASS38909','GASS41743','GASS39014',
            'GASS39082','GASS41869','GASS39119','GASS39120','GASS41924','GASS41863','GASS10218','GASS42015','GASS10211','GASS39270','GASS7813','GASS10297','GASS42084','GASS10292','GASS39211','GASS42025','GASS42020','GASS42017','GASS41969','GASS42140','GASS10367','GASS42013','GASS42141','GASS41970','GASS10358',
            'GASS10404','GASS42156','GASS39346','GASS42174','GASS10447','GASS42167','GASS39469','GASS46564','GASS39311','GASS39465','GASS25057','GASS39407','GASS39532','GASS39605','GASS39607','GASS39595','GASS39606','GASS39567','GASS39600','GASS46678','GASS28348','GASS28365','GASS42287','GASS26958','GASS25682',
            'GASS47221','GASS42402','GASS25721','GASS21023','GASS31156','GASS47405','GASS10817','GASS10827','GASS10831','GASS10856','GASS10813','GASS10836','GASS10850','GASS10863','GASS10841','GASS10885','GASS10889','GASS10890','GASS10844','GASS10872','GASS10896','GASS10918','GASS10884','GASS10942','GASS10943',
            'GASS10944','GASS10950','GASS10948','GASS10953','GASS10952','GASS122005','GASS10990','GASS10985','GASS11016','GASS122001','GASS11013','GASS11015','GASS11003','GASS122000','GASS11041','GASS122004','GASS11086','GASS11087','GASS11080','GASS11071','GASS11126','GASS11231','GASS11120','GASS11223','GASS11250',
            'GASS11249','GASS11257','GASS11268','GASS11267','GASS11314','GASS11312','GASS11311','GASS11340','GASS11193','GASS11298','GASS11297','GASS11295','GASS11192','GASS11303','GASS123006','GASS11285','GASS11280','GASS11284','GASS11292','GASS11291','GASS11347','GASS123010','GASS11462','GASS11349','GASS123003',
            'GASS11437','GASS11410','GASS11435','GASS11514','GASS11434','GASS11395','GASS11524','GASS11397','GASS11585','GASS11544','GASS11383','GASS11366','GASS11357','GASS11669','GASS11582','GASS11685','GASS11687','GASS11571','GASS11573','GASS11568','GASS11567','GASS11697','GASS123005','GASS123007','GASS123008',
            'GASS11754','GASS11791','GASS11783','GASS11808','GASS11845','GASS11824','GASS3318','GASS4037','GASS4045','GASS4054','GASS31592','GASS36169','GASS31775','GASS108144','GASS17135','GASS108106','GASS32257','GASS32568','GASS32619','GASS25852','GASS22436','GASS109112','GASS33039','GASS8335','GASS4914',
            'GASS23088','GASS48369','GASS5872','GASS29898','GASS111031','GASS112079','GASS30332','GASS13005','GASS25327','GASS947','GASS35649','GASS44718','GASS44710','GASS13666','GASS44942','GASS1115','GASS1137','GASS24973','GASS1221','GASS41728','GASS38809','GASS41904','GASS20790','GASS39548','GASS11112','GASS1977')

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

count = 1
for (galName in galList){ # loop through galaxies
  if(verb){cat(paste("\n",galName," (",count,"/",length(galList),")\n",sep=""))}
  for (band in bandList){ # loop through bands
    if(verb){cat(paste("  ",band,"\n",sep=""))}
    for (nComps in compList){ # loop through number of components.
      if(verb){cat(paste("    ",nComps,"\n",sep=""))}
      ### INPUTS ### -> otherwise looped
      #galName = "GASS8634"
      #band = "r"
      #nComps = 2
      
      ### Get image file ###
      if(verb){cat("INFO: Retrieving data.\n")}
      imgFilename = paste(galName,"_",band,".fits",sep="")
      imgFile = paste(GALS_DIR,galName,band,imgFilename,sep='/')
      image = readFITS(imgFile)$imDat
      header = readFITS(imgFile)$hdr
      dims = dim(image)
      # Check for NaN padding:
      padded = is.element(NaN,image)
      
      ### Get PSF file ###
      psfFilename = paste(galName,"_",band,"_PSF.fits",sep="")
      psfFile = paste(GALS_DIR,galName,band,psfFilename,sep='/')
      psf = readFITS(psfFile)$imDat
      
      ### Create outputs folder ###
      if(verb){cat("INFO: Creating output directories\n")}
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
      image = image - SOFT_BIAS
      psf = psf - SOFT_BIAS
      
      
      
      ###########################################################
      #####  Make Segmentation map with ProFit (/ProFound)  #####
      ###########################################################
      if(verb){cat("INFO: Creating Segmentation image.\n")}
      # @Robin Cook:
      #segmentation = profitProFound(image, sigma=2, skycut=1.5, tolerance=4, size=11,
      #                              magzero=ZERO_POINT, gain=GAIN, header=header,
      #                              stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      
      segmentation = profitProFound(image, sigma=2.75, skycut=1.0, tolerance=4, size=21,
                                    magzero=ZERO_POINT, gain=GAIN, header=header,
                                    stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)
      # Find the main (central) source
      if(verb){cat("INFO: Finding central source.\n")}
      mainID = find_main(segmentation$segstats,dims) # The main source is the one with the smallest separation from the centre
      
      # Expand Segmentation image
      segmentationDilated = profitMakeSegimDilate(image, segmentation$segim, plot=TRUE, size= 41, expand=mainID)
      
      
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
      
      
      ################################################
      ###   Determine Sky statistics with ProFit   ###
      ################################################
      if(verb){cat("INFO: Calculating sky statistics.\n")}
      # Sky Estimate:
      #mask = profitProFound(image, sigma=3, skycut=0.5, tolerance=4, size=9, magzero=ZERO_POINT, gain=GAIN, header=header)
      mask = segmentationDilated$objects
      skyEst = profitSkyEst(image,mask=mask,radweight=1) # Structure containing sky, skyErr, skyRMS, + ...
      sky = skyEst$sky
      if (output && outputSkyStats){
        statsFilename = paste(galName,"_",band,"_SkyStats.png",sep='')
        png(paste(outputDir,statsFilename,sep='/'),width=1000,height=800,pointsize=16)
        par(mfrow=c(2,1), mar=c(3.5,3.5,1,2))
        skyEst = profitSkyEst(image,mask=mask,radweight=1,plot=T,xlab='Sky')  
        magplot(skyEst$radrun,xlab='radius (pixels)',ylab='Sky values',pch=16,grid=TRUE)
        abline(h=skyEst$sky,col='red',lty='dashed',lwd=2)
        dev.off()
      }
      
      # Subtract sky:
      image = image - sky
      
      
      ##########################################
      #####   Make Sigma map with ProFit   #####
      ##########################################
      if(verb){cat("INFO: Making sigma map.\n")}
      sigma = profitMakeSigma(image,sky=0.0,skyRMS=skyEst$skyRMS,gain=GAIN) # sky level is defined as 0.0 as sky has already been subtracted.
      
      
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
        reInits = c(inits$maj[mainID]*1.0,inits$maj[mainID]*3)
        nSerInits = c(4,1)
      } else {
        magInits = c(inits$mag[mainID])
        reInits = c(inits$maj[mainID]*1.0)
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
        isoConverge = TRUE # IF isoConverge=TRUE THEN use these as initial conditions. ELSE run LaplacesApproximation()
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
      
      # IF the result from 1D isophotal fitting did not converge THEN attempt a LaplaceApproximation() fit:
      if (improveInits==TRUE && isoConverge==FALSE) { # LaplaceApproximation LM fit
        if(verb){cat("INFO: Attempting to improve inital conditions with LaplaceApproximation()")}
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
        noise = matrix( rnorm(dims[1]*dims[2],mean=0.0,sd=skyEst$skyRMS), dims[1], dims[2]) 
        initModelFilename = paste(baseFilename,"_ModelImage.png",sep='')
        png(paste(outputDir,initModelFilename,sep='/'),width=1200,height=400,pointsize = 20)
        par(mfrow=c(1,3), mar=c(0.4,0.4,1,1))
        magimage(image,axes=F,bad=0); text(0.15*dim(image)[1],0.95*dim(image)[2],"Image",col='white',cex= 2)
        magimage(initImage+noise,axes=F,bad=0); 
        text(0.25*dim(image)[1],0.95*dim(image)[2],"Initial Model",col='white',cex= 2)
        text(0.715*dim(image)[1],0.05*dim(image)[2],expression(paste("N: ",mu, "=0, ", sigma,"=")),col='white',cex= 1.5)
        text(0.9*dim(image)[1],0.05*dim(image)[2],format(skyEst$skyRMS,digits=3),col='white',cex= 1.5)
        magimage(optimImage+noise,axes=F,bad=0);
        text(0.3*dim(image)[1],0.95*dim(image)[2],"Optimised Model",col='white',cex= 2)
        text(0.715*dim(image)[1],0.05*dim(image)[2],expression(paste("N: ",mu, "=0, ", sigma,"=")),col='white',cex= 1.5)
        text(0.9*dim(image)[1],0.05*dim(image)[2],format(skyEst$skyRMS,digits=3),col='white',cex= 1.5)
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
        cat(paste("Description:   ",description,"\n",sep=""))
        cat(paste("Date: ",Sys.time(),"\n\n",sep=""))
        cat(paste("Elapsed time: ",elapsedTime,"\n",sep=""))
        
        cat("\n\n>> Inputs:\n")
        cat(paste("Base directory: ",outputDir,"\n",sep=""))
        cat(paste("Image: ",imgFile,"\n",sep=""))
        cat(paste(" Dimensions: ",dims[1]," x ",dims[2],"\n",sep=""))
        cat(paste(" Padding: ",padded,"\n",sep=""))
        cat(paste("PSF: ",psfFile,"\n",sep=""))
        cat(paste("Zero point: ",ZERO_POINT,"\n",sep=""))
        cat(paste("Gain: ",GAIN,"\n",sep=""))
        cat(paste("Soft-bias: ",SOFT_BIAS,"\n",sep=""))
        cat(paste("Sky: ",skyEst$sky,"\n",sep=""))
        cat(paste("Sky RMS: ",skyEst$skyRMS,"\n",sep=""))
        
        
        cat("\n\n>> Segmentation:\n")
        cat(paste("Num. Objects: ",length(segmentation$segstats[[1]]),"\n",sep=""))
        cat(paste("MainID: ",mainID,"\n",sep=""))
        
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
        for (key in names(optimModel$sersic)){
          for (n in seq(1,nComps)){
            if (optimModel$sersic[[key]][n] %in% intervals$sersic[[key]][[n]] && tofit$sersic[[key]][n]){
              cat(paste(" - ",key,n,' = ',optimModel$sersic[[key]][n],'\n',sep=''))
            }
          }
        }
        if(segmentation$segstats[mainID,]$edge_frac < 0.8){cat(paste('\n - Segmentation boundary = ',segmentation$segstats[mainID,]$edge_frac,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$edge_excess > 1.0){cat(paste('\n - Segmentation edge excess = ',segmentation$segstats[mainID,]$edge_excess,'\n',sep=''))}
        if(segmentation$segstats[mainID,]$asymm > 0.2){cat(paste('\n - Segmentation asymmetry = ',segmentation$segstats[mainID,]$asymm,'\n',sep=''))}
        
        sink()
      }
      
      ### Append to optimization archive
      if (prefix == ""){
        #resultFilename = paste(HOME,"/robincook/Google Drive/PhD/Fitting/Results/Results","_",band,"_",nComps,"comp",".csv",sep = "")
        resultFilename = paste("/home/rcook/Documents/PhD/Fitting/Results/Results","_",band,"_",nComps,"comp",".csv",sep = "")
      }else{
        #resultFilename = paste(HOME,"/robincook/Google Drive/PhD/Fitting/Results/Results-",prefix,"_",band,"_",nComps,"comp",".csv",sep = "")
        resultFilename = paste("/home/rcook/Documents/PhD/Fitting/Results/Results-",prefix,"_",band,"_",nComps,"comp",".csv",sep = "")
      }
      append_output(resultFilename,galName,nComps,modellist$sersic,optimModel$sersic)
      
      ### Save workspace to file
      if(output && outputWorkspace){  
        workspaceFilename = paste(baseFilename,"_WorkSpace.RData",sep='')
        save.image(paste(outputDir,workspaceFilename,sep='/'))
      }
      
    
    } # END nComps loop
  } # END band loop
  count = count + 1
} # END Galaxy loop

