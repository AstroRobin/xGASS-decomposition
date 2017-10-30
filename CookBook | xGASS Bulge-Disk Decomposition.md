
<h1><span style="color:mediumblue">Bulge/Disk Decomposition</span> <span style="color:navy">Cook</span><span style="color:mediumblue">book</span></h1><br>
<h3>Author: Robin H. W. <b>Cook</b></h3>
<h3>Date: 22 September, 2017</h3>

<h3><span style="color:navy"><em>Preamble</em></span></h3>
<span style="color:DimGray">
<em>In this <b>Cook</b>book, I will present the processes used to model the light profile of galaxies within the <i>extended GALEX Arecibo SDSS Survey</i> (xGASS) sample. We use various likelihood optimisation and image segmentation functions in [ProFit](https://github.com/ICRAR/ProFit "ProFit")/[ProFound](https://github.com/asgr/ProFound "ProFound") (developed by Aaron Robotham) to accurately disentangle the contributions of the compact inner region and extended outer region to the total light of a galaxy.
Most of the tools and scripts used to process the SDSS data and perform the bulge/disk decomposition can be found here:</em> https://github.com/AstroRobin/xGASS-decomposition
</span>

<h2><span style="color:navy"><u>Processing SDSS Data:</u></span></h2>

<h3><span style="color:navy">Data</span></h3>
The photometric data used for the bulge-disk decomposition of the xGASS galaxies comes from the Sloan Digital Sky Survey (SDSS) data release 7 (DR7) Data Access Server (DAS; http://classic.sdss.org/dr7/access/). This data release was chosen simply for the ease of access to the uncalibrated (in units of "counts") data products. One could also utilise later data releases if need be, however the image would need to be converted from the native SDSS flux units of ["nanomaggies" to photon counts](http://classic.sdss.org/dr5/algorithms/fluxcal.html).

The basic imaging products needed for the model fitting are the fpC-.fit (or drC-.fit for fits header supplemented images) files: these are the uncalibrated, corrected frames (bias subtracted, flat fielded, and with bad pixels replaced by interpolated values), as well as the psField-.fit files which contain the point spread function (PSF) fit for all positions in a particular field.

<h3><span style="color:navy">Fields (images)</span></h3>
The imaging files are retrievable from SDSS as individual frames. Here, fields are the individual data streams from a single CCD (one filter) in a scanline for a particular observing run. They measure 2048 x 1489 pixels (~13.5' x 9.8'), overlapping by 10% (128 row pixels). The frames of each of the 5 filters (u,g,r,i,z) for the same part of the sky are called a field. More information about the data products and terminology can be found here: http://www.sdss3.org/dr8/glossary.php.

The fields are uncalibrated, meaning that they are in units of Data Numbers (The SDSS equivalent of counts); ProFit can understand these units given that the Gain and Magnitude zeropoint are also provided such that it can convert easily between counts and fluxes/magnitudes. These numbers can be found in the header of DR7 drC-.fit files as well as in tsField-.fit and fpAtlas-.fit files. More information about converting between photometric units in SDSS can be found here: http://classic.sdss.org/dr5/algorithms/fluxcal.html


<h3><span style="color:navy">Cropping images</span></h3>
The SDSS frame images are far too large for use in ProFit and are most likely not centred upon the target of choice. For this reason, the user should crop the SDSS frames with their favourite cutout tool. For best use in ProFit One should ensure that the image is centred upon the galaxy and that it is large enough to include the diffuse outer regions of the object as well as including enough sky pixels to make accurate measurements of the background sky statistics.
I use the <span style="color:purple">astropy.Cutout2D()</span> function in Python, but might like to use <span style="color:purple">magcutoutWCS()</span> in the magicaxis package to maintain the workflow within R. Below, I will show an example of the cutout process for [GASS111029](http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id=588023669705146551) (RA = 175.56089, Dec = 20.09783) in a <b>python</b>.


```python
from astropy.io import fits as pyfits # FITS input/output
from astropy import wcs # World Coordinate System
from astropy.nddata import Cutout2D # Image cutout function
import warnings; warnings.filterwarnings('ignore') # Shush silly Astropy warnings

import matplotlib.pyplot as plt # plotting
import matplotlib.colors as mcol # colourbar scaling
```


```python
# Get the image data
fitsFile = pyfits.open("Data/drC-005225-r4-0076.fits") # Open the fits file
image = fitsFile[0].data - 1000 # There is a 1000 data number pedestal added to uncalibrated SDSS images
header = fitsFile[0].header # image header
w = wcs.WCS(header) # define the "World Coordinate System"
```


```python
# plot the image
print("\nSDSS field image:")

# plotted logarithmically for more contrast
plt.matshow(image,cmap='hot',norm=mcol.LogNorm())
plt.colorbar(orientation='vertical',shrink=0.8)
    
plt.show()
```

    
    SDSS field image:



![png](output_3_1.png)



```python


# define a position
ra = 175.56089
dec = 20.09783

# Get the x,y pixel location of the target
x, y = w.wcs_world2pix(ra,dec,0)
position = (x,y)
print("Pixel coordinates:\n  x = {0}\n  y = {1}".format(x,y))

# Define the cutout dimensions (2' x 2'):
cutSize = 2.0 * 60.0 # arcsec
pixDim = (cutSize/0.396, cutSize/0.396) # SDSS pixel scale: 0.396"/pixel

cutout = Cutout2D(image, position, pixDim, wcs = w)

print("\nGASS111029 - 1'x'1 cutout:")
plt.matshow(cutout.data,cmap='hot',norm=mcol.LogNorm())
plt.colorbar(orientation='vertical',shrink=0.8)

plt.show()
```

    Pixel coordinates:
      x = 1321.6083489900204
      y = 781.9720357360613
    
    GASS111029 - 1'x'1 cutout:



![png](output_4_1.png)



```python
# Save image to file:
newHeader = header
newHeader['naxis1'] = pixDim[0]
newHeader['naxis2'] = pixDim[1]

newhdu = pyfits.PrimaryHDU(cutout.data,header=newHeader)
hdulist = pyfits.HDUList([newhdu])
hdulist.writeto("Data/GASS111029_image.fits",overwrite=True)
```

A nice python wrapper for this tool can be found on the [github repo](https://github.com/AstroRobin/xGASS-decomposition).

<h3><span style="color:navy">Point Spread Functions (PSFs)</span></h3>
The psField-.fit files contain the information of the point-spread-function for a corresponding field (u,g,r,i,z) in an imaging run. However the format of these files is such that one needs to run a specialised tool from SDSS to extract the PSF fit at a particular location and filter within the SDSS field. This tool is part of a code to read atlas files found here: http://www.sdss.org/dr12/imaging/images/.
Once the tool has been installed and compiled, the user may extract a PSF from the psField-.fit files with a command like the following:

<code>
% read_PSF -h
Usage: read_PSF [options] input-file hdu row col output-file
Your options are:
    -? This message
    -h This message
    -i Print an ID string and exit
    -v Turn up verbosity (repeat flag for more chatter)
</code>

Here, <em>hdu</em> refers to the individual filters (1:5 for u,g,r,i,z) and row, col refer to the pixel coordinates from which the user wishes to extract the PSF. Bear in mind that the PSF will vary across any field and so it is important to ensure that the pixel locations correspond to the target object. Below is an example of a PSF around a particular galaxy in the r-band filter.

So, for example with GASS111029, we would extract the PSF image from the psField-.fit file like so:

<code>
% read_PSF Data/psField-005225-4-0076.fits 3 1322 782 Data/GASS111029_r_PSF.fits
</code>

where the 'row' and 'col' values were taken from the WCS conversion from sky coordinates to image pixels in the code block above.


```python
# Get the PSF image
PSFFile = pyfits.open("Data/GASS111029_r_PSF.fits") # Open the fits file
psf = PSFFile[0].data # There is a 1000 data number pedestal added to uncalibrated SDSS images

print("\nGASS111029 - Point-Spread Function:")
plt.matshow(psf,cmap='hot',norm=mcol.LogNorm(vmin=1000,vmax=psf.max())) # vmin = 1000 to account for the SOFTBIAS=1000
cbar = plt.colorbar(orientation='vertical',shrink=0.8)

plt.show()
```

    
    GASS111029 - Point-Spread Function:



![png](output_8_1.png)


<em>Note: Both uncalibrated drC-*.fit and psField-*.fit files have an added SOFT_BIAS of 1000 Data Numbers added to each pixel as a pedestal that must be subtracted prior to fitting.</em>

<h2><span style="color:navy"><u>Galaxy Profile Modelling:</u></span></h2>

<h3><span style="color:navy">Setting up R</span></h3>
<p>The main packages used to perform the galaxy modelling in R are: </p>
<ul> 
<li> <em>ProFit</em> - 2D Bayesian galaxy fitting tool
<li> <em>ProFound</em> - Source detection and image segmentation tool.
<li> <em>LaplacesDemon</em> - MCMC optimisation package.
<li> et al. (FITSio, magicaxis, etc.)
</ul>

<em>LaplacesDemon</em> is a package containing many MCMC optimsation routines which - given data, a specified model, and initial values - will maximise the logarithm of the marginal posterior distributions in the model parameters. <em>LaplacesDemon</em> is available of [CRAN](https://cran.r-project.org/web/packages/LaplacesDemon/index.html) and some a useful tutorial/beginners guide can be found here: https://cran.r-project.org/web/packages/LaplacesDemon/vignettes/LaplacesDemonTutorial.pdf.

<em>ProFit</em> and <em>ProFound</em> are both developed by Dr. Aaron Robowtham; they are available of CRAN but it is most useful to source them directly from their repsective github repositiories:
<ul>
<li> <em>ProFit</em>: https://github.com/ICRAR/ProFit
<li> <em>ProFound</em>: https://github.com/asgr/ProFound
</ul>

Below, I will show how these packages can be imported in R:



```python
### Setting up Jupyter to run R commands with rpy2
# Anything with a proceeding '%R' is an R command and should not be included if running within an R script.

import rpy2
%load_ext rpy2.ipython
```

    The rpy2.ipython extension is already loaded. To reload it, use:
      %reload_ext rpy2.ipython



```python
# ProFit/ProFound packages
%R library(devtools) # Needed for install_github() function which allows for packages to be sourced from GitHub repos
%R install_github("ICRAR/ProFit")
%R install_github("asgr/ProFound")

# MCMC optimisation routines
%R library(LaplacesDemon)

# Extras
%R library(FITSio) # .FITS file input/output
%R library(EBImage) # Image processing package
%R library(magicaxis) # "Magically Pretty Plots" ~ ASGR
```

<h3><span style="color:navy">Data Setup</span></h3>
Prior to running ProFit on these uncalibrated SDSS images, one must first subtract the SOFT BIAS (1000 Data Numbers) pedestal from both the image and PSF. Here is some <b>R code</b> which shows the importing and processing of the image data:

<code>
  library(ProFit)
  library(FITSio) \# Read/write FITS files

  \# Get image file
  imgFile = "GASS11029_r.fits"
  image = readFITS(imgFile)$imDat
  header = readFITS(imgFile)$hdr
  dims = dim(image) \# image dimensions

  \# Get PSF file
  psfFile = "GASS11029_r_PSF.fits"
  psf = readFITS(psfFile)$imDat

  \# Get information from FITS header
  SOFT_BIAS = as.numeric(header[which(header=="SOFTBIAS")+1]) \# SOFT_BIAS = 1000 DNs pedestal
  ZERO_POINT = as.numeric(header[which(header=="ZP")+1]) \# the zeropoint magnitude
  GAIN = as.numeric(header[which(header=="GAIN")+1]) \# Conversion between counts & photo-electrons

  \# Subtract SOFT_BIAS from image and PSF
  image = image - SOFT_BIAS
  psf = psf - SOFT_BIAS
  
</code>

The ZERO_POINT & GAIN will be used later for converting between counts and magnitudes. ProFit can handle images which are given in counts (here, Data Numbers by SDSS terminology).

If you are using SDSS data release 8+, you may find that your images are already calibrated and in units of SDSS specific nano-maggies. It is best that you convert your images into the appropriate units of counts (DNs) for ProFit (http://www.sdss.org/dr12/algorithms/fluxcal/).

<h3><span style="color:navy">Source Extraction and Image Segmentation</span></h3>
The ProFound package is used to perform the source extraction o the image which is then used to determine an appropriate segmentation region within which ProFit will model the target source. This is done in two stages:

- Firstly, I perform a source extraction on the image to find all distinct objects. Tuning certain parameters are necessary in order to robust segregate individual sources whilst also maintaining as much of the flux from a flocculent galaxy, for example.
- Secondly, I expand only the target object (identified as being within the centre of the image) such that it is pushed well within the sky of the image. If this second step is not performed, often the segmentation map does not anchor to the noise floor and so it is more difficult to constrain the light profile of the galaxy (particularly in the diffuse outer regions).

The <b>R code</b> code used for these steps is given below:

<code>
\# Get initial object list:
segmentation = profitProFound(image, sigma=2.75, skycut=1.0, tolerance=4, size=21,                             magzero=ZERO_POINT, gain=GAIN, header=header, stats=TRUE, rotstats=TRUE, boundstats=TRUE, plot=TRUE)

\# Get the target object (central-most)
mainID = find_main(segmentation$segstats,dims)

\# Expand Segmentation image
segmentationDilated = profitMakeSegimDilate(image, segmentation$segim, plot=TRUE, size= 41, expand=mainID)
</code>


<h3><span style="color:navy">Modelling the Surface Brightness Profiles of Galaxies</span></h3>

The goal of bulge-disk decomposition is to accurately separate the contributions of the bulge and disk to that of a galaxy's total flux. At minimum, this requires fitting a combination of two light profiles to the galaxy: one to for the spherically symmetrical and centrally concentrated bulge; the other for the extended and highly-flattened disk. However, it has been noted by many authors ([Laurikainen et al. 2004](http://adsabs.harvard.edu/abs/2004MNRAS.355.1251L), [2005](http://adsabs.harvard.edu/abs/2005MNRAS.362.1319L); [Gadotti 2008](http://adsabs.harvard.edu/abs/2008MNRAS.384..420G); [Kim et al. 2014](http://iopscience.iop.org/article/10.1088/0004-637X/782/2/64/meta)) that excluding secondary features (bars, lenses, rings, disk breaks, and spiral arms) from the models of galaxies where they are present, may introduce major uncertainties in the parameters of the derived bulge component ([Gao & Ho 2017](http://adsabs.harvard.edu/abs/2017ApJ...845..114G)). Most galaxy model-fitting tools typically provide researchers with a large toolbox of analytic functions (e.g. Sérsic, modified Ferrer, Core-Sersic, broken-exponential, Moffat, etc.; see [Robotham et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.466.1513R) for a description of all those available in <em>ProFit</em>) from which models can be constructed in which an arbirary combination of components may be used and for which the parameters of each may be free, fixed or constrained relative to one another (i.e. centres that are fixed upon one another).

For ideal galaxies, as the density of stars, and therefore light, decreases with increasing radius from the centre of a galaxy, the overall surface brightness of the galaxy also falls off. The manner by which this occurs is encapsulated in the surface brightness profile. In the case of a stellar disk, the surface brightness profile declines exponentially. In fact, this exponentially decaying profile is nothing more than a special case of the more generalised Sérsic profile [Sersic 1963](http://adsabs.harvard.edu/abs/1963BAAA....6...41S) in which the light intensity profile is given by the following equation:

$$I(r)=I_{e}e^{-b_{n}\left[\left(\frac{r}{r_{e}}\right)^{1/n}-1\right]}$$

where $I_{e}$ is the intensity of the profile at the radius $R_{e}$ (the radius within which half of the flux is contained), n is the Sérsic index. The Sérsic index has special cases of $n=0.5$ for a Gaussian, $n=1$ for an exponential and $n=4$ for the de Vaucouleurs profile (a good description of giant elliptical galaxies and central bulges). The $b_{n}$ term ensures the correct integration properties at $R_{e}$; explicitly, it is derived as the quantile at which the Gamma probability distribution function integrates to exactly 0.5 for a shape parameter of $2n$, i.e. $\Gamma(2n) = 2\gamma\,(2n,\,b_{n})$. For surface brightness profiles considered in two dimensions, i.e. as one would observe in the plane of the sky, the isophotal contours take a more complex form:

$$ R_{t} = \sqrt{(x-x_{0})^{2}+(y-y_{0})^{2}} $$


where $R_t$ is now the modified radius in two-dimensions; in the one-dimensional case, this simply meant the radius from the centre of the profile. Here, $x$, $y$ denote the position at which the light profile is being evaluated and $x_{0}$, $y_{0}$ indicate the the central position of the profile. As it is possible for a galaxy to be positioned in any orientation in the plane of the sky, a term describing the position angle ($\theta$) is introduced. In <em>ProFit</em>, $\theta$ is defined as being 0$^{\rm{\circ}}$ when the semi-major axis is aligned vertically and increases positively as the galaxy is rotated in a counter-clockwise direction. The modified radius then becomes:

$$ R_{m} = \sqrt{\left(R_{t}\sin{(\theta_{t})}A_{rat}\right)^{2} + \left(R_{t}\cos{(\theta_{t})}\right)^{2}} $$

where $R_{t}$ is the two-dimensional radius in the equation prior, $A_{rat}$ is the ratio between the semi-minor and semi-major axes and $\theta_{t}$ is given by $\theta_{t}=\tan^{-1}\left(\frac{x-x_{0}}{y-y_{0}}\right) + \theta$.
Finally, one may allow for the property of 'boxiness' which skews the unit circle either giving a diamond or box-like appearance. This alters $R_{m}$ in the above equation in the following way:

$$ R_{m} = \left[(R_{t}\sin{(\theta_{t})}A_{rat})^{2+B} + (R_{t}\cos{(\theta_{t})})^{2+B}\right]^{\frac{1}{2+B}} $$

here, $B$ is considered to be the 'boxiness' parameter and should vary between -2 and 2 for most practical uses.\\
Thus, a single Sérsic profile can be described in full by specifying <b>eight parameters</b>: two positional parameters $x_{0}$ and $y_{0}$; axial ratio $A_{rat}$; position angle $\theta$; boxiness $B$; some measure of the brightness (typically given as the magnitude $m$); the half-light radius $r_{e}$; and the Sérsic index $n$. I have provided this description of the Sérsic profile such that the reader may appreciate the level of complexity that can exist in modelling the light profiles of galaxies (even in the idealised two-component system) and to highlight the demand for effecient and robust optimisation tools, such as those available through <em>ProFit</em>.
