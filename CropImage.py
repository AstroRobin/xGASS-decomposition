"""
Crops an input FITS image to a specified pixel dimension about a position (ra, dec).

Author: Robin Cook
Date: 24/03/17

TO-DO: 
 - Add combine feature.
"""

# Astropy imports
from astropy.io import fits as pyfits # FITS input/output
from astropy import wcs # World Coordinate System
from astropy.nddata import Cutout2D # Image cutout function
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5 #, Ecliptic?
import astropy.units as u

# Shush Astropy warnings
import warnings; warnings.filterwarnings('ignore')

# Standard imports
import os
import subprocess
import sys #?
import math
from math import nan
import numpy as np

import matplotlib.pyplot as plt # Plotting
from argparse import ArgumentParser # Parsing Arguments from sys.argv

# Define constants
SDSS_PSCALE = 1.1e-4
HOME = str(subprocess.check_output("pwd")).split('/')[1]

def main():

	# Create argument parser
	parser = ArgumentParser(description="Tool for cropping astronomical images.")
	IOGroup = parser.add_argument_group('Input/Output',"Arguments for the input and output FITS images.")
	posGroup = parser.add_argument_group('Positional',"Arguments for the position of the cropped image.")
	controlGroup = parser.add_argument_group('Control',"Arguments for controlling the form of the output image.")
	miscGroup = parser.add_argument_group('Miscellaneous',"Miscellaneous arguments.")

	IOGroup.add_argument('-i','--input',
						action='store',dest='input',type=str,nargs='+',default=None,
						help="The path to the input FITS image file.",metavar="INPUT_IMAGE")
	IOGroup.add_argument('-o','--output',
						action='store',dest='output',type=str,default=None,
						help="The path of the output FITS image file.",metavar="OUTPUT_IMAGE")
	IOGroup.add_argument('-f','--file',
						action='store',dest='frameFile',type=str,default=None,
						help="The path to a file containing a list of inputs (format: 'ID,ra,dec,band,run,camcol,rerun,field')",metavar="")
	IOGroup.add_argument('-b','--band',
						action='store',dest='band',type=int,nargs='*',default=None,
						help="If '--file' specified, indicates which frequency bands to crop.",metavar="BANDS")
	IOGroup.add_argument('-g','--galaxy',
						action='store',dest='galaxy',type=str,default=None,
						help="The name of the galaxy. For use in '--log' only.",metavar="GALAXY")
	IOGroup.add_argument('-l','--log',
						action='store',dest='logFile',type=str,default=None,
						help="Specifies the location of a log file to append errors.",metavar="LOGFILE")
	posGroup.add_argument('-r','--ra',
						action='store',dest='ra',type=float,default=None,
						help="The Right Ascension for the centre of the output image (Default: ra = central)",metavar="RA")
	posGroup.add_argument('-d','--dec',
						action='store',dest='dec',type=float,default=None,
						help="The Declination for the centre of the output image (Default: dec = central)",metavar="DEC")
	posGroup.add_argument('-x','--xpos',
						action='store',dest='x',type=float,default=None,
						help="The x position in pixels for the centre of the output image (Default: x = central)",metavar="XPOS")
	posGroup.add_argument('-y','--ypos',
						action='store',dest='y',type=float,default=None,
						help="The y position in pixels for the centre of the output image (Default: y = central)",metavar="YPOS")
	posGroup.add_argument('-c','--coordinate',
						action='store',dest='frame',type=str,default='ICRS',
						help="The reference frame used for the sky coordinates (Default: ICRS)",metavar="FRAME")
	posGroup.add_argument('-s','--size',
						action='store',dest='size',type=str,default=None,
						help="The dimensions of the output FITS image (Default: xDim = 100 pixels; yDim = 100 pixels)",metavar="SIZE")
	posGroup.add_argument('-u','--units',
						action='store',dest='units',type=str,default='p',
						help="The units of the input dimensions ([ (\033[4md\033[0megrees) | (arc\033[4mm\033[0minutes) | (arc\033[4ms\033[0meconds) | (\033[4mp\033[0mixels) ]) (Default: Pixels)",metavar="UNITS")
	controlGroup.add_argument('-P','--pad',
						action='store',dest='pad',choices=['None','NaN','0','-1'],	default='None',
						help="Whether to pad any whitespace that may exist if the cropped image exceeds the edge of the frame.")
	controlGroup.add_argument('-C','--combine',
						action='store_true',dest='combine',default=False,
						help="Whether to combine multiple frames before cropping the image (Requires multiple arguments in '--image').")
	miscGroup.add_argument('-p','--plot',
						action='store_true',dest='plot',default=False,
						help="Plot cutout.")
	miscGroup.add_argument('-v','--verbose',
						action='store_true',dest='verbose',default=False,
						help="Print messages to stdout.")
	
	# Parse Arguments
	args = parser.parse_args()
	verbose = args.verbose

	# Check Input file(s)
	if (args.logFile != None):	
		log = open(args.logFile,'a')
		name = args.galaxy if args.galaxy else args.input[0].split('/')[-1]
	else:
		log = None

	if (args.input == None): print("ERROR: No input FITS file given.\n   -- ABORTING --"); exit()
	if (args.combine):
		inputPaths = []
		for ii in range(len(args.input)):
			if (os.path.exists(args.input[ii])): inputPaths.append(args.input[ii])
			else: print("ERROR: input FITS file \"{0}\" does not exist.".format(args.input[ii]))
		if (len(inputPaths) == 0): print("\nERROR: No valid input images given!\n  -- ABORTING --\n"); exit()
	else:
		if (len(args.input) > 1): print("\nWARNING: Multiple images given. Only using the first: \"{file}\"\n".format(file=args.input[0]))
		if (os.path.exists(args.input[0])): inputPaths = [args.input[0]]
		else:
			print("ERROR: input FITS file \"{0}\" does not exist.\n   -- ABORTING --".format(args.input[0]));
			if log: log.write("** {name:<10}: File \"{file}\" does not exist! **\n".format(name=name,file=args.input[0]))
			exit()
	
	# Check Output file
	if (args.output == None): outputPath = inputPath.replace('.fit','_crop.fit')
	else: outputPath = args.output

	# Open FITS file
	if (args.combine):
		# Swarp First
		# Then Open SWarped image
		print("\nERROR: --combine not implemented!\n  -- ABORTING --\n"); exit()
		pass
	else:
		fitsFile = pyfits.open(inputPaths[0])
		hdr = fitsFile[0].header
		w = wcs.WCS(hdr)

		if (verbose>1): print(w)

		imgDims = np.shape(fitsFile[0].data)

	# Get position
	if ((args.ra and args.dec) and (args.x and args.y)):
		print("\nWARNING: Both --ra/--dec and --xpos/--ypos have been given. Keeping --xpos/--ypos")
		args.ra = None; args.dec = None

	if not ((args.ra and args.dec) or (args.x and args.y)):
		print("\nWARNING: No position arguments have been given. Using centre value from FITS header.")
		ra = hdr['CRVAL1']
		dec = hdr['CRVAL2']

	if (args.ra and args.dec): 
		ra = args.ra # Check if outside of image
		dec = args.dec # Check if outside of image

		# Get Reference frame
		if (args.frame.lower() not in ["icrs","galactic","ecliptic","fk4","fk5"]):
			print("WARNING: Reference frame \"{0}\" does not exist".format(args.frame))
		else:
			frame = args.frame.lower()

		# Convert RA/Dec to x/y position
		#position = SkyCoord(ra, dec, frame=frame, unit='deg')

		if (hdr['CTYPE1'] == 'DEC--TAN' and hdr['CTYPE2'] == 'RA---TAN'):
			x, y = w.wcs_world2pix(dec,ra,0)
		elif (hdr['CTYPE1'] == 'RA---TAN' and hdr['CTYPE2'] == 'DEC--TAN'):
			x, y = w.wcs_world2pix(ra,dec,0)
		else:
			print("\nERROR: No valid coordinates in .fits header (check CTYPE(s) are set to 'RA---TAN' & 'DEC--TAN')\n  -- ABORTING --\n"); exit()

		position = (x,y)

	if (args.x and args.y):
		x = args.x
		y = args.y

		position = (x,y)
		frame = "None"
	
	# Get image dimensions
	if (args.size == None): xDim = 100; yDim = 100
	else:
		size = args.size.split(',')
		if (args.units == "s" or "arcsec" in args.units.lower()): pixelScale = 1.0/(3600*SDSS_PSCALE); units = "\""
		elif (args.units.lower() == "m" or "arcmin" in args.units.lower()): pixelScale = 1.0/(60*SDSS_PSCALE); units = "\'"
		elif (args.units.lower() == "d" or "deg" in args.units.lower()): pixelScale = 1.0/(SDSS_PSCALE); units = " deg"
		elif (args.units.lower() == "p" or "pixel" in args.units.lower()): pixelScale = 1.0; units = " pixels"
		else: print("WARNING: units \"{0}\" is not appropriate. Assuming units as pixels (\"p\")".format(args.units))

		size = [float(val) for val in size]
		if (len(size) == 1): 
			pixelDims = (math.ceil(size[0]*pixelScale),math.ceil(size[0]*pixelScale))
			dims = ("{0}{1}".format(size[0],units),"{0}{1}".format(size[0],units))
		else:
			pixelDims = (math.ceil(size[0]*pixelScale),math.ceil(size[1]*pixelScale))
			dims = ("{0}{1}".format(size[0],units),"{0}{1}".format(size[1],units))

	# Validate Padding:
	if   (args.pad == 'None'): padding = None
	elif (args.pad.lower() == 'nan'): padding = nan
	elif (args.pad in ['0','0.0','0.']): padding = 0
	elif (args.pad in ['-1','-1.0','-1.']): padding = -1
	else: print("\nWARNING: Padding value \"{pad}\" is not valid. Using \"NaNs\" instead.\n".format(pad=args.pad)); padding = nan

	if (verbose):
		if (frame.lower() == 'galactic'):
			print("\n  -- CropImage.py -- \n\nParameters:\n > Image I/O:\n   - input image: \"{0}\"\n   - output image: \"{1}\"\n > Position:\n   - Long = {2.l:.4f}\n   - Lat = {2.b:.4f}\n > Image Dimensions:\n   - xDim = {3[0]} ({4[0]} pixels)\n   - yDim = {3[1]} ({4[1]} pixels)\n".format(','.join(inputPaths),outputPath,position,dims,pixelDims,units))
		else:
			print("\n  -- CropImage.py -- \n\nParameters:\n > Image I/O:\n   - input image: \"{0}\"\n   - output image: \"{1}\"\n > Position:\n   - {2}\n > Image Dimensions:\n   - xDim = {3[0]} ({4[0]} pixels)\n   - yDim = {3[1]} ({4[1]} pixels)\n".format(','.join(inputPaths),outputPath,position,dims,pixelDims,units))

		print("Pixel dims: {0}\nPositions\n x = {1}\n y = {2}".format(pixelDims,x,y))

	if (verbose and args.pad != None): print("\nINFO: Using \"{0}\" as padding for overlapping edges.".format(padding))
	try:
		cutout = Cutout2D(fitsFile[0].data.astype(float), position, pixelDims, wcs = w, mode = 'trim' if padding == None else 'partial', fill_value = padding)
	except: # Coordinates do not line up with frame
		if log: log.write("** {name:<10}: Coordinates given ({x:.3f} , {y:.3f}) do not overlap with image (dims = {xdim} , {ydim}) **\n".format(name=name,x=float(x),y=float(y),xdim=hdr['NAXIS1'],ydim=hdr['NAXIS2']))
		print("\nERROR: Coordinates given ({x:.3f},{y:.3f}) do not overlap with image (dims = {xdim},{ydim})\n  -- ABORTING --\n".format(x=float(x),y=float(y),xdim=hdr['NAXIS1'],ydim=hdr['NAXIS2'])); exit()

	
	if(verbose): print("\nINFO Cutout array:\n",cutout.data)
	
	if (args.plot):
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.imshow(cutout.data,origin='lower')
		plt.show()

	outDims = np.array(cutout.data).shape
	if (outDims[0] != pixelDims[0] or outDims[1] != pixelDims[1]):
		if (log): log.write("** {name:<10}: Pixel dimensions trimmed at edge ({x},{y}) **\n".format(name=name,x=outDims[0],y=outDims[1]))
		print("\n *** WARNING: Pixel dimensions trimmed at edge ({x},{y}) ***\n".format(x=outDims[0],y=outDims[1]))
	
	hdrCutout = hdr
	hdrCutout['naxis1'] = pixelDims[0]
	hdrCutout['naxis2'] = pixelDims[1]
	
	# Convert reference pixel from origin to that of the cutout
	oldRefPixels = (hdr['CRPIX1'],hdr['CRPIX2'])
	newRefPixels = cutout.to_cutout_position(oldRefPixels)

	#print oldRefPixels
	#print newRefPixels

	hdrCutout['CRPIX1'] = newRefPixels[0]
	hdrCutout['CRPIX2'] = newRefPixels[1]

	newhdu = pyfits.PrimaryHDU(cutout.data,header=hdrCutout)
	hdulist = pyfits.HDUList([newhdu])
	hdulist.writeto(outputPath,overwrite=True)

if __name__ == "__main__": main()