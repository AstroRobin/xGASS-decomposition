"""
A (wrapper) script for extracting the PSF from psField-***.fit SDSS images.

Author: Robin H. Cook
Date: 05/04/17
"""

# Standard imports
import os
import subprocess
from argparse import ArgumentParser # Argument/Option parsing
from astropy.io import fits as pyfits
from astropy import wcs # World Coordinate System

# Shush Astropy warnings
import warnings; warnings.filterwarnings('ignore')

import numpy as np

# Constants
HOME = str(subprocess.check_output("pwd")).split('/')[1]
PSF_PATH = "\"/{home}/robincook/Google Drive/PhD/SDSS/Programs/readAtlasImages-v5_4_11/read_PSF\"".format(home=HOME)

def find_edges(array):
	"""
	Find the regions of trim outside of the PSF image.

	<param: array [2D array (float)]> - The input image array

	<return: bounds [list ([x1,y1],[x2,y2])]> - the x/y positions of the upper-left and lower-right bounds fopr the trim region.
	"""

	aveCol = np.mean(array,axis=0)
	aveRow = np.mean(array,axis=1)

	# Get image dimensions
	dims = np.shape(array)
	rows = dims[0]; cols = dims[1]

	# forwards through columns
	for ii in range(0,cols,+1):
		if (aveCol[ii+1] != aveCol[ii]):
			x1 = ii+1
			break

	# backwards through columns
	for ii in range(cols-1,-1,-1): # from last element to zero-th element backwards
		if (aveCol[ii-1] != aveCol[ii]):
			x2 = ii-1
			break

	# forwards through rows
	for ii in range(0,rows,+1):
		if (aveRow[ii+1] != aveRow[ii]):
			y1 = ii+1
			break

	# backwards through rows
	for ii in range(rows-1,-1,-1): # from last element to zero-th element backwards
		if (aveRow[ii-1] != aveRow[ii]):
			y2 = ii-1
			break

	return ([x1,y1],[x2,y2])

def main():

	parser = ArgumentParser()

	parser.add_argument('input',help="The input file path.",default=None)
	parser.add_argument('band',help="The frequency band of the image.",default=None)
	parser.add_argument('-r','--ra',
						action='store',dest='ra',type=float,default=None,
						help="The Right Ascension of the region of interest (Requires: -i/--image).",metavar="RA")
	parser.add_argument('-d','--dec',
						action='store',dest='dec',type=float,default=None,
						help="The Declination of the region of interest (Requires: -i/--image).",metavar="DEC")
	parser.add_argument('-i','--image',
						action='store',dest='image',type=str,default=None,
						help="If --ra/--dec specified, use this image to find the corresponding pixel coordinates.",metavar="IMAGE")
	parser.add_argument('-x','--col',
						action='store',dest='x',type=float,default=None,
						help="The column at which to extract the PSF.",metavar="XPOS")
	parser.add_argument('-y','--row',
						action='store',dest='y',type=float,default=None,
						help="The row at which to extract the PSF.",metavar="XPOS")
	parser.add_argument('-o','--output',
						action='store',dest='output',type=str,default=None,
						help="The output FITS file path (Default = <inputFile>_PSF.fits)",metavar="OUTPUT")
	parser.add_argument('-t','--notrim',
						action='store_false',dest='trim',default=True,
						help="Whether to trim the excess padding from the SDSS output PSF (Default = True).")
	parser.add_argument('-v',
						action='count',dest='verbosity',default=0,
						help="The level of verbosity to print to stdout.")

	args = parser.parse_args() # Parse arguments

	# Validate verbosity
	global vrb
	if (args.verbosity < 0): vrb = 0
	elif (args.verbosity > 2): vrb = 2
	else: vrb = args.verbosity

	# Get input filename
	if (os.path.exists(args.input) and "psField" in args.input): inputFile = args.input
	else: print("\nERROR: file \"{0}\" is not a valid \"psField-***.fit\" file!\n -- ABORTING --".format(args.input)); exit()

	# Get file path information
	filename = inputFile.split('/')[-1]
	directory = '/'.join(filename.split('/')[0:-1])
	if directory == "": directory = "."

	# Get output filename
	if (args.output == None): outputFile = inputFile.replace('psField','PSF')
	else:
		if (len(args.output.split('/')) > 1 and HOME in args.output.split('/')): # if full path given
			outputFile = args.output
		else: # absolute path reference
			outputFile = "{0}/{1}".format(directory,args.output)


	# Print input parameters to stdout
	if(vrb>0):print("\nINFO: Input file  = {0}\nINFO: Output file = {1}".format(inputFile,outputFile))

	# Get positions (Ensure pixels are within image)
	#if (args.col >= 0 and args.col < hdr['NAXIS1']): col = args.col
	#else: print("ERROR: Invalid column position."); exit()
	#if (args.row >= 0 and args.row < hdr['NAXIS2']): row = args.row
	#else: print("ERROR: Invalid row position."); exit()

	if (args.x != None and args.y != None):
		if (args.x >= 0): col = args.x
		if (args.y >= 0): row = args.y

	elif(args.ra != None and args.dec != None):
		if (args.image == None): print("\nERROR: argument \"-i\"/\"--image\" required with \"--ra\" and \"--dec\".\n  -- ABORTING --"); exit()
		
		# Open .FITS file
		if (os.path.exists(args.image)):
			fitsFile = pyfits.open(args.image)
			img = np.array(fitsFile[0].data)
			hdr = fitsFile[0].header
			dims = np.shape(img)
		else:
			print("\nERROR: file \"{0}\" does not exist\n  -- ABORTING --".format(args.image)); exit()

		# Define World-Coordinate-System (wcs):
		w = wcs.WCS(hdr)

		if (hdr['CTYPE1'] == 'DEC--TAN' and hdr['CTYPE2'] == 'RA---TAN'):
			col, row = w.wcs_world2pix(args.dec,args.ra,1)
		elif (hdr['CTYPE1'] == 'RA---TAN' and hdr['CTYPE2'] == 'DEC--TAN'):
			col, row = w.wcs_world2pix(args.ra,args.dec,1)
		else:
			print("\nERROR: No valid coordinates in .fits header (check CTYPE(s) are set to 'RA---TAN' & 'DEC--TAN')\n  -- ABORTING --\n"); exit()

		#print("row: {0}\ncol: {1}".format(col,row))

		if(vrb>0):print("INFO: Converting RA/DEC to pixel coordinates.\n       RA = {0} -> {1}\n       Dec = {2} -> {3}".format(args.ra,col,args.dec,row))

	else: 
		print("\nERROR: No valid position (x/y | ra/dec) coordinates given.\n -- ABORTING --"); exit()

	if(vrb>0):print("\nINFO: row = {0}\nINFO: col = {1}".format(row,col))

	# Get appropriate frequency band
	bandRef = {"u":1,
			   "g":2,
			   "r":3,
			   "i":4,
			   "z":5}

	try:
		bandID = bandRef[args.band]
	except(KeyError):
		print("\nERROR: Invalid band \"{0}\" given.".format(args.band))


	cmd = "{0} \"{1}\" {2} {3} {4} \"{5}\"".format(PSF_PATH,inputFile,bandID,row,col,outputFile)
	if(vrb):print("\nINFO: ATLAS read_PSF command:\n  {0}\n".format(cmd))
	os.system(cmd)

	if (args.trim):
		if(vrb):print("INFO: Trim PSF = True")
		if(vrb):print("INFO: Reading in created PSF \"{0}\"".format(outputFile))
		
		# Open created PSF file
		PSFFile = pyfits.open(outputFile)
		img = np.array(PSFFile[0].data)
		hdr = PSFFile[0].header

		if(vrb):print("INFO: Finding PSF padding edges.")
		edges = find_edges(img)
		if(vrb):print("INFO: Found edges at ({0[0]},{0[1]})".format(edges))

		cropImg = img[edges[0][0]:edges[1][0]+1,edges[0][1]:edges[1][1]+1]

		if(vrb):print("INFO: Saving trimmed PSF as \"{0}\"".format(outputFile))
		newhdu = pyfits.PrimaryHDU(cropImg)
		hdulist = pyfits.HDUList([newhdu])
		hdulist.writeto(outputFile,overwrite=True)



if __name__ == "__main__": main()