"""
Script to create a Results-***.csv file from a set of optimisation runs

Author: Robin H. W. Cook
Date: 17/05/18
"""

import os
from argparse import ArgumentParser

import pandas as pd
from astropy.io import fits
import numpy as np

from IPython.display import display

def get_colnames(bands,comps,inits,sfxs,extras): # Determine the column names for the output table
	"""
	<param: bands (list)> - The list of bands
	<param: comps (list)> - The list of component numbers
	<param: inits (bool)> - Whether there should be initial guesses in the file
	<param: sfxs (list)> - The suffixes to be used for the bulge and disk component 
	<param: extras (list)> - The list of extra column names
	
	<return: colNames (list)> - The list of column names for the table
	"""

	paramNames = ["xcen","ycen","mag","re","nser","ang","axrat","box"]
		
	colNames = ["GASSID"]
	for band in bands:
		for comp in comps:
			if (inits == True): # First add the initial guess (if specified)
				for par in paramNames:
					if (comp==1):
						colNames.append("{p}_{b}0".format(p=par,b=band)) # single component
					if (comp==2):
						colNames.append("{p}{s}_{b}0".format(p=par,s=sfxs[0],b=band)) # bulge
						colNames.append("{p}{s}_{b}0".format(p=par,s=sfxs[1],b=band)) # disk

			for par in paramNames: # Add the optimised parameters
					if (comp==1):
						colNames.append("{p}_{b}".format(p=par,b=band)) # single component
					if (comp==2):
						colNames.append("{p}{s}_{b}".format(p=par,s=sfxs[0],b=band)) # bulge
						colNames.append("{p}{s}_{b}".format(p=par,s=sfxs[1],b=band)) # disk

			for ext in extras: # add extra column names
				colNames.append("{e}_{b}_{n}C".format(e=ext,b=band,n=comp))

	return(colNames)


def get_result(model,inits,comp,extras): # Extra the necessary values from an individual galaxy result file
	"""
	<param: model (pandas.DataFrame)> - The data frame containing a galaxy fitting result
	<param: inits (bool)> - Whether there should be initial guesses in the file
	<param: comp (list)> - The number of component fitted in this run
	<param: extras (list)> - The list of extra column names
	
	<return: gal (list)> - The list of column names for the table
	"""


	# [val for pair in zip(l1, l2) for val in pair]

	## Get the list of the extra column values that are present in the file, else set them to NaN
	extraVals = [model[[ext]].values.tolist()[0][0] if ext in list(model) else np.nan for ext in extras] # check which extra columns actually exist in this table

	params = ["x","y","mag","re","nser","ang","axrat","box"]

	## Get the model parameters
	paramsList = []
	if (comp == 1): # 1comp model
		if (inits == True): paramsList = paramsList + ["{0}_in".format(par) for par in params]
		paramsList = paramsList + ["{0}_out".format(par) for par in params]
		
	else: # 2comp model
		if (inits == True): paramsList = paramsList + [val for pair in zip(["{0}1_in".format(par) for par in params], ["{0}2_in".format(par) for par in params]) for val in pair]
		paramsList = paramsList + [val for pair in zip(["{0}1_out".format(par) for par in params], ["{0}2_out".format(par) for par in params]) for val in pair]

	result = model[paramsList].values.tolist()[0] + extraVals

	return(result)

def get_num_params(inits,comp,extras): # Get the number of parameters (optimised/[initial]) for the specified fit

	"""
	<param: inits (bool)> - Whether there should be initial guesses in the file
	<param: comp (list)> - The number of component fitted in this run
	<param: extras (list)> - The list of extra column names

	<return: numParams (int)> - The number of parameters for the specified model fit
	"""

	# 8 [parameters] x [num. components] x [initial values or not] + [num. extra columns]

	numParams = 8 * comp * (2 if inits else 1) + len(extras)
	return(numParams)


def main():

	GALS_DATA = "/home/robincook/Google Drive/PhD/GASS/Catalogues/xGASS_SDSSPhoto.csv"
	GALS_DIR = "/home/robincook/Documents/PhD/GASS/Galaxies"

	parser = ArgumentParser(description="Searches all directories for a particular optimisation run and concatenates results into a table.")

	inputGroup = parser.add_argument_group(title="Input Options.",description="Which galaxies to write to the results table.")
	selectGroup = parser.add_argument_group(title="Model Selection Options.",description="Selections for which model fits to be written to the table.")
	outputGroup = parser.add_argument_group(title="Output Options.",description="Options for the output of the results table.")
	miscGroup = parser.add_argument_group(title="Miscellaneous Options.",description="Miscellaneous controls for the program.")

	inputGroup.add_argument('-d','--galaxy-dir',
						action='store',dest='galsDir',type=str,default=GALS_DIR,
						help="The directory containing subdirectories for each galaxy (Default: \"{0}\").".format(GALS_DIR),metavar="PATH")
	inputGroup.add_argument('-f','--file',
						action='store',dest='file',type=str,default=None,
						help="A file containing a list of galaxy IDs. If none given, use all galaxies in --galaxy-dir.",metavar="FILE")
	inputGroup.add_argument('-g','--galaxies',
						action='store',dest='galaxies',nargs='+',default=None,
						help="A list of galaxies to be displayed",metavar="GAL")

	selectGroup.add_argument('-r','--run',
						type=str,action='store',dest='run',default=None,
						help="The optimisation run to be searched for and concatenated into a table.",metavar="RUN")
	selectGroup.add_argument('-b','--bands',
						action='store',dest='bands',type=str,nargs='*',default=['r'],choices=['u','g','r','i','z'],
						help="The SDSS filter (['u','g','r','i','z']).",metavar="BAND")
	selectGroup.add_argument('-c','--components',
						action='store',dest='comps',type=int,nargs='*',default=[1,2],choices=[1,2],
						help="The number of components in the fit.",metavar="COMPS")

	outputGroup.add_argument('-i','--inits',
						action='store_true',dest='inits',default=False,
						help="Whether to write initial guesses in the output file.")
	outputGroup.add_argument('-e','--extra',
						action='store',dest='extras',type=str,nargs='*',default=[],
						help="The names of the extra columns to be added to the output table.",metavar="EXTRA")
	outputGroup.add_argument('-s','--suffix',
						action='store',dest='sfxs',type=str,nargs=2,default=["1","2"],
						help="The suffix to be used for the bulge and disk component (Default: Bulge=1; Disk=2)",metavar="SUFFIX")
	outputGroup.add_argument('-t','--tables',
						action='store',dest='tables',type=str,nargs="+",default=[GALS_DATA],
						help="The additional data tables to crossmatch the results table with (Default: \"{0}\").".format(GALS_DATA),metavar="TABLE")
	outputGroup.add_argument('-o','--output',
						action='store',dest='output',type=str,default=None,
						help="The name of the output file.",metavar="PATH")

	miscGroup.add_argument('-v',
						action='count',dest='verbosity',default=0,
						help="The level of verbosity to print to stdout. if -v > 1; the HTML file will be opened after creation.")

	### PARSE ARGUMENTS ###
	args = parser.parse_args()

	## Validate verbosity
	global vrb;	vrb = args.verbosity if (args.verbosity >=0 and args.verbosity <=3) else 3

	## Validate galaxy directory
	if not (os.path.exists(args.galsDir)):
		print("\nERROR: Galaxy directory: \"{0}\" does not exist!\n  -- ABORTING --\n".format(args.galsDir)); exit()
	else:
		galsDir = args.galsDir if args.galsDir[-1] != '/' else args.galsDir[:-1]
		if (vrb>0): print("INFO: Galaxy working directory:\n  \"{0}\"".format(galsDir))

	## Get the list of all galaxies that exist in this directory
	allGals = [gal for gal in os.listdir(galsDir) if '.' not in gal]

	## Input choices:
	if (args.run):
		run = args.run
	else:
		print("\nERROR: No run name given!\n  -- ABORTING --\n"); exit()

	bands = args.bands
	comps = args.comps

	if (vrb>0): 
		print("INFO: Working with run \"{0}\".".format(run))
		print("INFO: SDSS band: {0}".format(bands))
		print("INFO: Number of components: {0}".format(",".join([str(c) for c in comps])))

	## Validate component suffixes
	if (args.sfxs[0] == args.sfxs[1]):
		print("\nERROR: Component suffixes cannot be the same!\n  -- ABORTING --\n"); exit()
	else:
		sfxs = args.sfxs
	if(vrb>0): print("INFO: Component suffixes:\n Bulge = \"{0}\"\n Disk = \"{1}\"".format(sfxs[0],sfxs[1]))

	## Output file:
	if args.output != None:
		outFilename = args.output
	else:
		bandStr = "".join(bands)
		compStr = "_{n}comp".format(n=comps[0]) if len(comps) == 1 else ""
		outFilename = '/home/robincook/Google Drive/PhD/Fitting/Results/Results-{run}_{band}{comp}.csv'.format(run=run,band=bandStr,comp=compStr)

	## Get list of galaxies
	galNames = []
	if (args.galaxies and args.file): print("\nWARNING: Both '--galaxies' and '--inputfile' were specified. Keeping --galaxies."); args.file = None
	if (args.galaxies == None and args.file == None): galNames = allGals # If no --inputfile OR --galaxies given, assume all galaxies in the --galaxy-dir

	if (args.galaxies): # If --galaxies specified
		if (vrb>0): print("INFO: Using the list of galaxies given.")
		for gal in args.galaxies:
			galNames.append(gal)

	if (args.file): # If --inputfile was given
		if (vrb>0): print("INFO: Getting list of galaxies from \"{0}\"".format(args.file))
		if not os.path.exists(args.file): print("\nERROR: file \"{file}\" does not exist!\n  -- ABORTING --\n".format(file=args.file)); exit()
		
		# Get galaxy names and check if they exist within the gicen galaxy directory
		galFile = open(args.file,'r')
		galNames = [gal.replace("\n","") for gal in galFile.readlines()] 
		galFile.close()

	## Get the list of galaxies that aren't in the galaxy directory
	noGal = list(set(allGals) - set(galNames))

	## Validate whether these galaxies have the appropriate directory/files
	noFitting = [] # No Fitting file
	noRun = [] # Does not contain the fitting Runile

	for gal in galNames:
		if ('Fitting' not in os.listdir('{dir}/{gal}'.format(dir=galsDir,gal=gal))): # Exclude galaxies without fits
			print("WARNING: Directory \".../{gal}/Fitting\" does not exist!".format(gal=gal))
			noFitting.append(gal)
			continue
		if (run not in os.listdir('{dir}/{gal}/Fitting'.format(dir=galsDir,gal=gal))): # Exclude those witout Run folder
			print("WARNING: Run \".../{gal}/Fitting/{run}\" does not exist!".format(gal=gal,run=run))
			noRun.append(gal)
			continue

	if (vrb>1 and len(noFitting)>0): print("\nWARNING: Galaxies without the \".../Fitting\" directory: {0}".format(",".join(noFitting)))
	if (vrb>1 and len(noRun)>0): print("\nWARNING: Galaxies without a \".../Fitting/{0}\" run directory: {1}\n".format(run,",".join(noRun)))

	## Remove galaxies with no fitting or run directory
	galList = [gal for gal in galNames if gal not in noGal+noFitting+noRun]
	if (vrb>1): print("INFO: Num galaxies: {0}".format(len(galNames)))

	## Check if at least one valid galax was given and if not, Abort.
	if (len(galList) == 0): print("ERROR: No valid galaxies given.\n  -- ABORTING --"); exit()

	## Get column names
	extras = args.extras
	colNames = get_colnames(bands,comps,args.inits,sfxs,extras)
	if (vrb>1): print("\nINFO: Column names:\n  ",", ".join(colNames))

	## Create table data frame
	df = pd.DataFrame(columns=colNames)

	## Loop over all galaxies and create table rows
	for gal in galList:
		if (vrb>2): print(gal+": ",end='')
		galRow = [gal] # initialise galaxy row with the galaxy ID
		for band in bands: # Loop over bands
			if (vrb>2): print(band,end='-')
			for comp in comps: # Loop over components
				if (vrb>2): print(comp,end=' ')
				resultFile = "{d}/{g}/Fitting/{r}/{g}-{r}_{b}_{n}comp_Output.csv".format(d=galsDir,g=gal,r=run,b=band,n=comp)
				if (os.path.exists(resultFile)):
					# Add this particular model fit (band & nComps) for this galaxy to the current galaxy row
					resultData = pd.read_csv(resultFile) # Read the second line
					galRow = galRow + get_result(resultData,args.inits,comp,extras)
				else:
					# Add a set of nans for the appropriate length of the model
					if(vrb>1): print("\n WARNING: \"{0}\" does not exist. Appending NaNs.".format(resultFile))
					galRow = galRow + [np.nan]*get_num_params(args.inits,comp,extras) # 8 parameters in each component + the extra columns
		
		## Append this galaxy to the data frame
		if(vrb>2):print(flush=True)
		df.loc[len(df)] = galRow


	## Cross match result table with additional tables
	if (len(args.tables)>0 and vrb): print("\nINFO: Cross-matching additional tables.")
	for tab in args.tables:
		if (os.path.exists(tab)):
			## Open the table
			if (".fits" in tab): # If a .fits file
				with fits.open(tab) as tabData:
					dfTable = pd.DataFrame(tabData[0].data)
			elif (".csv" in tab): # If a .csv file
				dfTable = pd.read_csv(tab)
			else: # All other files
				dfTable = pd.read_table(tab)

			## Merge results table with given data table matching on the GASSIDs of the results table.
			if(vrb): print("INFO: Cross-matching with \"{0}\".".format(tab))
			df = pd.merge(df,dfTable,left_on="GASSID",right_on="GASSID" if "GASSID" in list(dfTable) else list(dfTable)[0],how='left')

		else:
			print("\nWARNING: Table \"{0}\" does not exist!\n  -- Ignoring this table --\n".format(tab))


	## Save the data frame to file
	if(vrb): print("INFO: Saving results table to \"{0}\".".format(args.output))
	df.to_csv(args.output,na_rep='nan',index=False)


if __name__ == "__main__": main()