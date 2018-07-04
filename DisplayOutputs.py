"""
A script to write a nicely formatted HTML page of a series of outputs from running ProFit on galaxies.

Author: Robin H. Cook
Date: (<) 31/05/17
"""

import time
import os
import math
from astropy.table import Table
from argparse import ArgumentParser

import subprocess

HOME = str(subprocess.check_output("pwd")).split('/')[1]
GALS_DATA = "/{home}/robincook/Google Drive/PhD/GASS/Catalogues/xgass_xcoldgass_master_09052017.fits".format(home=HOME)

# Define colours for frequency bands
global colorRef
colorRef = {'u':'cyan','g':'ForestGreen','r':'red','i':'magenta','z':'purple','bulge':'Firebrick','disk':'MidnightBlue'}

def format_filename(name,dir=None,id=None,run=None,band=None,nComps=None,ext="png"):
	return("{d}/{r}/{id}-{r}_{b}_{n}comp_{name}.{ext}".format(d=dir,id=id,r=run,b=band,n=nComps,name=name,ext=ext))

class Galaxy:
	""" An object defining a Galaxy and its outputs """

	def __init__(self,ID,SDSSname,ra,dec,bands,comps,extras,runs=['']):
		self.id = ID
		self.name = SDSSname
		self.ra = ra
		self.dec = dec
		self.bands = bands
		self.comps = comps
		self.runs = self.get_runs(runs)
		self.nBands = len(self.bands)
		self.nComps = len(self.comps)
		self.nRuns = len(self.runs)
		self.nRows = self.nBands*self.nComps*self.nRuns
		self.dir = "{0}/{1}/Fitting".format(galsDir,self.id)
		self.extras = extras

	def get_runs(self,runs):
		self.runs = []
		for run in runs:
			if (run != ''): self.runs.append(run)

		return self.runs

	def write_row(self):

		showSDSS = True # Whether to show an SDSS colour image of the galaxy
		imSize = 225 # the size of the SDSS cutout image
		scaleIm = True # Whether to scale the SDSS cutout image by galaxy size

		line = "" # Initialise the the line of .html for this galaxy
		
		# Make left-most galaxy information panel
		line = line + ('\n\n<!-- {G.id} -->\n' 
						   '<tr align="center">\n'
						   '\t<td valign="top" rowspan="{G.nRows}">\n\n'
						   '\t\t<!-- Galaxy Information -->\n'
						   '\t\t<table bordercolor="black" border="2" align="center">\n'
						   '\t\t\t<tr align="center">\n'
						   '\t\t\t\t<td><b>Galaxy ID:</b> {G.id}</td>\n'
						   '\t\t\t</tr>\n'
						   '\t\t\t<tr align="center">\n'
						   '\t\t\t\t<td> <b>Name:</b> {G.name}</td>\n'
						   '\t\t\t</tr>\n'
						   '\t\t\t<tr align="center">\n'
						   '\t\t\t\t<td> <b>RA:</b> {G.ra:.5f} <b>Dec:</b> {G.dec:.5f}</td>\n'
						   '\t\t\t</tr>\n'.format(G=self))

		if (showSDSS): # Whether to show an SDSS colour image
			if (scaleIm == True): # Attempt to make image scale based on a 1-component model of the galaxy
				imScale = 0.792 # Assume default scaling
				for run in self.runs: # Try to find a run with a single-component fit
					if os.path.exists(format_filename("Output",ext="csv",dir=self.dir,id=self.id,run=run,band=self.bands[0],nComps=1)):
						tempModel = Table.read(format_filename("Output",ext="csv",dir=self.dir,id=self.id,run=run,band=self.bands[0],nComps=1))
						#tempModel = Table.read("{G.dir}/{run}/{G.id}-{run}_r_1comp_Output.csv".format(G=self,run=run))
						imScale = math.sqrt(tempModel["re_out"][0])/12.0
						break # found one! exit loop.
			
			line = line + ('\t\t\t<tr><!-- SDSS image -->\n'
						   '\t\t\t\t<td><div class="rotate">\n'
						   '\t\t\t\t\t<iframe frameborder="2" width="{size}" height="{size}" src="http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra={G.ra}&dec={G.dec}&width={size}&height={size}&scale={scale:.3f}" name="{G.id}image" id="{G.id}image"> </iframe>\n'
						   '\t\t\t\t</div></td>\n'
						   '\t\t\t</tr>\n'
						   '\t\t</table>\n'.format(G=self,size=imSize,scale=imScale))

		# Add extra lines from the input file if specified
		if (len(self.extras) > 0):
			for ext in self.extras:
				line = line + ('\t\t<br>\n'
					           '\t\t<b>{0}</b>'.format(ext))


		# Add a judgement table if specified
		if (judge):
			line = line + ('\t\t<br>\n'
					   '\t\t<form action="process.php" method="POST" align="left">\n'
					   '<br><br>\n'
					   '\t\t\t<b>No. comps: </b>\n'
					   '\t\t\t<input type="radio" name="{G.id}_nComps" value="1">1\n'
					   '\t\t\t<input type="radio" name="{G.id}_nComps" value="2">2\n'
					   '\t\t\t<input type="radio" name="{G.id}_nComps" value="3">3\n'
					   '\t\t\t<input type="radio" name="{G.id}_nComps" value="0">?\n'
					   '\t\t\t<br>\n'
					   '\t\t\t<b>Quality: </b>\n'
					   '\t\t\t<input type="radio" name="{G.id}_qlty" value="good"><u>G</u>ood\n'
					   '\t\t\t<input type="radio" name="{G.id}_qlty" value="okay"><u>O</u>kay\n'
					   '\t\t\t<input type="radio" name="{G.id}_qlty" value="bad"><u>B</u>ad\n'
					   '\t\t\t<br><b> Problems: </b><br>\n'
					   '\t\t\t&nbsp&nbsp<input type="checkbox" name="{G.id}_comp"> # of <u>C</u>omponents <br>\n'
					   '\t\t\t&nbsp&nbsp<input type="checkbox" name="{G.id}_seg"> <u>S</u>egmentation <br>\n'
					   '\t\t\t&nbsp&nbsp<input type="checkbox" name="{G.id}_sky"> <u>B</u>ackground subtraction <br>\n'
					   '\t\t\t&nbsp&nbsp<input type="checkbox" name="{G.id}_sky"> <u>D</u>ifficult galaxy <br>\n'
					   '\t\t\t&nbsp&nbsp<input type="checkbox" name="{G.id}_sky"> <u>I</u>mage problems <br>\n'
					   '\t\t</form>\n'.format(G=self))

		line = line + '\t</td>\n' # Close first table entry

		for run in self.runs: # For each run specified

			runStr = "({run})".format(run=run)
			runDir = "/{run}".format(run=run)

			for band in self.bands:
				# <!-- **************** -->
				for comp in self.comps:
					# Attempt to find initial/optimised model file
					#if os.path.exists("{G.dir}{rdir}/{G.id}-{run}_{band}_{comp}comp_Output.csv".format(G=self,comp=comp,band=band,run=run,rdir=runDir)):
					if os.path.exists(format_filename("Output",ext="csv",dir=self.dir,id=self.id,run=run,band=band,nComps=comp)):
						modelFound = True
						model = Table.read(format_filename("Output",ext="csv",dir=self.dir,id=self.id,run=run,band=band,nComps=comp))
					else:
						modelFound = False

					line =line+('\n\t<td>\n'
								'\t\t<center><font size="+2" color="{col}"><b>{run}<br>{n}comp<br>{b} band</b></font></center>\n'
								'\t</td>\n'.format(col=colorRef[band],n=comp,run=run,b=band))


					for disp in displays:
						if (disp in ["ModelInital","ModelOptimised"] and modelFound):
							# begin table
							line = line + '\n\t<td align="center" valign="top"> <!-- {G.id}: {col} ({comp}-component) -->\n'.format(G=self,rdir=runDir,comp=comp,band=band,run=run,col=colNames[disp],disp=disp)
							line = line + ('\t\t<table bgcolor=\"lightgray\" bordercolor=\"black\" align=\"center\" border=\"2\" cellspacing=\"1\" cellpadding=\"5\">\n'
										   '\t\t\t<tr>\n'
										   '\t\t\t\t<th colspan="12"><font color=\"{col}\" size="+1"> <b>{run}</b> </font></th>\n' # Add the run name
										   '\t\t\t</tr>\n'
										   '\t\t\t<tr>\n'
										   '\t\t\t\t<th> comp. </th>\n' # which comp
										   '\t\t\t\t<th> x<sub>cen</sub> </th>\n' # x pos
										   '\t\t\t\t<th> y<sub>cen</sub> </th>\n' # y pos
										   '\t\t\t\t<th> mag </th>\n' # magnitude
										   '\t\t\t\t<th> r<sub>e</sub> </th>\n' # effective radius
										   #'\t\t\t\t<th> &lt;&#956;<sub>e</sub>&gt;' # surface brightness within Re
										   '\t\t\t\t<th> n </th>\n' # Sersic Index
										   '\t\t\t\t<th> &#952; </th>\n' # Position angle (\theta)
										   '\t\t\t\t<th> b/a </th>\n' # Axial ratio
										   '\t\t\t\t<th> box </th>\n' # Boxiness
										   '\t\t\t\t<th> &#967;<sup>2</sup> </th>\n' # Chi squared value
										   '\t\t\t\t<th align=\"center\"> time (hr) </th>\n' # Chi squared value
										   '\t\t\t\t<th align=\"center\"> stationary </th>\n' # Stationarity
										   '\t\t\t</tr>\n'.format(run=run,col=colorRef[band]))	
							if (comp == 1):
								# Add first (and only) component
								line = line + '\t\t\t<tr>\n'
								line = line + "\t\t\t\t<td> <i><b>total</b></i> </td>\n"
								for col in model.colnames[9:-3]:
									line = line + "\t\t\t\t<td align=\"center\"> {0:.2f} </td>\n".format(model[col][0])

								line = line + "\t\t\t\t<td align=\"center\"> {0:.2f} </td>\n".format(model["chisq"][0])
								line = line + "\t\t\t\t<td align=\"center\"> {0:.2f} </td>\n".format(model["elapsed_time"][0]/3600)
								line = line + "\t\t\t\t<td align=\"center\"> {0} </td>\n".format("True" if model["stationarity"][0] == "1" else "False")
								line = line + '\t\t\t</tr>\n'

							elif (comp == 2):
								# Add first component
								line = line + '\t\t\t<tr>\n'
								line = line + "\t\t\t\t<td> <i><font color=\"{col}\"> <b>inner</b> </font></i> </td>\n".format(col=colorRef['bulge'])
								for col in model.colnames[17:-3:2]: # From index 17 (x1_out) to the end by steps of 2
									line = line + "\t\t\t\t<td align=\"center\"><font color=\"{col}\"> {param:.2f} </font></td>\n".format(param=model[col][0],col=colorRef['bulge'])

								line = line + "\t\t\t\t<td rowspan=\"2\" align=\"center\"> {param:.2f} </td>\n".format(param=model["chisq"][0])
								line = line + "\t\t\t\t<td rowspan=\"2\" align=\"center\"> {param:.2f} </td>\n".format(param=model["elapsed_time"][0]/3600)
								line = line + "\t\t\t\t<td rowspan=\"2\" align=\"center\"> {param} </td>\n".format(param="True" if model["stationarity"][0] == 1 else "False")
								line = line + '\t\t\t</tr>\n'

								# Add second component
								line = line + '\t\t\t<tr>\n'
								line = line + "\t\t\t\t<td> <i><font color=\"{col}\"> <b>outer</b> </font></i> </td>\n".format(col=colorRef['disk'])
								for col in model.colnames[18:-3:2]: # From index 18 (x2_out) to the end by steps of 2
									line = line + "\t\t\t\t<td align=\"center\"><font color=\"{col}\"> {param:.2f} </font></td>\n".format(param=model[col][0],col=colorRef['disk'])

								line = line + '\t\t\t</tr>\n'

								# Add ratios row for Bulge-to-Total and Re1/Re2
								line = line + ('\t\t\t<tr>\n'
											  '\t\t\t\t<td> <i><b>diff</b></i> </td>\n'
											  '\t\t\t\t<td align=\"center\"> &#916;: {dx:0.2f} </td>\n'
											  '\t\t\t\t<td align=\"center\"> &#916;: {dy:0.2f} </td>\n'
											  '\t\t\t\t<td align=\"center\"> &#8453;: {b2t:0.2f} </td>\n'
											  '\t\t\t\t<td align=\"center\"> &#189;: {reRat:0.2f} </td>\n'
											  '\t\t\t\t<td align=\"center\" colspan="7" style="border-left:none;border-bottom:none;border-right:none"> </td>\n'
											  '\t\t\t</tr>\n').format(dx = model["x1_out"][0] - model["x2_out"][0],
											  						 dy = model["y1_out"][0] - model["y2_out"][0],
											  						 b2t = calc_b2t(model["mag1_out"][0],model["mag2_out"][0]),
											  						 reRat = model["re1_out"][0] / model["re2_out"][0])

								
								
							line = line + ('\t\t</table>\n'
										   '\t</td>\n')
										   
							# name,dir=None,id=None,run=None,band=None,nComps=None
						else:
							line = line + ('\n\t<td align="center"> <!-- {G.id}: {col} ({comp}-component) -->\n'
										   '\t\t<a target="_blank" href="{file}">\n'
										   '\t\t\t<img src="{file}" height="300">\n'
										   '\t\t</a>\n'
										   '\t</td>\n'
										   ''.format(G=self,col=colNames[disp],comp=comp,file=format_filename(disp,self.dir,id=self.id,run=run,band=band,nComps=comp)))



					line = line + '\n</tr>'
					
			
		if(vrb>2):print(line)
		return line

	def write_model(self,model,nComps,state):
		""" Writes the model into hmtl table rows
		<param: model [Table]> - The model data table
		<param: nComps [int]> - the number of components
		<param: state [int]> - whether the model is the initial (0) or optimised (1)

		<return: line> - the HTML formatted line for the component.
		"""
		line = "    <tr>\n"

		nCols = len(model.colnames[2:])
		for col in model.colnames[2:]: # Write column names
			line = line + "      <th> {0} </th>\n".format(col)

		line = line + "    </tr>\n"
		
		for comp in range(nComps):
			line = line + "    <tr>\n"
			for col in model.colnames[2:]:
				line = line + "      <th> {0:.2f} </th>\n".format(model[col][comp+nComps*state])			

			line = line + "    </tr>\n"

		return line

def calc_b2t(mag1,mag2):
	"""
	Caclulates the Bulge-to-Total ratio of two components.

	<param: mag1 [float]> - The magnitude of the first component.
	<param: mag2 [float]> - The magnitude of the second component.

	<return: b2t [float]> - The bulge-to-total ratios of the fluxes.
	"""

	# convert to fluxes using standard equation.
	flux1 = 10**(-0.4*mag1)
	flux2 = 10**(-0.4*mag2)

	# calculate B/T ratio.
	b2t = flux1/(flux1+flux2)

	return(b2t)


def main():

	showInput = True
	showInit = True
	showOptim = False
	showEllipse = True

	parser = ArgumentParser()
	parser.add_argument('-g','--galaxies',
						action='store',dest='galaxies',nargs='+',default=None,
						help="A list of galaxies to be displayed (If none given, print list and exit).",metavar="GALAXY")
	parser.add_argument('-i','--inputfile',
						action='store',dest='file',type=str,default=None,
						help="A file containing a list of galaxy IDs.",metavar="FILE")
	parser.add_argument('-b','--bands',
						action='store',dest='bands',nargs='+',default=['r'],
						help="A list of bands to be displayed for each galaxy [u|g|r|i|z].",metavar="BAND")
	parser.add_argument('-c','--components',
						action='store',dest='components',nargs='+',type=int,default=[1],
						help="A list of the numbers of components for each fit (typically 1 and/or 2 components)",metavar="COMPONENT")
	parser.add_argument('-r','--runs',
						action='store',dest='runs',nargs='+',default=[''],
						help="A list of different runs of ",metavar="RUN")
	parser.add_argument('-d','--displays',
						action='store',dest='displays',type=str,nargs='*',default=["Inputs","LikelihoodOptimised","EllipseOptimised","ModelOptimised"],#choices=["Inputs","EllipseInitial","EllipseOptimised","LikelihoodInitial","LikelihoodOptimised","ModelInitial","ModelOptimised","ChiSquaredInitial","ChiSquaredOptimised","SkyStats","CornerPlot"],
						help="The types of plots to be displayed given as a list of filename extensions (choices = Inputs EllipseInitial EllipseOptimised LikelihoodInitial LikelihoodOptimised ModelInitial ModelOptimised ChiSquaredInitial ChiSquaredOptimised SkyStats CornerPlot).",metavar="EXTENSION")
	parser.add_argument('-e','--extra',
						action='store',dest='extra',type=str,nargs='*',default=None,
						help="Specifies any extra columns in the input table (--inputfile/-i) to be written in the first panel.",metavar="EXTRA")
	parser.add_argument('-j','--judge',
						action='store',dest='judge',type=str,default=None,
						help="A judgment form will be created for each galaxy and the results sent to this file.",metavar="FILENAME")
	parser.add_argument('-o','--outputfile',
						action='store',dest='output',type=str,default="/{home}/robincook/Google Drive/PhD/GASS/Optimisations/GASSFittingSummary.html".format(home=HOME),
						help="The filename of the output file",metavar="FILENAME")
	parser.add_argument('-t','--title',
						action='store',dest='title',type=str,default="Galaxy Model Fitting Summary",
						help="The title to be displayed on the HTML page.",metavar="TITLE")
	parser.add_argument('-s','--subtitle',
						action='store',dest='subtitle',type=str,default="",
						help="The subtitle or explanation paragraph to be displayed on the HTML page.",metavar="SUBTITLE")
	parser.add_argument('-n','--number',
						action='store',dest='num',type=int,default=1,
						help="The number of lists to split up the larger list between (Default=1).",metavar="NUMBER")
	parser.add_argument('-D','--galaxy-dir',
						action='store',dest='galsDir',type=str,default="/home/robincook/Documents/PhD/GASS/Galaxies",
						help="A file containing a list of galaxy IDs.",metavar="FILE")
	parser.add_argument('-l','--list',
						action='store_true',dest='list',default=False,
						help="Print available galaxies and exit.")
	parser.add_argument('-v',
						action='count',dest='verbosity',default=0,
						help="The level of verbosity to print to stdout. if -v > 1; the HTML file will be opened after creation.")

	
	### PARSE ARGUMENTS ###
	args = parser.parse_args()

	## Validate verbosity
	global vrb;	vrb = args.verbosity if (args.verbosity >=0 and args.verbosity <=2) else 2

	## Get the galaxy directory
	global galsDir; galsDir = args.galsDir

	## Get all galaxy IDs
	galData = Table.read(GALS_DATA)
	allGals = ["{0}".format(ID).replace(" ","") for ID in galData[galData.colnames[0]]] # Assumes first column (index=0) contains the galaxy IDs

	## Print list of galaxies and exit.
	if (args.list): print("\nINFO: Galaxy IDs for available galaxies:\n{0}\n".format(' '.join(allGals))); exit()

	## Get galaxy list
	if (args.galaxies != None and args.extra != None): print("\nERROR: Cannot specify both --galaxies and --extra.\n  -- ABORTING --\n"); exit()
	if (args.galaxies == None and args.file == None): print("\nERROR: No '--galaxies' or '--file' specified by user!\n  -- ABORTING --\n"); exit()
	if (args.galaxies and args.file): print("\nWARNING: Both '--galaxies' and '--inputfile' were specified. Keeping --inputfile."); args.galaxies = None

	galNames = []
	if (args.galaxies): # If galaxy specified
		for gal in args.galaxies:
			if (gal not in allGals): print("WARNING: Galaxy \"{0}\" does not exist.".format(gal))
			else: galNames.append(gal)

		extras = [] # No file was given, so there cannot be extra columns

	if (args.file):
		if not os.path.exists(args.file): print("\nERROR: file \"{file}\" does not exist!\n  -- ABORTING --\n".format(file=args.file)); exit()
		galInfo = Table.read(args.file,format='ascii.csv')
		galNames = [gal for gal in galInfo[galInfo.colnames[0]]]
		if (args.extra):
			extras = []
			for extra in args.extra: # Validate extra columns
				if (extra not in galInfo.colnames):
					print("\nWARNING: Extra column \"{0}\" does not exist in input table.\n  -- Removing from extras list --\n".format(extra))
				else:
					extras.append(extra)
		else:
			extras = []

	## Check if no valid galaxies given: exit
	if (len(galNames) == 0): print("ERROR: No valid galaxies given.\n  -- ABORTING --"); exit()

	## Get fit parameters
	compList = args.components
	bandList = args.bands
	runList = args.runs

	## Get the file to output judgements of the fits.
	global judge; judge = True if args.judge != None else False
	global judgeFile; judgeFile = args.judge

	### DEFINE GALAXY OBJECT LIST ###
	galList = []
	for ii, gal in enumerate(galNames):
		try:
			index = allGals.index(gal)
			galList.append(Galaxy(gal,galData['SDSS'][index],float(galData['RA'][index]),float(galData['DEC'][index]),bandList,compList,[galInfo[ext][ii] for ext in extras],runList))
		except (ValueError):
			print("\nWARNING: {gal} not found!\n".format(gal=gal))


	if(vrb>1):
		print("\nINFO: Galaxy list ({0})".format(len(galList)))
		for gal in galList: print("  - {0}".format(gal.id))
		print("\nINFO: Frequency bands")
		for band in bandList: print("  - {0}".format(band))
		print("\nINFO: Components")
		for comp in compList: print("  - {0}".format(comp))
		print("\nINFO: Optimisation runs")
		for run in runList: print("  - {0}".format(run))


	# Set predefined column names:
	global displays; displays = args.displays
	global colNames; colNames = {}
	for disp in displays:
		if   ("Inputs" in disp): colNames[disp] = "Inputs"
		elif ("EllipseInitial" in disp): colNames[disp] = "Initial Ellipse Plot"
		elif ("EllipseOptimised" in disp): colNames[disp] = "Optimised Ellipse Plot"
		elif ("LikelihoodInitial" in disp): colNames[disp] = "Initial Likelihood Plot"
		elif ("LikelihoodOptimised" in disp): colNames[disp] = "Optimised Likelihood Plot"
		elif ("ModelInitial" in disp): colNames[disp] = "Initial Model"
		elif ("ModelOptimised" in disp): colNames[disp] = "Optimised Model"
		elif ("SkyStats" in disp): colNames[disp] = "Sky Statistics"
		elif ("CornerPlot" in disp): colNames[disp] = "MCMC Posterior Distributions"
		elif ("MCMCPosteriors" in disp): colNames[disp] = "MCMC Posterior Distributions"
		elif ("SegMap" in disp): colNames[disp] = "Segmentation Map"
		elif ("Sigma" in disp): colNames[disp] = "Sigma Map"
		elif ("PSF" in disp): colNames[disp] = "PSF"
		else: colNames[disp] = disp


	# Divide main list into sublists
	nGals = len(galList)
	sizeArr = [int(nGals/args.num)+1 if kk < nGals%args.num else int(nGals/args.num) for kk in range(args.num)]

	subLists = []; index = 0
	for size in sizeArr:
		subLists.append(galList[index:index+size])
		index = index+size

	# base filename
	filename = args.output if ('.html' in args.output) else "{0}.html".format(args.output)
	for ii, subList in enumerate(subLists):

		### OPEN OUTPUT FILE ###
		if (args.num == 1):
			file = open(filename,'w')
			if(vrb>0):print("\nINFO: Output filename = \"{file}\"".format(file=filename))
		else:
			file = open(filename.replace('.html','_{0}.html'.format(ii+1)),'w')
			if(vrb>0):print("\nINFO: Output filename = \"{file}\"".format(file=filename.replace('.html','_{0}.html'.format(ii+1))))

		### WRITE HTML FILES ###
		## Write preamble.
		file.write("<!DOCTYPE html>\n"
				   "<html>\n"
				   "<head>\n"
				   "\t<title> {title} </title>\n".format(title=args.title.replace(" ","_")))

		file.write("\t\n <!-- Page breaking for table rows -->\n"
				   "\t<style type=\"text/css\">\n"
				   "\t\ttable { page-break-inside:auto;\n"
				   "\t\ttable-layout: fixed;}\n\n"
				   "\t\ttr    { page-break-inside:avoid; page-break-after:auto }\n"
				   "\t\tthead { display:table-header-group }\n"
				   "\t\ttfoot { display:table-footer-group }\n\n"
			   	   "\t\t.formcol { width: 1%; }\n\n"
				   "\t\tdiv.vertical\n"
				   "\t\t{\n"
				   "\t\t\ttransform: rotate(-90deg);\n"
				   "\t\t\t-webkit-transform: rotate(-90deg); /* Safari/Chrome */\n"
				   "\t\t\t-moz-transform: rotate(-90deg); /* Firefox */\n"
				   "\t\t\t-o-transform: rotate(-90deg); /* Opera */\n"
				   "\t\t\t-ms-transform: rotate(-90deg); /* IE 9 */\n"
				   "\t\t}\n\n"
				   "\t\tdiv.rotate\n"
				   "\t\t{\n"
				   "\t\t\ttransform: rotate(+90deg);\n"
				   "\t\t\t-webkit-transform: rotate(+90deg); /* Safari/Chrome */\n"
				   "\t\t\t-moz-transform: rotate(+90deg); /* Firefox */\n"
				   "\t\t\t-o-transform: rotate(+90deg); /* Opera */\n"
				   "\t\t\t-ms-transform: rotate(+90deg); /* IE 9 */\n"
				   "\t\t}\n"
				   "\t</style>\n"
				   "\t<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"/home/robincook/Google Drive/Media/Pictures/Icons/favicon.ico\"/>\n")

		file.write("</head>\n"
				   "<body>\n"
				   "\t<h1> {title} </h1>\n"
				   "\t<h2>\n"
				   "\t\tAuthor: Robin H. W. Cook<br>\n"
				   "\t\tDate:&nbsp;{date}&nbsp;-&nbsp;{time}\n"
				   "\t</h2>\n".format(title=args.title,date=time.strftime("%d/%m/%Y"),time=time.strftime("%I:%M:%S")))

		file.write("\t<br>\n"
				   "\t<h3>\n"
				   "\t\t<font color=\"navy\">{subtitle}</font>\n"
				   "\t</h3>\n"
				   "\t<br>\n".format(subtitle=args.subtitle.replace("\\n","<br>\n\t\t")))


		## Begin TABLE
		file.write("\n\n<table bgcolor=\"white\" bordercolor=\"black\" align=\"center\" border=\"5\" cellspacing=\"1\" cellpadding=\"10\">\n")
		file.write("<tr>\n")
		file.write("\t<th><font size=\"+2.5\">Object</font></th>\n")
		file.write("\t<th class=\"formcol\" valign=\"center\" align=\"center\"><font size=\"+2.0\">Run</font></th>")
		for disp in displays:
			file.write("\t<th><font size=\"+2.5\">{name}</font></th>\n".format(name=colNames[disp]))
		
		file.write('</tr>\n\n')


		## Write each Galaxy to table row
		for gal in subList:
			file.write(gal.write_row())

		## End Table
		file.write('\n\n</table>')


		## Write closing statements
		file.write("\n</body>\n"
				   "\n</html>\n")


		file.close()

		#Open with Google-Chrome directly afterwards!
		if (vrb > 1 and ii == 0): os.system("google-chrome \"{0}\"".format(filename))


if __name__ == "__main__": main()