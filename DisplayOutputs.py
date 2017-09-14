"""
A script to write a nicely formatted HTML page of a series of outputs from running ProFit on GASS galaxies

Author: Robin H. Cook
Date: (<) 31/05/17
"""

import time
import os
from astropy.table import Table
from argparse import ArgumentParser

import subprocess

HOME = str(subprocess.check_output("pwd")).split('/')[1]
GALS_DIR = "/{home}/robincook/Documents/PhD/GASS/Galaxies".format(home=HOME)
GALS_DATA = "/{home}/robincook/Google Drive/PhD/GASS/Catalogues/xgass_xcoldgass_master_09052017.fits".format(home=HOME)

# Define colours for frequency bands
global colorRef
colorRef = {'u':'magenta','g':'ForestGreen','r':'red','i':'cyan','z':'teal'}

class Galaxy:
	""" An object defining a Galaxy and its outputs """

	def __init__(self,GASSID,SDSSname,ra,dec,bands,comps,forms=['']):
		self.id = GASSID
		self.name = SDSSname
		self.ra = ra
		self.dec = dec
		self.bands = bands
		self.comps = comps
		self.forms = self.get_forms(forms)
		self.nbands = len(self.bands)
		self.ncomps = len(self.comps)
		self.nforms = len(self.forms)
		self.nrows = self.nbands*self.ncomps
		self.dir = "{0}/GASS{1}/Fitting".format(GALS_DIR,self.id)

	def get_forms(self,forms):
		self.forms = []
		for form in forms:
			if (form != ''): self.forms.append('-{0}'.format(form))

		return self.forms

	def write_row(self):

		line = ""
		for form in self.forms:
			formStr = "({form})".format(form=form[1:])
			formDir = "/{form}".format(form=form[1:])
			line = line + ('\n\n<!-- GASS{G.id} -->\n'
						   '<tr align="center">\n'
						   '\t<td valign="top" rowspan="{G.nrows}">\n\n'
						   '\t\t<!-- Galaxy Information -->\n'
						   '\t\t<font size="+3"><b>{form}</b></font>\n'
						   '\t\t<br><br>\n'
						   '\t\t<table bordercolor="black" border="2" align="center">\n'
						   '\t\t\t<tr align="center">\n'
						   '\t\t\t\t<td><b>GASS ID:</b> GASS{G.id}</td>\n'
						   '\t\t\t\t</tr>\n'
						   '\t\t\t<tr align="center">\n'
						   '\t\t\t\t<td> <b>SDSS:</b> {G.name}</td>\n'
						   '\t\t\t</tr>\n'
						   '\t\t\t<tr>\n'
						   '\t\t\t\t<td>\n'
						   '\t\t\t\t\t<iframe frameborder="2" width="225" height="225" src="http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra={G.ra}&dec={G.dec}&width=225&height=225" name="GASS{G.id}image" id="GASS{G.id}image"> </iframe>\n'
						   '\t\t\t\t</td>\n'
						   '\t\t\t</tr>\n'
						   '\t\t</table>\n'.format(G=self,form=formStr))

			if (judge): 
				line = line + ('\t\t<br>\n'
						   '\t\t<form action="process.php" method="POST" align="left">\n'
						   '\t\t\t<b>No. comps: </b>\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_nComps" value="1">1\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_nComps" value="2">2\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_nComps" value="3">3\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_nComps" value="0">?\n'
						   '\t\t\t<br>\n'
						   '\t\t\t<b>Quality: </b>\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_qlty" value="good">good\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_qlty" value="okay">okay\n'
						   '\t\t\t<input type="radio" name="GASS{G.id}_qlty" value="bad">bad\n'
						   '\t\t\t<br><b> Problems: </b><br>\n'
						   '\t\t\t&nbsp&nbsp<input type="checkbox" name="GASS{G.id}_seg"> segmentation <br>\n'
						   '\t\t\t&nbsp&nbsp<input type="checkbox" name="GASS{G.id}_comp"> no. components <br>\n'
						   '\t\t\t&nbsp&nbsp<input type="checkbox" name="GASS{G.id}_sky"> sky subtraction <br>\n'
						   '\t\t</form>\n'.format(G=self))

			line = line + '\t</td>\n'

			for band in self.bands:
				line = line + ('\n\t<td rowspan={G.ncomps} align="center"> <!-- GASS{G.id}: ProFit Inputs -->\n'
							   '\t\t<font size="+3" color="{col}" align="center"><b>({band}-band)</b></font><br><br>\n'
							   '\t\t<a target="_blank" href="{G.dir}{fdir}/GASS{G.id}_{band}_Inputs.png">\n'
							   '\t\t\t<img src="{G.dir}{fdir}/GASS{G.id}_{band}_Inputs.png" height="300">\n'
							   '\t\t</a>\n'
							   '\t</td>\n'
							   ''.format(G=self,fdir=formDir,band=band,col=colorRef[band]))

				# <!-- **************** -->
				
				for comp in self.comps:
					if os.path.exists("{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_Output.csv".format(G=self,comp=comp,band=band,form=form,fdir=formDir)):
						found = True
						model = Table.read("{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_Output.csv".format(G=self,comp=comp,band=band,form=form,fdir=formDir),format='ascii.csv')
					else:
						found = False


					for disp in displays:
						#print(disp)
						#line = line + ('\n\t<!-- GASS{G.id}: Initial Model ({comp}-component) -->\n'
						#			   '\t<td align="center">\n'
						#			   '\t\t<font size="+2" color="{col}"><b>{comp}-Component (initial) </b></font><br><br>\n'
						#			   '\t\t<table bordercolor="black" align="center" border="2" cellspacing="1" cellpadding="2">\n'
						#			   ''.format(G=self,comp=comp,form=form,col=colorRef[band]))

						#line = line + self.write_model(model,comp,state=0)

						#line = line + ('  </table>\n'
						#			'  <br><br>\n\n'
						#			'  <a target="_blank" href="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_LikelihoodInitial.png">\n'
						#			'    <img src="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_LikelihoodInitial.png" width="500">\n'
						#			'  </a>\n'
						#			'</td>'
						#			''.format(G=self,fdir=formDir,comp=comp,band=band,form=form,col=colorRef[band]))

						### OPTIMISED MODEL ###
						#line = line + ('\n  <!-- GASS{G.id}: Output Model ({comp}-component) -->\n'
						#			'<td align="center">\n'
						#			'  <font size="+2" color="{col}"><b>{comp}-Component (optimised) </b></font><br><br>\n'
						#			'  <table bordercolor="black" align="center" border="2" cellspacing="1" cellpadding="2">\n'
						#			''.format(G=self,comp=comp,form=form,col=colorRef[band]))

						#line = line + self.write_model(model,comp,state=1)

						#line = line + ('  </table>\n'
						#			'  <br><br>\n\n'
						#			'  <a target="_blank" href="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_LikelihoodOptimised.png">\n'
						#			'    <img src="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_LikelihoodOptimised.png" width="500">\n'
						#			'  </a>\n'
						#			'</td>\n'
						#			''.format(G=self,fdir=formDir,comp=comp,band=band,form=form))

						
						line = line + ('\n\t<td align="center"> <!-- GASS{G.id}: {col} ({comp}-component) -->\n'
									   '\t\t<a target="_blank" href="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_{disp}.png">\n'
									   '\t\t\t<img src="{G.dir}{fdir}/GASS{G.id}{form}_{band}_{comp}comp_{disp}.png" height="300">\n'
									   '\t\t</a>\n'
									   '\t</td>\n'
									   ''.format(G=self,fdir=formDir,comp=comp,band=band,form=form,col=colNames[disp],disp=disp))

					line = line + '\n</tr>'
					
			
		if(vrb>1):print(line)
		return line

	def write_model(self,model,nComps,state):
		""" Writes the model into hmtl table rows
		<param: model [Table]> - The model data table
		<param: nComps [int]> - the number of components
		<param: state [int]> - whether the model is the initial (0) or optimised (1)


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
	parser.add_argument('-f','--forms',
						action='store',dest='forms',nargs='+',default=[''],
						help="A list of different forms of ",metavar="FORM")
	parser.add_argument('-d','--displays',
						action='store',dest='displays',type=str,nargs='*',default=['LikelihoodOptimised'],
						help="The types of plots to be displayed given as a list of filename extensions (e.g. _LikelihoodOptimised).",metavar="EXTENSION")
	parser.add_argument('-j','--judge',
						action='store',dest='judge',type=str,default=None,
						help="A judgment form will be created for each galaxy and the results sent to this file.",metavar="FILENAME")
	parser.add_argument('-o','--output',
						action='store',dest='output',type=str,default="/{home}/robincook/Documents/PhD/GASS/Optimisations/GASSFittingSummary.html".format(home=HOME),
						help="The filename of the output file",metavar="FILENAME")
	parser.add_argument('-t','--title',
						action='store',dest='title',type=str,default="GASS Galaxies Optimisation Summary",
						help="The title to be displayed on the HTML page.",metavar="TITLE")
	parser.add_argument('-l','--list',
						action='store_true',dest='list',default=False,
						help="Print available galaxies and exit.")
	parser.add_argument('-v',
						action='count',dest='verbosity',default=0,
						help="The level of verbosity to print to stdout.")

	### PARSE ARGUMENTS ###
	args = parser.parse_args()

	## Validate verbosity
	global vrb
	if (args.verbosity < 0): vrb = 0
	elif (args.verbosity > 2): vrb = 2
	else: vrb = args.verbosity

	## Get all galaxy IDs
	galData = Table.read(GALS_DATA)
	allGals = ["GASS{0}".format(ID) for ID in galData['GASS']]

	## Print list of galaxies and exit.
	if (args.list): print("\nINFO: GASS IDs for available galaxies:\n{0}\n".format(' '.join(allGals))); exit()

	## Get galaxy list
	if (args.galaxies == None and args.file == None): print("\nERROR: No '--galaxies' or '--file' specified by user!\n  -- ABORTING --\n"); exit()
	if (args.galaxies and args.file): print("\nWARNING: Both '--galaxies' and '--file' were specified. Keeping --file."); args.galaxies = None

	galNames = []
	if (args.galaxies): # If galaxy specified
		for gal in args.galaxies:
			# galName = gal if ('GASS' in gal) else 'GASS{id}'.format(id=gal) # correct for explixit inclusion of 'GASS' or not
			if (gal not in allGals): print("WARNING: Galaxy \"{0}\" does not exist.".format(gal))
			else: galNames.append(gal)

	if (args.file):
		if not os.path.exists(args.file): print("\nERROR: file \"{file}\" does not exist!\n  -- ABORTING --\n".format(file=args.file)); exit()
		galFile = open(args.file,'r')
		galNames = [gal.replace('\n','') for gal in galFile.readlines()]

		galFile.close()

	## Check if no valid galaxies given: exit
	if (len(galNames) == 0): print("ERROR: No valid galaxies given.\n  -- ABORTING --"); exit()

	## Get fit parameters
	compList = args.components
	bandList = args.bands
	formList = args.forms

	## Get the file to output judgements of the fits.
	global judge; judge = True if args.judge != None else False
	global judgeFile; judgeFile = args.judge

	### DEFINE GALAXY OBJECT LIST ###
	galList = []
	for gal in galNames:
		index = allGals.index(gal)
		galList.append(Galaxy(gal.replace('GASS',''),galData['SDSS'][index],float(galData['RA'][index]),float(galData['DEC'][index]),bandList,compList,formList))

	if(vrb>1):
		print("\nINFO: Galaxy list ({0})".format(len(galList)))
		for gal in galNames: print("  - {0}".format(gal))
		print("\nINFO: Frequency bands")
		for band in bandList: print("  - {0}".format(band))
		print("\nINFO: Components")
		for comp in compList: print("  - {0}".format(comp))
		print("\nINFO: Optimisation forms")
		for form in formList: print("  - {0}".format(form))
	

	### OPEN OUTPUT FILE ###
	filename = args.output if ('.html' in args.output) else "{0}.html".format(args.output)
	if(vrb>0):print("\nINFO: Output filename = \"{file}\"".format(file=filename))
	file = open(filename,'w')

	# Set predefined column names:
	global displays; displays = args.displays
	global colNames; colNames = {}
	for disp in displays:
		if   ("EllipseInitial" in disp): colNames[disp] = "Initial Ellipse Plot"
		elif ("EllipseOptimised" in disp): colNames[disp] = "Optimised Ellipse Plot"
		elif ("LikelihoodInitial" in disp): colNames[disp] = "Initial Likelihood Plot"
		elif ("LikelihoodOptimised" in disp): colNames[disp] = "Optimised Likelihood Plot"
		else: colNames[disp] = disp


	### WRITE HTML FILE ###
	## Write preamble.
	file.write("<!DOCTYPE html>\n"
			   "<html>\n"
			   "<head>\n"
			   "\t<title> {title} </title>\n".format(title=args.title))

	file.write("\t<style type=\"text/css\">\n"
			   "\t\ttable { page-break-inside:auto }\n"
			   "\t\ttr    { page-break-inside:avoid; page-break-after:auto }\n"
			   "\t\tthead { display:table-header-group }\n"
			   "\t\ttfoot { display:table-footer-group }\n"
			   "\t</style>\n")

	file.write("</head>\n"
			   "<body>\n"
			   "\t<h1> {title} </h1>\n"
			   "\t<h2>\n"
			   "\t\tAuthor: Robin H. Cook<br>\n"
			   "\t\tDate:&nbsp;{date}&nbsp;-&nbsp;{time}\n"
			   "\t</h2>\n".format(title=args.title,date=time.strftime("%d/%m/%Y"),time=time.strftime("%I:%M:%S")))



	## Begin TABLE
	file.write("\n\n<table bgcolor=\"white\" bordercolor=\"black\" align=\"center\" border=\"5\" cellspacing=\"1\" cellpadding=\"10\">\n")
	file.write("<tr>\n")
	file.write("\t<th><font size=\"+2.5\">Object</font></th>\n")
	file.write("\t<th><font size=\"+2.5\">Inputs</font></th>\n")
	for disp in displays:
		file.write("\t<th><font size=\"+2.5\">{name}</font></th>\n".format(name=colNames[disp]))
	
	file.write('</tr>\n\n')


	## Write each Galaxy to table row
	for gal in galList:
		file.write(gal.write_row())

	## End Table
	file.write('\n\n</table>')


	## Write closing statements
	file.write("\n</body>\n"
			   "\n</html>\n")


	file.close()


if __name__ == "__main__": main()