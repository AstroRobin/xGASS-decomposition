""" GetSDSSFile.py
A script to automatically download SDSS images and files using wget commands.

Author: Robin H. Cook
Date: 09/05/17
"""

from argparse import ArgumentParser
import os
import time
from astropy.table import Table

import subprocess
HOME = str(subprocess.check_output("pwd")).split('/')[1]

#[RUN]-[CAMCOL]-[RERUN]-[FIELD]


def main():

	parser = ArgumentParser()
	parser.add_argument('-r','--run',
						action='store',dest='run',nargs='+',type=int,default=None,
						help="The run number(s) for the frame.",metavar="RUN")
	parser.add_argument('-c','--camcol',
						action='store',dest='camcol',nargs='+',type=int,default=None,
						help="The camcol number(s) used for the frame.",metavar="CAMCOL")
	parser.add_argument('-R','--rerun',
						action='store',dest='rerun',nargs='+',type=int,default=None,
						help="The rerun number(s) for the frame.",metavar="RERUN")
	parser.add_argument('-f','--field',
						action='store',dest='field',nargs='+',type=int,default=None,
						help="The field number(s) for the frame.",metavar="FIELD")
	parser.add_argument('-i','--input',
						action='store',dest='input',type=str,default=None,
						help="The path to the input file listing the run, camcol, rerun and frame numbers.",metavar="INPUT")
	parser.add_argument('-b','--bands',
						action='store',dest='bands',nargs='*',type=str,default='r',
						help="The frequency band(s) to be included [all|u|g|r|i|z]",metavar="BANDS")
	parser.add_argument('-t','--types',
						action='store',dest='types',nargs='*',type=str,default='D',
						help="The files types to be downloaded [D (drC)| P (psField) | T (tsObj)].",metavar="TYPES")
	parser.add_argument('-o','--output',
						action='store',dest='output',nargs='+',type=str,default=None,
						help="The output filename(s) to be saved to.",metavar="OUTPUT")
	parser.add_argument('-g','--galaxy',
						action='store',dest='galaxy',nargs='+',type=str,default=None,
						help="The galaxy ID(s) to be used in the output filenames.",metavar="GALAXY")
	parser.add_argument('-d','--directory',
						action='store',dest='directory',type=str,default="/{home}/robincook/Google Drive/PhD/SDSS/Data/Frames".format(home=HOME),
						help="Specifies the base directory to store SDSS files in. ( \"<directory>/<run>-<camcol>-<field>/\" )",metavar="DIR")
	parser.add_argument('-l','--log',
						action='store',dest='logFile',type=str,default=None,
						help="Specifies the location of a log file to append errors.",metavar="LOGFILE")
	parser.add_argument('-v',
						action='count',dest='verbosity',default=0,
						help="The level of verbosity to print to stdout.")

	args = parser.parse_args()

	# Validate verbosity
	global vrb
	if (args.verbosity < 0): vrb = 0
	elif (args.verbosity > 2): vrb = 2
	else: vrb = args.verbosity

	if (args.bands[0].lower() == 'all'): bandList = ['u','g','r','i','z']
	else:
		bandList = []
		for f in args.bands:
			if (f in ['u','g','r','i','z']):
				bandList.append(f)
			else:
				print("WARNING: Filter \"{0}\" is not valid.".format(f))


	if (bandList == []): print("ERROR: No valid filters given!\n -- ABORTING --"); exit()
	if(vrb>0): print("INFO: Filters =     {0}".format(', '.join(bandList)))

	if (args.input != None):
		if not os.path.exists(args.input): print("\nERROR: File \"{0}\" does not exist!\n  -- ABORTING --\n".format(args.input)); exit()
		data = Table.read(args.input)
		IDList     = [kk for kk in data['name']]
		runList    = [int(kk) for kk in data['run']]
		camcolList = [int(kk) for kk in data['camcol']]
		rerunList  = [int(kk) for kk in data['rerun']]
		fieldList  = [int(kk) for kk in data['field']]
	else:
		# Parse Run values
		if (args.run == None): print("\nERROR: No \"run\" value(s) specified!\n  -- ABORTING --"); exit()
		runList = args.run
		
		# Parse camcol values
		if (args.camcol == None): print("\nERROR: No \"camcol\" value(s) specified!\n  -- ABORTING --"); exit()
		if (len(args.camcol) == len(args.run)):
			camcolList = []
			for num in args.camcol:
				if (num>0 and num<=6): camcolList.append(num)
				else: print("ERROR: camcol number \"{0}\" not valid.\n -- ABORTING --".format(num)); exit()
		else: print("ERROR: Size of camcol list ({0}) must be equal to size of run list ({1}).\n -- ABORTING --".format(len(args.camcol),len(args.run))); exit()
		
		# Parse rerun values
		if (args.rerun == None): rerunList = [40]*len(runList)
		else:
			if (len(args.rerun) == len(args.run)): rerunList = args.rerun
			else: print("ERROR: Size of rerun list ({0}) must be equal to size of run list ({1}).\n -- ABORTING --".format(len(args.rerun),len(args.run))); exit()

		# Parse field values
		if (args.field == None): print("\nERROR: No \"field\" value(s) specified!\n  -- ABORTING --"); exit()
		if (len(args.field) == len(args.run)): fieldList = args.field
		else: print("ERROR: Size of field list ({0}) must be equal to size of run list ({1}).\n -- ABORTING --".format(len(args.field),len(args.run))); exit()

	typeList = []
	for ii in range(len(args.types)):
		if (args.types[ii].lower() in ['d','p','t']): typeList.append(args.types[ii].lower())
		else: print("WARNING: File type code \"{0}\" not recognised. Choose from [D (drC)| P (psField) | T (tsObj)]".format(args.types[ii].lower()))

	if (len(typeList) == 0): print("\nERROR: No valid file types given.\n  -- ABORTING --"); exit()

	if(vrb>0):print("INFO: Run list =    {0}".format(', '.join([str(kk) for kk in runList])))
	if(vrb>0):print("INFO: Camcol list = {0}".format(', '.join([str(kk) for kk in camcolList])))
	if(vrb>0):print("INFO: Rerun list =  {0}".format(', '.join([str(kk) for kk in rerunList])))
	if(vrb>0):print("INFO: Field list =  {0}".format(', '.join([str(kk) for kk in fieldList])))
	if(vrb>0):print("INFO: File types =  {0}".format(', '.join([str(kk) for kk in typeList])))
	
	# Parse output filenames
	if (args.galaxy and args.output):
		print("WARNING: Both \"galaxy\" and \"output\" arguments specified! Using \"output\" only.")
		args.galaxy == None

	# Validate output filename list sizes
	if (args.galaxy):
		if (len(args.galaxy) == len(args.run)): galList = args.galaxy
		else: print("ERROR: Size of galaxy list ({0}) must be equal to size of run list ({1}).\n -- ABORTING --".format(len(args.galaxy),len(args.run))); exit()
	elif (args.output):
		if (len(typeList)>1): print("\nERROR: Cannot specify output filenames when multiple file types given.\n  -- ABORTING --")
		if (len(args.output) == len(args.run)): outputList = args.output
		else: print("ERROR: Size of output list ({0}) must be equal to size of run list ({1}).\n -- ABORTING --".format(len(args.output),len(args.run))); exit()


	# Validate Directory for outputs
	if (args.directory):
		if (os.path.exists(args.directory)): directory = args.directory
		else: print("ERROR: Directory \"{0}\" does not exist!\n  -- ABORTING --".format(args.directory)); exit()

	# Validate Log file
	if (args.logFile):
		if (not os.path.exists(args.logFile)):
			print("\nWARNING: Log file \"{file}\" does not exists! Printing to stdout instead.\n".format(file=args.logFile))
			args.logFile = None
		else: errFile = open(args.logFile,'a')
	
	"""
	if(vrb>0):"INFO: Making {0} directories.".format(len(galList)*len(bandList))
							for gal in galList:
								os.system('mkdir {0}'.format(gal))
								for band in bandList:
									os.system('mkdir {0}/{1}'.format(gal,band))
	"""

	#wget "http://das.sdss.org/cgi-bin/drC?RUN=1043&RERUN=40&CAMCOL=6&FIELD=190&FILTER=r"
	#wget "http://das.sdss.org/raw/5360/40/objcs/2/psField-005360-2-0098.fit"
	#wget "http://das.sdss.org/raw/1035/40/calibChunks/6/tsObj-001035-6-40-0154.fit"

	errorList = []

	for ii in range(len(runList)):
		if not os.path.exists("{dir}/{run:06d}-{camcol}-{field:04d}".format(dir=directory,run=runList[ii],camcol=camcolList[ii],field=fieldList[ii])):
			os.system("mkdir \"{dir}/{run:06d}-{camcol}-{field:04d}\"".format(dir=directory,run=runList[ii],camcol=camcolList[ii],field=fieldList[ii]))
		for jj in range(len(typeList)):
			if (typeList[jj] == 'd'):
				for kk in range(len(bandList)):
					cmd = "wget \"http://das.sdss.org/cgi-bin/drC?RUN={run}&RERUN={rerun}&CAMCOL={camcol}&FIELD={field}&FILTER={band}\"".format(run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii],band=bandList[kk])
					if   (args.output): cmd = cmd + " -O \"{0}\"".format(outputList[ii].replace(".fit","_"+bandList[kk]+".fit"))
					else: cmd = cmd + " -O \"{dir}/{run:06d}-{camcol}-{field:04d}/drC-{run:06d}-{band}{camcol}-{field:04d}.fits\"".format(dir=directory,band=bandList[kk],run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii])
					if (vrb>0): cmd = cmd + " --verbose"

					print(cmd)
					os.system(cmd)

					if not os.path.exists("{dir}/{run:06d}-{camcol}-{field:04d}/drC-{run:06d}-{band}{camcol}-{field:04d}.fits".format(dir=directory,band=bandList[kk],run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii])):
						errorList.append("drC-{run:06d}-{band}{camcol}-{field:04d}.fits".format(dir=directory,band=bandList[kk],run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii]))

			elif(typeList[jj] == 'p'):
				cmd = "wget \"http://das.sdss.org/raw/{run}/{rerun}/objcs/{camcol}/psField-{run:06d}-{camcol}-{field:04d}.fit\"".format(run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii])
				if   (args.output): cmd = cmd + " -O \"{0}\"".format(outputList[ii].replace(".fit","_"+bandList[kk]+"_PSF.fits"))
				else: cmd = cmd + " -O \"{dir}/{run:06d}-{camcol}-{field:04d}/psField-{run:06d}-{camcol}-{field:04d}.fits\"".format(dir=directory,run=runList[ii],camcol=camcolList[ii],field=fieldList[ii])
				if (vrb>0): cmd = cmd + " --verbose"

				os.system(cmd)

				if not os.path.exists("{dir}/{run:06d}-{camcol}-{field:04d}/psField-{run:06d}-{camcol}-{field:04d}.fits".format(dir=directory,run=runList[ii],camcol=camcolList[ii],field=fieldList[ii])):
					errorList.append("psField-{run:06d}-{camcol}-{field:04d}.fits".format(dir=directory,run=runList[ii],camcol=camcolList[ii],field=fieldList[ii]))

			elif(typeList[jj] == 't'):
				cmd = "wget \"http://das.sdss.org/raw/{run}/{rerun}/calibChunks/{camcol}/tsObj-{run:06d}-{camcol}-{rerun}-{field:04d}.fit\"".format(run=runList[ii],rerun=rerunList[ii],camcol=camcolList[ii],field=fieldList[ii])
				if   (args.output): cmd = cmd + " -O \"{0}\"".format(outputList[ii].replace(".fit","_"+bandList[kk]+"_Obj.fits"))
				else: cmd = cmd + " -O \"{dir}/{run:06d}-{camcol}-{field:04d}/tsObj-{run:06d}-{camcol}-{rerun}-{field:04d}.fits\"".format(dir=directory,run=runList[ii],camcol=camcolList[ii],rerun=rerunList[ii],field=fieldList[ii])
				if (vrb>0): cmd = cmd + " --verbose"

				os.system(cmd)

				if not os.path.exists("{dir}/{run:06d}-{camcol}-{field:04d}/tsObj-{run:06d}-{camcol}-{rerun}-{field:04d}.fits".format(dir=directory,run=runList[ii],camcol=camcolList[ii],rerun=rerunList[ii],field=fieldList[ii])):
					errorList.append("tsObj-{run:06d}-{camcol}-{rerun}-{field:04d}.fits".format(dir=directory,run=runList[ii],camcol=camcolList[ii],rerun=rerunList[ii],field=fieldList[ii]))					


	if (len(errorList) > 0):
		if (args.logFile != None):
			errFile.write('\n({0}  {1})\n'.format(time.strftime("%x"),time.strftime("%X")))
			errFile.write('** Files not able to be downloaded from SDSS DR7 DAS: **\n')
			for err in errorList: errFile.write("{0}\n".format(err))
		else:
			print("\n\nINFO: Error list:\n")
			for err in errorList: print(err)
	
	if (args.logFile != None): errFile.close()		


if __name__ == "__main__": main()