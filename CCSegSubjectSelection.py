#!/usr/bin/python

import os

import CCSegUtils

def criticalError(errorString):
	print "Error: " + str(errorString)
	quit()

def printHelpText():
	print "SUBJECT SELECTION - use these options to select which images in <input directory> are processed"
	print "Each of these options take an option <file> that point to a file containing a newline delimited list of image prefixes"
	print
	print "\t--single-subject=<nifti file prefix>: only this subject will be processed"
	print "\t--subjects-include=<comma separated list>: comma separated list of subjects to process"
	print "\t--subjects-include-file=<file>: file name that contains a newline-delimited list of subjects to process"
	print "\t--subjects-exclude=<comma separated list>: comma separated list of subjects to exclude"
	print "\t--subjects-exclude-file=<file>: file name that contains a newline-delimited list of subjects to exclude"

def getOptions():
	return ['subjects-include=', 'subjects-exclude=', 'subjects-include-file=', 'subjects-exclude-file=', 'single-subject=']
		
def selectSubjectsFromOpts(inputDirectory, opts = None):
	
	singleSubject = None
	subjectsIncludeFile = None
	subjectsExcludeFile = None
	subjectsIncludeListArg = None
	subjectsExcludeListArg = None
	
	for o, a in opts:
		if o == '--subjects-include':
			subjectsIncludeListArg = a.split(',')
		if o == '--subjects-exclude':
			subjectsExcludeListArg = a.split(',')
		if o == '--subjects-include-file':
			subjectsIncludeFile = a
		if o == '--subjects-exclude-file':
			subjectsExcludeFile = a
		if o == '--single-subject':
			singleSubject = a
	
	#print "inputDirectory"
	#print inputDirectory
	#print "subjectsIncludeFile"
	#print subjectsIncludeFile
	#print "subjectsExcludeFile"
	#print subjectsExcludeFile
	#print "subjectsIncludeListArg"
	#print subjectsIncludeListArg
	#print "subjectsExcludeListArg"
	#print subjectsExcludeListArg
	#print "singleSubject"
	#print singleSubject

	subjectsInclude = None
	if subjectsIncludeFile != None:
		print "using include file: " + subjectsIncludeFile
		if not os.path.isfile(subjectsIncludeFile):
			criticalError("The subjects include file was given, but the file was not found")
		fp = open(subjectsIncludeFile, 'r')
		subjectsInclude = fp.readlines()
		subjectsInclude = map(lambda s: s.strip(), subjectsInclude)
		fp.close()
	
	subjectsExclude = None
	if subjectsExcludeFile != None:
		print "using exclude file: " + subjectsExcludeFile
		if not os.path.isfile(subjectsExcludeFile):
			criticalError("The subjects exclude file was given, but the file was not found")
		fp = open(subjectsExcludeFile, 'r')
		subjectsExclude = fp.readlines()
		subjectsExclude = map(lambda s: s.strip(), subjectsExclude)
		fp.close()
	
	if subjectsIncludeListArg != None:
		if subjectsInclude != None:
			for curSubject in subjectsIncludeListArg:
				if curSubject in subjectsInclude:
					subjectsInclude.append(curSubject)
		else:
			subjectsInclude = subjectsIncludeListArg[:]
	
	if subjectsExcludeListArg != None:
		if subjectsExclude != None:
			for curSubject in subjectsExcludeListArg:
				if not curSubject in subjectsExclude:
					subjectsExclude.append(curSubject)
		else:
			subjectsExclude = subjectsExcludeListArg[:]
	
	
	#print "subjectsInclude = "
	#print subjectsInclude
	#print "subjectsExclude = "
	#print subjectsExclude

	# if we are using an include AND exclude lists, remove all elements of the include list that are in the exclude list
	if subjectsInclude != None and subjectsExclude != None:
		print "Warning, using both include and exclude lists, removing excluded subjects from include list"
		subjectsInclude = [x for x in subjectsInclude if not x in subjectsExclude]

	# get a list of nifti files in the input directory
	inputDirectoryFileList = os.listdir(inputDirectory)
	
	# contains a list of nifti file prefixes for all nii.gz files found in inputDirectory
	#inputNIFTIFiles = [x for x in sorted(inputDirectoryFileList) if x.endswith(".nii.gz")]
	# remove the .nii.gz 
	#inputNIFTIFiles = map(lambda s: s[:-7], inputNIFTIFiles)
	
	inputDirectoryFileList = [os.path.join(inputDirectory, x) for x in inputDirectoryFileList]
	inputNIFTIFiles = CCSegUtils.imglob(inputDirectoryFileList)
	# remove None elements
	inputNIFTIFiles = [x for x in inputNIFTIFiles if x != None]
	#print "inputNIFTIFiles = "
	#print inputNIFTIFiles
		
	# strip the input directory
	for z in range(len(inputNIFTIFiles)):
		(head, tail) = os.path.split(inputNIFTIFiles[z])
		inputNIFTIFiles[z] = tail
		del tail
		del head
	
	#inputNIFTIFiles = map(lambda s: s[:-7], inputNIFTIFiles)

	# replace with CCSegUtils.imglob
	
	if subjectsInclude != None:
		inputNIFTIFiles = [x for x in inputNIFTIFiles if x in subjectsInclude]
	
	if subjectsExclude != None:
		inputNIFTIFiles = [x for x in inputNIFTIFiles if x not in subjectsExclude]
	
	if singleSubject != None:
		if singleSubject in inputNIFTIFiles:
			actualInputFiles = [singleSubject]
			print "processing single subject: " + singleSubject
		else:
			criticalError("Single subject was given as: " + singleSubject + ", but the NIFTI file was not found in input directory")
	else:
		actualInputFiles = inputNIFTIFiles[:]
	
	return sorted(actualInputFiles)
	
