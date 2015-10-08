#!/usr/bin/python

import os
import shutil
# loadLongSubjectFile(fileName)
# fileName is a list of subjects
# each line contains a tab-delimited list of image names
# the first file in the list is the first timepoint, the others are other timepoints
# returns a dictionary with the first timepoint file as the key and for that key the value is a list of the remaining files in that line, i.e. the other timepoints for that subject

def loadLongSubjectFile(fileName):
	if not os.path.isfile(fileName):
		return None
	else:
		fp = open(fileName)

		mainAndTimePoints = dict()
		timePointsToReference = dict()

		for line in fp:
			splitLine = line.strip().split()
			referenceImage = splitLine[0]
			otherTimePoints = splitLine[1:]
			if referenceImage in mainAndTimePoints:
				print "Warning: reference image appears twice in file - " + referenceImage
			else:
				if len(otherTimePoints) > 0:
					mainAndTimePoints[referenceImage] = otherTimePoints[:]
					for otherImage in otherTimePoints:
						if otherImage in timePointsToReference or otherImage in mainAndTimePoints:
							print "Warning: image appears twice - " + otherImage
						else:
							timePointsToReference[otherImage] = referenceImage
				else:
					timePointsToReference[referenceImage] = None
		fp.close()
		return (mainAndTimePoints, timePointsToReference)

def longOptionText():
	return "\t\t--long-subject-timepoint-file=<file>: the file with the subject ids and timepoints, see LONGITUDINAL FILE FORMAT for details"

def longOptionFileFormat():
	return "\tLONGITUDINAL FILE FORMAT\n\t\tFile is a newline-delimited text file with each row being: <tp1> <tp2> <tp3> ... each tpX is the name of a NIFTI file in the input directory"

if __name__ == "__main__":
	mainAndTimePoints, timePointsToReference = loadLongSubjectFile('long_file.txt')
	#print mainAndTimePoints
	#print timePointsToReference
	
	inDir = 'T1NeckCroppedRegToFirstTimepoint'
	
	for curFile in timePointsToReference.iterkeys():
		if timePointsToReference[curFile] != None:

			extensions = ['.mat', '.nii.gz']
			
			for curExtension in extensions:
				curInFile = os.path.join(inDir, curFile + curExtension)
				curOutFile = os.path.join(inDir, curFile + "_to_" + timePointsToReference[curFile] + curExtension)
				if os.path.isfile(curInFile):
					print curInFile + " -> " + curOutFile
					shutil.move(curInFile, curOutFile)
