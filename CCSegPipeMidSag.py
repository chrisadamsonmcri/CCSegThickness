#!/usr/bin/python

import tempfile
import os
import numpy
import nibabel
import scipy
import math
import h5py
import shutil
import re
import gzip
import subprocess

import errno
import time

import CCSegUtils
import Otsu

import numpy.linalg

# (IOrientation, JOrientation, KOrientation) = nifti_orientations(SForm)
# 
# DESCRIPTION
#	Returns the most likely orientations for the X, Y, Z axes given a
#	transformation matrix
# niftiio/nifti_io.c
# 1772    Input:  4x4 matrix that transforms (i,j,k) indexes to (x,y,z) coordinates,
# 1773            where +x=Right, +y=Anterior, +z=Superior.
# 1774            (Only the upper-left 3x3 corner of R is used herein.)
# 1775    Output: 3 orientation codes that correspond to the closest "standard"
# 1776            anatomical orientation of the (i,j,k) axes.
# 1777    Method: Find which permutation of (x,y,z) has the smallest angle to the
# 1778            (i,j,k) axes directions, which are the columns of the R matrix.
# 1779    Errors: The codes returned will be zero.
# 1780 
# 1781    For example, an axial volume might get return values of
# 1782      IOrientation = NIFTI_R2L   (i axis is mostly Right to Left)
# 1783      JOrientation = NIFTI_P2A   (j axis is mostly Posterior to Anterior)
# 1784      KOrientation = NIFTI_I2S   (k axis is mostly Inferior to Superior)
# 1785    </pre>
# 1786 

def niftiCodeToString(Code):
	if Code == 1:
		return 'NIFTI_L2R';
	elif Code == -1:
		return 'NIFTI_R2L';
	elif Code == 2:
		return 'NIFTI_P2A';
	elif Code == -2:
		return 'NIFTI_A2P';
	elif Code == 3:
		return 'NIFTI_I2S';
	elif Code == -3:
		return 'NIFTI_S2I';
	else:
		return 'NIFTI_INVALID'

def niftiOrientations(SForm):
	
	import numpy.linalg
	#print SForm
	xi = SForm[0, 0]; xj = SForm[0, 1]; xk = SForm[0, 2];
	yi = SForm[1, 0]; yj = SForm[1, 1]; yk = SForm[1, 2];
	zi = SForm[2, 0]; zj = SForm[2, 1]; zk = SForm[2, 2];
	#print str(xi) + " " + str(xj) + " " + str(xk)
	#print str(yi) + " " + str(yj) + " " + str(yk)
	#print str(zi) + " " + str(zj) + " " + str(zk)
	# normalise i axis
	val = numpy.sqrt(xi * xi + yi * yi + zi * zi);
	if val == 0:
		return ('NIFTI_INVALID', 'NIFTI_INVALID', 'NIFTI_INVALID');

	xi = xi / val; yi = yi / val; zi = zi / val;
	
	val = numpy.sqrt(xj * xj + yj * yj + zj * zj);
	if val == 0:
		return ('NIFTI_INVALID', 'NIFTI_INVALID', 'NIFTI_INVALID');

	xj = xj / val; yj = yj / val; zj = zj / val;

	# orthogonalize j axis to i axis, if needed
	val = xi * xj + yi * yj + zi * zj;
	if (numpy.abs(val) > 1e-4):
		xj = xj - val * xi; yj = yj - val * yi; zj = zj - val * zi;
#		
#		% renormalise j axis
		val = numpy.sqrt(xj * xj + yj * yj + zj * zj);
		if val == 0:
			return ('NIFTI_INVALID', 'NIFTI_INVALID', 'NIFTI_INVALID');

		xj = xj / val; yj = yj / val; zj = zj / val;
	#if (numpy.abs(val) > 1e-4)

	# normalise k axis
	val = numpy.sqrt(xk * xk + yk * yk + zk * zk);
	#	if it is zero, make it the cross product between i and j
	if(val == 0):
		xk = yi * zj - zi * yj; yk = zi * xj - zj * xi; zk = xi * yj - yi * xj;
	else:
		xk = xk / val; yk = yk / val; zk = zk / val;

	# orthogonalize k axis to j axis, if needed
	val = xj * xk + yj * yk + zj * zk;
	if(numpy.abs(val) > 1e-4):
		xk = xk - val * xj; yk = yk - val * yj; zk = zk - val * zj;
		
		# renormalise k axis
		val = numpy.sqrt(xk * xk + yk * yk + zk * zk);
		if(val == 0):
			return ('NIFTI_INVALID', 'NIFTI_INVALID', 'NIFTI_INVALID');
		xk = xk / val; yk = yk / val; zk = zk / val;
	#if(numpy.abs(val) > 1e-4)
	
	# form Q matrix
	Q = numpy.matrix([[xi, xj, xk], [yi, yj, yk], [zi, zj, zk]]);
	#print Q
	DetQ = numpy.linalg.det(Q);
	if(DetQ == 0):
		return ('NIFTI_INVALID', 'NIFTI_INVALID', 'NIFTI_INVALID');

	vbest = -666.0;
	ibest = 1; pbest = 1; qbest = 1; rbest = 1;
	jbest = 2;
	kbest = 3;

	for ColI in range(3):
		for ColJ in range(3):
			if not ColJ == ColI:
				for ColK in range(3):
					if not (ColI == ColK or ColJ == ColK):
						P = numpy.matrix(numpy.zeros([3, 3]))
						for ColP in [-1, 1]:
							for ColQ in [-1, 1]:
								for ColR in [-1, 1]:
									P[0, ColI] = ColP; P[1, ColJ] = ColQ; P[2, ColK] = ColR;
									DetP = numpy.linalg.det(P)
									if DetP * DetQ <= 0:
										continue;
									M = P * Q;
									val = numpy.trace(M)
									if val > vbest:
										vbest = val;
										ibest = ColI; jbest = ColJ; kbest = ColK;
										pbest = ColP; qbest = ColQ; rbest = ColR;
										#print str(vbest) + " " + str(ibest) + " " + str(jbest) + " " + str(kbest) + " " + str(pbest) + " " + str(qbest) + " " + str(rbest)
									#if val > vbest:
								#for ColR in [-1, 1]:
							#for ColQ in [-1, 1]:
						#for ColP in [-1, 1]:
					#if not (ColI == ColK or ColJ == ColK):
				#for ColK in range(3):
			#if not ColJ == ColI:
		#for ColJ in range(3):
	#for ColI in r3yyange(3):
	
	IOrientation = niftiCodeToString((ibest + 1) * pbest)
	JOrientation = niftiCodeToString((jbest + 1) * qbest)
	KOrientation = niftiCodeToString((kbest + 1) * rbest)
	
	return (IOrientation, JOrientation, KOrientation)

def midsagExtract(inputBase, outputBase, MSPMethod, doGraphics = False):
	
	# check for the presence of FSLDIR in the enviornment variables
	assert('FSLDIR' in os.environ),"FSLDIR not set, FSL must be set up"
	
	NIFTIFileName = CCSegUtils.findNIFTIFromPrefix(inputBase)

	if NIFTIFileName == None:
		print("Something went wrong, the NIFTI file doesnt exist")
		quit()
	NII = nibabel.load(NIFTIFileName)
	#print NII
	#print inputBase
	# find out if we have a 2D or a 3D image
	NIISize = NII.get_shape()
	assert(len(NIISize) == 2 or len(NIISize) == 3),"The input NIFTI file is not 2D or 3D: " + inputBase
	
	if len(NIISize) == 2:
		
		#FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
		#FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
		#FID.create_dataset("MSPMethod", data=MSPMethod)
		#FID.create_dataset("originalOrientationString", data=origOrientationString)
		#FID.create_dataset("originalNativeFile", data=(outputBase + "_native.nii.gz"))
		#FID.create_dataset("skullCrop", data=skullCrop)
		#FID.create_dataset("originalNativeCroppedFile", data=(outputBase + "_native_cropped.nii.gz"))
		#FID.create_dataset("flirtMAT", data=flirtMAT)
		#FID.create_dataset("flirtTemplateFile", data=flirtTemplateFile)
		#FID.create_dataset("flirtCropZerosRows", data=flirtCropZerosRows)
		#FID.create_dataset("flirtCropZerosCols", data=flirtCropZerosCols)
		print "2D input not supported yet"
		quit()
		# 2D image, so it is already the midsagittal plane
	else:
		# 3D image

		outputMAT = outputBase + "_midsag.hdf5"
		(head, tail) = os.path.split(outputBase)
		if doGraphics:
			PNGDirectory = os.path.join(head, "midsag")
			try:
				os.makedirs(PNGDirectory)
			except OSError as exc: # Python >2.5
				if exc.errno == errno.EEXIST and os.path.isdir(PNGDirectory):
					pass
				else:
					raise Exception

			outputPNG = os.path.join(PNGDirectory, tail + "_midsag.png")
		del head; del tail;
		# use FSLORIENT to get the RADIOLOGICAL/NEUROLOGICAL orientation of the image
		# we may not need this
		
		NIIPixdims = NII.get_header().get_zooms()[1:3]
		if MSPMethod == 'long':
			# testing
			(head, subjectID) = os.path.split(outputBase)
			stdMat = os.path.join(head, subjectID + "_to_std.mat")
			assert(os.path.isfile(stdMat)),"FLIRT MAT file not found, need to run CCSegLongPreprocess: " + stdMat
			
			flirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate()
			flirtTemplateFileBrainMask = flirtTemplateFile[:-7] + "_brain_mask.nii.gz"

			NIITempDir = tempfile.mkdtemp()
			toStdNII = os.path.join(NIITempDir, subjectID + "_to_std.nii.gz")
			flirtInterp = 'trilinear'

			commandString = [os.path.join(os.environ['FSLDIR'], 'bin', 'flirt'), '-interp', flirtInterp, '-applyxfm', '-init', stdMat, '-ref', flirtTemplateFile, '-out', toStdNII, '-in', NIFTIFileName]
			#print " ".join(commandString)
			subprocess.call(commandString)

			flirtMAT = numpy.loadtxt(stdMat)

			#toStdNII = CCSegUtils.findNIFTIFromPrefix(os.path.join(head, subjectID + "_to_std"))
			#assert(toStdNII != None),"standard image not found, need to run CCSegLongPreprocess: " + os.path.join(head, subjectID + "_to_std")
			NIIBrainMask = nibabel.load(flirtTemplateFileBrainMask)

			NII = nibabel.load(toStdNII)
			NIIPixdims = NII.get_header().get_zooms()[1:3]
			NIISize = NII.get_shape()
			
			NIIData = numpy.rot90(NII.get_data(), 1)
			NIIBrainMaskData = numpy.rot90(NIIBrainMask.get_data(), 1)

			T = math.floor(NIIData.shape[1] / 2)
			# extract the midsagittal slice
			if (NIIData.shape[1] % 2 == 0):
				# even number of slices
				midSagAVW = (numpy.double(NIIData[:, T - 1]) + numpy.double(NIIData[:, T])) / 2.0
				midSagAVWBrainMask = (numpy.double(NIIBrainMaskData[:, T - 1]) + numpy.double(NIIBrainMaskData[:, T])) / 2.0
			else:
				# odd number of slices
				midSagAVW = numpy.double(NIIData[:, T])
				midSagAVWBrainMask = numpy.double(NIIBrainMaskData[:, T])
			
			midSagAVW[numpy.where(midSagAVWBrainMask < 0.5)] = numpy.nan

			midSagAVW = numpy.rot90(midSagAVW, 1)
			midSagAVW = numpy.array(midSagAVW[:, ::-1])
			
			shutil.rmtree(NIITempDir)

		else:
			Transform = NII.get_sform(coded=True)
			if Transform[0] == None:
				print "Trying qform in " + inputBase
				Transform = NII.get_qform(coded=True)
				assert(Transform[0] != None),"No transformation information in NIFTI file" + inputBase
			det = numpy.linalg.det(Transform[0][0:3, 0:3])
			(icode, jcode, kcode) = niftiOrientations(Transform[0])

			del Transform

			if numpy.abs(det) < 1e-12:
				OrientationString = None
			elif det < 0.0:
				OrientationString = "RL PA IS"
			else:
				OrientationString = "LR PA IS"
			
			origOrientationString = list()

			pat = re.compile("NIFTI_([LRASIP])2([LRASIP])")
			mat = pat.match(icode)
			if mat != None:
				origOrientationString.append(mat.group(1) + mat.group(2))
			mat = pat.match(jcode)
			if mat != None:
				origOrientationString.append(mat.group(1) + mat.group(2))
			mat = pat.match(kcode)
			if mat != None:
				origOrientationString.append(mat.group(1) + mat.group(2))
			
			origOrientationString = " ".join(origOrientationString)
			#print icode + " " + jcode + " " + kcode
			#assert(not (icode == "NIFTI_INVALID" or jcode == "NIFTI_INVALID" or kcode == "NIFTI_INVALID")),"Invalid NIFTI orientations in " + inputBase

#		if icode == "NIFTI_L2R" or jcode == "NIFTI_L2R" or kcode == "NIFTI_L2R":
#			OrientationString = "LR PA IS"
#		elif icode == "NIFTI_R2L" or jcode == "NIFTI_R2L" or kcode == "NIFTI_R2L":
#			OrientationString = "RL PA IS"
#		else:
#			OrientationString = None
#		del icode
#		del jcode
#		del kcode
#
			assert(OrientationString != None),"Could not find RADIOLOGICAL/NEUROLOGICAL orientation in " + inputBase
			
			#print "OrientationString = " + OrientationString

			NIITempDir = tempfile.mkdtemp()
			#print NIITempDir
			NIIFileForMidSag = os.path.join(NIITempDir, 'in.nii.gz')

			os.system(os.environ['FSLDIR'] + '/bin/fslswapdim ' + inputBase + ' ' + OrientationString + ' ' + NIIFileForMidSag)
			
			NII = nibabel.load(NIIFileForMidSag)
			shutil.copyfile(NIIFileForMidSag, outputBase + "_native.nii.gz")
			NIIData = numpy.array(NII.get_data())
			NIIData = numpy.rot90(NIIData, 1)
			NIIPixdims = NII.get_header().get_zooms()
			#print "pixdims: " + str(NIIPixdims)
			#print NII
			# perform 3-class otsu thresholding
			#print numpy.min(NIIData)
			#print numpy.max(NIIData)
			OtsuSeg = Otsu.robustOtsu(NIIData, [0.02, 0.98], NumberClasses = 2)

			L, NumLabels = scipy.ndimage.measurements.label(OtsuSeg, structure = numpy.ones([3, 3, 3]))
			LAreas = scipy.ndimage.measurements.labeled_comprehension(OtsuSeg, L, numpy.arange(1, NumLabels + 1), numpy.size, numpy.uint32, 0)
			MaxLabel = numpy.argmax(LAreas) + 1
				
			OtsuSeg[numpy.where(L != MaxLabel)] = 0
			del L; del MaxLabel; del NumLabels
			
			# cut off the neck, find the bounding box
			I = numpy.nonzero(OtsuSeg)

			minI = numpy.min(I[0]); maxI = numpy.max(I[0]);
			minJ = numpy.min(I[1]); maxJ = numpy.max(I[1]);
			minK = numpy.min(I[2]);	maxK = numpy.max(I[2]);
			
			minK = numpy.int32(numpy.maximum(numpy.floor(maxK - 180 / NIIPixdims[2]), 0))
			
			skullCrop = numpy.array([[minI, maxI], [minJ, maxJ], [minK, maxK]])

			#print str(minI) + " " + str(maxI) + ", " + str(minJ) + " " + str(maxJ) + ", " + str(minK) + " " + str(maxK)
			
			NIIData = numpy.take(NIIData, numpy.arange(minI, maxI + 1), axis=0)
			NIIData = numpy.take(NIIData, numpy.arange(minJ, maxJ + 1), axis=1)
			NIIData = numpy.take(NIIData, numpy.arange(minK, maxK + 1), axis=2)

			if numpy.max(NIIData) > 32767:
				NIIData = numpy.double(NIIData)
				NIIData = (NIIData - numpy.min(NIIData)) / (numpy.max(NIIData) - numpy.min(NIIData))
				NIIData = numpy.round(NIIData * 1000.0)
			
			NIIData = numpy.int16(NIIData)

			NIIData = numpy.rot90(NIIData, -1)

			NIISaving = nibabel.Nifti1Image(NIIData, NII.get_affine(), NII.get_header())
			NIISaving.set_data_dtype(numpy.int16)	

			NIIFileForART = os.path.join(NIITempDir, 'in_art.nii')
			NIIFileForARTOutput = os.path.join(NIITempDir, 'in_art_output.nii')
			nibabel.save(NIISaving, NIIFileForART)
			nibabel.save(NIISaving, outputBase + "_native_cropped.nii.gz")
			del NIISaving

			if MSPMethod == 'acpcdetect':
				scriptPath = os.path.realpath(__file__)
				(head, tail) = os.path.split(scriptPath)

				CommandString = 'ARTHOME=' + os.path.join(head, 'ART') + " " + os.path.join(head, 'ART', 'acpcdetect') + ' -i ' + NIIFileForART + ' -o ' + NIIFileForARTOutput
				print CommandString
				os.system(CommandString)

				inF = file(NIIFileForARTOutput, 'rb')
				s = inF.read()
				inF.close()

				outF = gzip.GzipFile(NIIFileForARTOutput + ".gz", 'wb')
				outF.write(s)
				outF.close()
				
				os.unlink(NIIFileForARTOutput)
				
				flirtTemplateFile = outputBase + "_template.nii.gz"
				shutil.copyfile(NIIFileForARTOutput + ".gz", flirtTemplateFile)
				
				# get the aligned image and register the original image to it to get the transformation
				flirtOutputFile = None

			elif MSPMethod == 'flirt':
				
				flirtTemplateFile = '***'
				flirtOutputFile = NIIFileForARTOutput

			flirtCost = 'mutualinfo'
			flirtInterp = 'trilinear'

			NIIFileARTOutputAffineMat = os.path.join(NIITempDir, 'in_art_output.mat')

			D = 15
			if flirtTemplateFile == "***":
				realFlirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate()
			else:
				realFlirtTemplateFile = flirtTemplateFile

			CommandString = os.environ['FSLDIR'] + '/bin/flirt -in ' + NIIFileForART + ' -ref ' + realFlirtTemplateFile + ' -dof 6 -searchrx ' + str(-D) + ' ' + str(D) + ' -searchry ' + str(-D) + ' ' + str(D) + ' -searchrz ' + str(-D) + ' ' + str(D) + ' -omat ' + NIIFileARTOutputAffineMat + ' -cost ' + flirtCost + " -interp " + flirtInterp

#		interpTypes = ['trilinear', 'nearestneighbour']
#		costFns = ['mutualinfo','corratio','normcorr','normmi','leastsq']
#
#		for flirtInterp in interpTypes:
#			for flirtCost in costFns:
#				CommandString = os.environ['FSLDIR'] + '/bin/flirt -in ' + NIIFileForART + ' -ref ' + flirtTemplateFile + ' -dof 6 -searchrx ' + str(-D) + ' ' + str(D) + ' -searchry ' + str(-D) + ' ' + str(D) + ' -searchrz ' + str(-D) + ' ' + str(D) + ' -omat ' + NIIFileARTOutputAffineMat + ' -cost ' + flirtCost + ' -interp ' + flirtInterp
#				
#				start_time = time.time()
#				os.system(CommandString)
#				elapsed_time = time.time() - start_time
#				print "interp: " + flirtInterp + ", cost: " + flirtCost + ", time: " + str(elapsed_time)	
#		quit()
			if flirtOutputFile != None:
				CommandString = CommandString + " -out " + NIIFileForARTOutput 
			del D
			#print CommandString
			os.system(CommandString)

			shutil.copyfile(NIIFileForARTOutput + ".gz", outputBase + "_native_cropped_to_template.nii.gz")

			#flirtMAT = open(NIIFileARTOutputAffineMat, 'r')
			flirtMAT = numpy.loadtxt(NIIFileARTOutputAffineMat)
			shutil.copyfile(NIIFileARTOutputAffineMat, outputBase + "_native_cropped_to_template.mat")

			
			# find out whether the output file is a nifti or compressed nifti
			#print NIIFileForARTOutput
			NII = nibabel.load(NIIFileForARTOutput + ".gz")
		
			NIIData = NII.get_data()
			NIIData = numpy.rot90(NIIData, 1)
			
			T = math.floor(NIIData.shape[1] / 2)
			# extract the midsagittal slice
			if (NIIData.shape[1] % 2 == 0):
				# even number of slices
				midSagAVW = (numpy.double(NIIData[:, T - 1]) + numpy.double(NIIData[:, T])) / 2.0
			else:
				# odd number of slices
				midSagAVW = numpy.double(NIIData[:, T])
			midSagAVW = numpy.rot90(midSagAVW, 1)
			midSagAVW = numpy.array(midSagAVW[:, ::-1])

			# crop out the zeros, sometimes FLIRT shrinks the image, causes problems downstream with the registration
			I = numpy.nonzero(midSagAVW)

			flirtCropZerosRows = numpy.arange(numpy.min(I[0]), numpy.max(I[0]) + 1)
			flirtCropZerosCols = numpy.arange(numpy.min(I[1]), numpy.max(I[1]) + 1)

			midSagAVW = numpy.take(midSagAVW, flirtCropZerosRows, axis=0)
			midSagAVW = numpy.take(midSagAVW, flirtCropZerosCols, axis=1)
			#print "flirtCropZerosRows: " + str(flirtCropZerosRows[0]) + " " + str(flirtCropZerosRows[-1])
			#print "flirtCropZerosCols: " + str(flirtCropZerosCols[0]) + " " + str(flirtCropZerosCols[-1])
			del I
			shutil.rmtree(NIITempDir)

	if MSPMethod != 'long':
		FID = h5py.File(outputMAT, 'w')

		FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
		FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
		FID.create_dataset("MSPMethod", data=MSPMethod)
		FID.create_dataset("originalOrientationString", data=origOrientationString)
		FID.create_dataset("originalNativeFile", data=(outputBase + "_native.nii.gz"))
		FID.create_dataset("skullCrop", data=skullCrop)
		FID.create_dataset("originalNativeCroppedFile", data=(outputBase + "_native_cropped.nii.gz"))
		FID.create_dataset("flirtMAT", data=flirtMAT)
		FID.create_dataset("flirtTemplateFile", data=flirtTemplateFile)
		FID.create_dataset("flirtCropZerosRows", data=flirtCropZerosRows)
		FID.create_dataset("flirtCropZerosCols", data=flirtCropZerosCols)
		FID.close()
	else:
		FID = h5py.File(outputMAT, 'w')

		FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
		FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
		FID.create_dataset("MSPMethod", data=MSPMethod)
		FID.create_dataset("flirtMAT", data=flirtMAT)
		FID.create_dataset("flirtTemplateFile", data=flirtTemplateFile)
		FID.close()

		#print NII.get_header().get_zooms()
		#print NIIPixdims
		#ylab.imshow(midSagAVW)
		#pylab.set_cmap(pylab.cm.gray)
		#pylab.show()

	#if len(NIISize) == 2:
	
	if doGraphics:
		T = numpy.double(midSagAVW)
		T = (T - numpy.min(T)) / (numpy.max(T) - numpy.min(T))
		T = numpy.uint8(numpy.round(T * 255))
			
		scipy.misc.imsave(outputPNG, T)
		del T

#def midsagExtract(inputBase, outputBase, MSPMethod):


