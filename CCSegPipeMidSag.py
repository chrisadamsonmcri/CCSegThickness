#!/usr/bin/env python3

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

import matplotlib.pyplot as plt
import pylab
import CCSegPipeMidSagSymmetric


#import nifti

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

# apply an orientation transformation according to the ornt
# ornt is a [3, 2] array
# the format is
# [newx, flipx] 
# [newy, flipy] 
# [newz, flipz] 
# each row, the index is the new column
def applyOrntToNIIAffine(NII, ornt_transform):
	NIIAffine = NII.get_affine()
	
	# use fsl's method for 
	# make a transformation affine matrix
	transformAffine = numpy.zeros_like(NIIAffine)
	transformAffine[3, 3] = 1

	for curDim in range(3):
		newDim = int(ornt_transform[curDim, 0])
		transformAffine[curDim, newDim] = ornt_transform[curDim, 1]
		#print str(curDim) + " " + str(newDim)
		if ornt_transform[curDim, 1] < 0:
			transformAffine[curDim, 3] = (NII.shape[newDim] - 1) * NII.header.get_zooms()[newDim]

	pixDimsVector = numpy.concatenate((numpy.array(NII.header.get_zooms()), [1]))

	return numpy.matrix(NIIAffine) * numpy.diag(1.0 / pixDimsVector) * numpy.matrix(transformAffine) * numpy.diag(pixDimsVector)
	
def midsagExtract(inputBase, outputBase, MSPMethod, doGraphics = False, skullStripped = False, useMNIBrainMask = True):
	
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
	NIIShape = NII.shape
	assert(len(NIIShape) == 2 or len(NIIShape) == 3),"The input NIFTI file is not 2D or 3D: " + inputBase
	
	if len(NIIShape) == 2:
		
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
		print("2D input not supported yet")
		quit()
		# 2D image, so it is already the midsagittal plane
	else:
		# 3D image

		outputMAT = outputBase + "_midsag.hdf5"
		(head, tail) = os.path.split(outputBase)
		subjectID = tail
		if doGraphics:
			PNGDirectory = os.path.join(head, "midsag")
			CCSegUtils.mkdirSafe(PNGDirectory)

		del head; del tail;
		# use FSLORIENT to get the RADIOLOGICAL/NEUROLOGICAL orientation of the image
		# we may not need this
		
		NIIPixdims = NII.get_header().get_zooms()[1:3]
		NIIAffine = NII.get_affine()
		#print NII.shape

		if MSPMethod == 'long':
			# testing
			(head, subjectID) = os.path.split(outputBase)
			stdMat = os.path.join(head, subjectID + "_to_std.mat")
			assert(os.path.isfile(stdMat)),"FLIRT MAT file not found, need to run CCSegLongPreprocess: " + stdMat
			
			flirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate()
			flirtTemplateFileBrainMask = flirtTemplateFile[:-7] + "_brain_mask.nii.gz"
			
			# use the to std file if it is already there
			toStdNII = CCSegUtils.findNIFTIFromPrefix(os.path.join(head, subjectID + "_to_std"))
			if toStdNII != None:
				print(("using already standard file: " + toStdNII))
				deleteStdNII = False
			else:
				NIITempDir = tempfile.mkdtemp()
				toStdNII = os.path.join(NIITempDir, subjectID + "_to_std.nii.gz")
				flirtInterp = 'trilinear'

				commandString = [os.path.join(os.environ['FSLDIR'], 'bin', 'flirt'), '-interp', flirtInterp, '-applyxfm', '-init', stdMat, '-ref', flirtTemplateFile, '-out', toStdNII, '-in', NIFTIFileName]
				#print " ".join(commandString)
				subprocess.call(commandString)
				deleteStdNII = True
				
			#shutil.copy(toStdNII, os.path.join(head, subjectID + "_to_std.nii.gz"))
			flirtMAT = numpy.loadtxt(stdMat)

			#toStdNII = CCSegUtils.findNIFTIFromPrefix(os.path.join(head, subjectID + "_to_std"))
			#assert(toStdNII != None),"standard image not found, need to run CCSegLongPreprocess: " + os.path.join(head, subjectID + "_to_std")
			NIIBrainMask = nibabel.load(flirtTemplateFileBrainMask)

			NII = nibabel.load(toStdNII)
			NIIPixdims = NII.get_header().get_zooms()[1:3]
			NIIShape = NII.shape
			
			NIIData = numpy.rot90(NII.get_data(), 1)
			NIIBrainMaskData = numpy.rot90(NIIBrainMask.get_data(), 1)

			T = int(math.floor(NIIData.shape[1] / 2))
			# extract the midsagittal slice
			if (NIIData.shape[1] % 2 == 0):
				# even number of slices
				midSagAVW = (numpy.double(NIIData[:, T - 1]) + numpy.double(NIIData[:, T])) / 2.0
				midSagAVWBrainMask = (numpy.double(NIIBrainMaskData[:, T - 1]) + numpy.double(NIIBrainMaskData[:, T])) / 2.0
			else:
				# odd number of slices
				midSagAVW = numpy.double(NIIData[:, T])
				midSagAVWBrainMask = numpy.double(NIIBrainMaskData[:, T])
			pylab.clf()
			pylab.imshow(midSagAVWBrainMask)
			pylab.show()
	#print NIIPixdims
	#	
			midSagAVW[numpy.where(midSagAVWBrainMask < 0.5)] = numpy.nan

			midSagAVW = numpy.rot90(midSagAVW, 1)
			midSagAVW = numpy.array(midSagAVW[:, ::-1])
			
			if deleteStdNII == True:
				shutil.rmtree(NIITempDir)

		else: # not if MSPMethod == 'long':
			
			NIITempDir = tempfile.mkdtemp()
			# nibabel reorientation 
			NIIAXCodes = nibabel.aff2axcodes(NII.get_affine())
			NIIOrnt = nibabel.orientations.axcodes2ornt(NIIAXCodes)

			#print NIIAXCodes
			#print NIIOrnt
				
			MNITemplate = CCSegUtils.MNI152FLIRTTemplate(skullStripped = skullStripped)
			MNITemplateNII = nibabel.load(MNITemplate)
			MNIAXCodes = nibabel.aff2axcodes(MNITemplateNII.get_affine())
			MNIOrnt = nibabel.orientations.axcodes2ornt(MNIAXCodes)
			#print MNIAXCodes
			#print "Affine of template"
			#print MNITemplateNII.get_affine()
			# gets the transformation from start_ornt to end_ornt
			NIIToMNITransformOrnt = nibabel.orientations.ornt_transform(NIIOrnt, MNIOrnt)
			#print "Transform ornt"
			#print NIIToMNITransformOrnt
			
			axialNIIAffine = applyOrntToNIIAffine(NII, NIIToMNITransformOrnt)

			# transforms the image array
			axialNIIIMG = nibabel.orientations.apply_orientation(NII.get_data(), NIIToMNITransformOrnt)
			#rint NII.shape
			#rint D.shape
			#rint "Affine of original image"
			#rint NII.get_affine()

			axialNII = nibabel.Nifti1Image(axialNIIIMG, affine = axialNIIAffine)
			nibabel.save(axialNII, outputBase + "_native.nii.gz")
			
			axialNIIShape = axialNII.shape

			#axialNIIIMG = numpy.rot90(axialNIIIMG, 1)
			axialNIIPixdims = axialNII.get_header().get_zooms()
			#print NIIPixdims
			#print NewNII.header.get_zooms()
			#quit()
			#print "pixdims: " + str(NIIPixdims)
			#print NII
			# perform 3-class otsu thresholding
			#print numpy.min(NIIData)
			#print numpy.max(NIIData)
			### NECK CROPPING ###
			OtsuSeg = Otsu.robustOtsu(axialNIIIMG, [0.02, 0.98], NumberClasses = 2)

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
			
			# find the slice 180mm inferior of the top of the scalp
			minK = numpy.int32(numpy.maximum(numpy.floor(maxK - 180 / axialNIIPixdims[2]), 0))
			
			skullCrop = numpy.array([[minI, maxI], [minJ, maxJ], [minK, maxK]])

			#print str(minI) + " " + str(maxI) + ", " + str(minJ) + " " + str(maxJ) + ", " + str(minK) + " " + str(maxK)
				
			axialNIIIMG = numpy.take(axialNIIIMG, numpy.arange(minI, maxI + 1), axis=0)
			axialNIIIMG = numpy.take(axialNIIIMG, numpy.arange(minJ, maxJ + 1), axis=1)
			axialNIIIMG = numpy.take(axialNIIIMG, numpy.arange(minK, maxK + 1), axis=2)
			### NECK CROPPING ###

			# scale the image if the intensities are too high
			if numpy.max(axialNIIIMG) > 32767:
				axialNIIIMG = numpy.double(axialNIIIMG)
				axialNIIIMG = (axialNIIIMG - numpy.min(axialNIIIMG)) / (numpy.max(axialNIIIMG) - numpy.min(axialNIIIMG))
				axialNIIIMG = numpy.round(axialNIIIMG * 1000.0)
			
			axialNIIIMG = numpy.int16(axialNIIIMG)
			#axialNIIIMG = numpy.rot90(axialNIIIMG, -1)

			NIISaving = nibabel.Nifti1Image(axialNIIIMG, axialNII.get_affine(), axialNII.get_header())
			axialCroppedNIIShape = NIISaving.shape
			#print "axialNIIShape"
			#print axialNIIShape
			NIISaving.set_data_dtype(numpy.int16)
			
			if not MSPMethod == 'symmetric':
				NIIFileForART = os.path.join(NIITempDir, 'in_art.nii')
				NIIFileForARTOutput = os.path.join(NIITempDir, 'in_art_output.nii')
				nibabel.save(NIISaving, NIIFileForART)
			nibabel.save(NIISaving, outputBase + "_native_neckcropped.nii.gz")
			#del NIISaving
			#rint NIISaving.header.get_zooms()
			if MSPMethod == 'symmetric':
				# the 
#				realFlirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate()
#				NIIFileARTOutputAffineMat = os.path.join(NIITempDir, 'in_art_output.mat')
#				flirtCost = 'mutualinfo'
#				flirtInterp = 'trilinear'
#
#				# perform a neck crop
#				flirtROTDegrees = 15
#				CMD = [os.environ['FSLDIR'] + '/bin/flirt', 
#				'-in', NIIFileForART, 
#				'-ref', realFlirtTemplateFile, 
#				'-dof', str(12),
#				'-searchrx', str(-flirtROTDegrees), str(flirtROTDegrees),
#				'-searchry', str(-flirtROTDegrees), str(flirtROTDegrees),
#				'-searchrz', str(-flirtROTDegrees), str(flirtROTDegrees),
#				'-omat', NIIFileARTOutputAffineMat, '-cost', flirtCost, "-interp", flirtInterp]
#				subprocess.call(CMD)	
#				NeckCropAffineMAT = numpy.loadtxt(NIIFileARTOutputAffineMat)
#				NeckCropAffineMATINV = numpy.linalg.inv(NeckCropAffineMAT)
#
#				shutil.rmtree(NIITempDir)
#				quit()
				OtsuSeg = Otsu.robustOtsu(axialNIIIMG, [0.02, 0.98], NumberClasses = 2)
				print("Midsagittal extraction, symmetry method iterations")				
				IMG = axialNIIIMG * OtsuSeg
				I = numpy.where(OtsuSeg)
				
				IMGCOG = numpy.array([numpy.mean(I[0]), numpy.mean(I[1]), numpy.mean(I[2])]) * axialNIIPixdims
				del I
				
				IMGxx = numpy.arange(axialNIIIMG.shape[0]) * NIISaving.header.get_zooms()[0] - IMGCOG[0]
				IMGyy = numpy.arange(axialNIIIMG.shape[1]) * NIISaving.header.get_zooms()[1] - IMGCOG[1]
				IMGzz = numpy.arange(axialNIIIMG.shape[2]) * NIISaving.header.get_zooms()[2] - IMGCOG[2]
				
				SIMGxx = IMGxx[::2]
				SIMGyy = IMGyy[::2]
				SIMGzz = IMGzz[::2]
				# downsample the image 

				# smooth it
				SIMG = scipy.ndimage.filters.gaussian_filter(IMG, 2.0 / numpy.array(axialNIIPixdims), mode = 'constant', cval = 0)
				
				SIMG = SIMG[::2, ::2, ::2]
				#print SIMG.shape
				#print IMG.shape
				#CCSegUtils.imshow(SIMG[:, :, 50])
				#3plt.show()
				#IMGCoords = numpy.meshgrid(IMGxx, IMGyy, IMGzz, indexing = 'ij')
				
				#print IMGxx.shape
				#print IMGyy.shape
				#print IMGzz.shape
				# optimise the parameters in an interleaved fashion

				paramRange = 10

				curRotY = 0
				curRotZ = 0
				curTransX = 0
				
				#lastParams = numpy.array([curRotY, curRotZ, curTransX])
				lastCost = CCSegPipeMidSagSymmetric.corrCoefCost(SIMG)
				for curIter in range(5):
					initialTransX = numpy.arange(curTransX - paramRange * axialNIIPixdims[0], curTransX + (paramRange + 1) * axialNIIPixdims[0], axialNIIPixdims[0])
					initialTransXCosts = numpy.zeros(initialTransX.size)
					initialTransXCosts.fill(-numpy.inf) 
					# save time, evaluate every second one
					for z in range(0, initialTransX.size, 2):
						initialTransXCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, curRotZ, initialTransX[z])
					# then do the neighbours of the best one thus far
					bestSoFar = numpy.argmax(initialTransXCosts)
					if bestSoFar > 0:
						initialTransXCosts[bestSoFar - 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, curRotZ, initialTransX[bestSoFar - 1])
					if bestSoFar < initialTransX.size - 1:
						initialTransXCosts[bestSoFar + 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, curRotZ, initialTransX[bestSoFar + 1])

					#for z in range(initialTransX.size):
						#print "transx "  + str(initialTransX[z])
					#	initialTransXCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, curRotZ, initialTransX[z])
					
					# positive values for translation will push the image up
					# negative values for translation will push the image down
					bestInitialTransXIDX = numpy.argmax(initialTransXCosts)
					curTransX = initialTransX[bestInitialTransXIDX]
					#print "Best TransX: " + str(curTransX)
					if doGraphics:
						pylab.clf()
						TIMG = CCSegPipeMidSagSymmetric.transformIMG(IMG, IMGxx, IMGyy, IMGzz, curRotY, curRotZ, curTransX) 
						CCSegPipeMidSagSymmetric.showMidSag(TIMG)
						pylab.gcf().set_size_inches((20, 10), forward = True)
						outputPNG = os.path.join(PNGDirectory, subjectID + "_sym_iter_" + str(curIter + 1) + "_1.png")
						pylab.savefig(outputPNG)
						CCSegUtils.cropAutoWhitePNG(outputPNG)

	
					initialRotZ = numpy.arange(curRotZ - paramRange, curRotZ + paramRange + 1)
					initialRotZCosts = numpy.zeros(initialRotZ.size)
					initialRotZCosts.fill(-numpy.inf)
					for z in range(0, initialRotZ.size, 2):
						initialRotZCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, initialRotZ[z], curTransX)

					#for z in range(initialRotZ.size):
					#	initialRotZCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, initialRotZ[z], curTransX)
					bestSoFar = numpy.argmax(initialRotZCosts)
					if bestSoFar > 0:
						initialRotZCosts[bestSoFar - 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, initialRotZ[bestSoFar - 1], curTransX)
					if bestSoFar < initialRotZ.size - 1:
						initialRotZCosts[bestSoFar + 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, initialRotZ[bestSoFar + 1], curTransX)
						
					#plt.plot(initialRotZ, initialRotZCosts)
						
					bestInitialRotZIDX = numpy.argmax(initialRotZCosts)
					curRotZ = initialRotZ[bestInitialRotZIDX]
					#print "Best initial RotZ: " + str(curRotZ)
					if doGraphics:
						pylab.clf()
						TIMG = CCSegPipeMidSagSymmetric.transformIMG(IMG, IMGxx, IMGyy, IMGzz, curRotY, curRotZ, curTransX) 
						CCSegPipeMidSagSymmetric.showMidSag(TIMG)
						pylab.gcf().set_size_inches((20, 10), forward = True)
						outputPNG = os.path.join(PNGDirectory, subjectID + "_sym_iter_" + str(curIter + 1) + "_2.png")
						pylab.savefig(outputPNG)
						CCSegUtils.cropAutoWhitePNG(outputPNG)

	
					initialRotY = numpy.arange(curRotY - paramRange, curRotY + paramRange + 1)
					initialRotYCosts = numpy.zeros(initialRotY.size)
					initialRotYCosts.fill(-numpy.inf)
					
					for z in range(0, initialRotY.size, 2):
						initialRotYCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, initialRotY[z], curRotZ, curTransX)

					#for z in range(initialRotY.size):
					#	initialRotYCosts[z] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, initialRotY[z], curRotZ, curTransX)
					
					#plt.plot(initialRotY, initialRotYCosts)
					bestSoFar = numpy.argmax(initialRotYCosts)
					
					if bestSoFar > 0:
						initialRotYCosts[bestSoFar - 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, initialRotY[bestSoFar - 1], curRotZ, curTransX)
					if bestSoFar < initialRotY.size - 1:
						initialRotYCosts[bestSoFar + 1] = CCSegPipeMidSagSymmetric.transformCost(SIMG, SIMGxx, SIMGyy, SIMGzz, initialRotY[bestSoFar + 1], curRotZ, curTransX)
					
					bestInitialRotYIDX = numpy.argmax(initialRotYCosts)
					curRotY = initialRotY[bestInitialRotYIDX]
					#print "Best initial RotY: " + str(curRotY)
					#if doGraphics:
					#	pylab.clf()
					#	TIMG = CCSegPipeMidSagSymmetric.transformIMG(IMG, IMGxx, IMGyy, IMGzz, curRotY, curRotZ, curTransX) 
					#	CCSegPipeMidSagSymmetric.showMidSag(TIMG)
					#	pylab.gcf().set_size_inches((20, 10), forward = True)
					#	outputPNG = os.path.join(PNGDirectory, subjectID + "_sym_iter_" + str(curIter + 1) + "_3.png")
					#	pylab.savefig(outputPNG)
					#	CCSegUtils.cropAutoWhitePNG(outputPNG)

					TIMG = CCSegPipeMidSagSymmetric.transformIMG(SIMG, SIMGxx, SIMGyy, SIMGzz, curRotY, curRotZ, curTransX) 
#					SR = 2
#					SC = 2
#					pylab.subplot(SR, SC, 1)
#					CCSegUtils.imshow(TIMG[:, :, TIMG.shape[2] // 2])
#					plt.plot([0, TIMG.shape[1]], [TIMG.shape[0] / 2, TIMG.shape[0] / 2], 'w-')
#					pylab.subplot(SR, SC, 2)
#					T = numpy.rot90(TIMG[TIMG.shape[0] // 2], 1)
#					CCSegUtils.imshow(T)
#					pylab.subplot(SR, SC, 3)
#					T = numpy.rot90(TIMG[:, TIMG.shape[1] // 2], 1)
#					CCSegUtils.imshow(T)
#					plt.plot([T.shape[1] / 2, T.shape[1] / 2], [0, TIMG.shape[0]], 'w-')
#					#plt.plot([0, TIMG.shape[1]], [TIMG.shape[0] / 2, TIMG.shape[0] / 2], 'w-')
#					print initialRotYCosts[bestInitialRotYIDX]
#					plt.show()
					print(("iteration " + str(curIter) + " params: " + str(curRotY) + " " + str(curRotZ) + " " + str(curTransX)))
					
					if initialRotYCosts[bestInitialRotYIDX] == lastCost:
						lastCost = initialRotYCosts[bestInitialRotYIDX]
						del TIMG
						del initialRotY
						del initialRotZ
						del initialTransX
						del bestInitialRotYIDX
						del bestInitialRotZIDX
						del bestInitialTransXIDX
						del SIMG
						break
					lastCost = initialRotYCosts[bestInitialRotYIDX]
				
				print(("final params: " + str(curRotY) + " " + str(curRotZ) + " " + str(curTransX)))
				
				#MGxx = numpy.arange(axialNIIIMG.shape[0]) * NIISaving.header.get_zooms()[0] - IMGCOG[0]
				#MGyy = numpy.arange(axialNIIIMG.shape[1]) * NIISaving.header.get_zooms()[1] - IMGCOG[1]
				#MGzz = numpy.arange(axialNIIIMG.shape[2]) * NIISaving.header.get_zooms()[2] - IMGCOG[2]
				
				TIMG = CCSegPipeMidSagSymmetric.transformIMG(IMG, IMGxx, IMGyy, IMGzz, curRotY, curRotZ, curTransX) 
				#CCSegPipeMidSagSymmetric.showMidSag(TIMG)
				#plt.show()
				NIISaving = nibabel.Nifti1Image(TIMG, axialNII.get_affine(), axialNII.get_header())
				nibabel.save(NIISaving, outputBase + "_native_midsag_sym.nii.gz")
				flirtCost = 'mutualinfo'
				flirtInterp = 'trilinear'
				
				#NIIFileARTOutputAffineMat = os.path.join(NIITempDir, 'in_art_output.mat')
				
				# register to our 2mm MNI template, project brain mask from template space back to native to perform
				# a cropping to remove neck
				flirtROTDegrees = 15
				
				if skullStripped == False: #$not os.path.isfile(outputBase + "_native_midsag_sym_to_template.mat"):
					CMD = [os.environ['FSLDIR'] + '/bin/flirt', 
					'-in', outputBase + "_native_midsag_sym", 
					'-out', outputBase + "_native_midsag_sym_to_template", 
					'-ref', CCSegUtils.MNI152FLIRTSymNeckCropTemplate(),
					'-dof', str(12),
					'-searchrx', str(-flirtROTDegrees), str(flirtROTDegrees),
					'-searchry', str(0), str(0),
					'-searchrz', str(0), str(0),
					'-omat', outputBase + "_native_midsag_sym_to_template.mat",
					'-cost', flirtCost,
					'-interp', flirtInterp]
					subprocess.call(CMD)
					
					T = numpy.loadtxt(outputBase + "_native_midsag_sym_to_template.mat")
					invT = numpy.linalg.inv(T)
					numpy.savetxt(outputBase + "_native_midsag_template_to_sym.mat", invT)

					CMD = [os.environ['FSLDIR'] + '/bin/flirt', 
					'-in', CCSegUtils.MNI152FLIRTSymNeckCropTemplate() + "_brain_mask", 
					'-out', outputBase + "_MNI_brain_mask", 
					'-ref', outputBase + "_native_midsag_sym",
					'-applyxfm',
					'-init', outputBase + "_native_midsag_template_to_sym.mat",
					'-interp', 'nearestneighbour',
					'-datatype', 'char']
					subprocess.call(CMD)

					brainMaskNII = nibabel.load(outputBase + "_MNI_brain_mask.nii.gz")
					I = numpy.where(brainMaskNII.get_data() > 0)
				#CCSegPipeMidSagSymmetric.showMidSag(TIMG)
				#brainMaskBoundingBox = 	
				#brainMaskCrop = numpy.array([[numpy.min(I[0]), numpy.max(I[0])], [numpy.min(I[1]), numpy.max(I[1])], [numpy.min(I[2]), numpy.max(I[2])]])
				midSlice = int(math.floor(TIMG.shape[0] / 2))
				
				if skullStripped == False: #$not os.path.isfile(outputBase + "_native_midsag_sym_to_template.mat"):
					if (TIMG.shape[0] % 2 == 0):
						midSagBrainMask = (numpy.double(brainMaskNII.get_data()[midSlice - 1, :, :]) + numpy.double(brainMaskNII.get_data()[midSlice, :, :])) / 2.0
					else:
						midSagBrainMask = numpy.double(brainMaskNII.get_data()[midSlice, :, :])

				# extract the midsagittal slice
				if (TIMG.shape[0] % 2 == 0):
					# even number of slices
					midSagAVW = (numpy.double(TIMG[midSlice - 1, :, :]) + numpy.double(TIMG[midSlice, :, :])) / 2.0
				else:
					# odd number of slices
					midSagAVW = numpy.double(TIMG[midSlice, :, :])
				
				midSagAVW = numpy.rot90(midSagAVW, 1)
				
				if skullStripped == False:

					midSagBrainMask = numpy.rot90(midSagBrainMask, 1)
					midSagAVW[numpy.logical_not(midSagBrainMask > 0)] = 0
				#parasagittalSlices, parasagittalFX, parasagittalFY, parasagittalFZ = CCSegUtils.parasagittalSlicesAndGradients(TIMG, axialNIIPixdims, numSlices = 3)

				#pylab.subplot(1, 2, 1)
				#pylab.imshow(midSagAVW)
				#pylab.subplot(1, 2, 2)
				#pylab.imshow(midSagBrainMask)
				#pylab.show()
				#quit()
				#midSagAVW = numpy.array(midSagAVW[:, ::-1])
				if doGraphics:
					pylab.clf()
					CCSegPipeMidSagSymmetric.showMidSag(TIMG, midSagAVW)
					pylab.gcf().set_size_inches((20, 10), forward = True)
					outputPNG = os.path.join(PNGDirectory, subjectID + ".png")
					pylab.savefig(outputPNG)
					CCSegUtils.cropAutoWhitePNG(outputPNG)

			elif MSPMethod == 'acpcdetect':
				scriptPath = os.path.realpath(__file__)
				(head, tail) = os.path.split(scriptPath)

				CommandString = 'ARTHOME=' + os.path.join(head, 'ART') + " " + os.path.join(head, 'ART', 'acpcdetect') + ' -i ' + NIIFileForART + ' -o ' + NIIFileForARTOutput
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
				flirtDOF = 7

			elif MSPMethod == 'flirt_affine':
				
				flirtTemplateFile = '***'
				flirtOutputFile = NIIFileForARTOutput
				flirtDOF = 12
			
			if not MSPMethod == 'symmetric':
				flirtCost = 'mutualinfo'
				flirtInterp = 'trilinear'

				NIIFileARTOutputAffineMat = os.path.join(NIITempDir, 'in_art_output.mat')

				if flirtTemplateFile == "***":
					realFlirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate(skullStripped = skullStripped)
				else:
					realFlirtTemplateFile = flirtTemplateFile
				
				flirtROTDegrees = 15
				CMD = [os.environ['FSLDIR'] + '/bin/flirt', 
				'-in', NIIFileForART, 
				'-ref', realFlirtTemplateFile, 
				'-dof', str(flirtDOF),
				'-searchrx', str(-flirtROTDegrees), str(flirtROTDegrees),
				'-searchry', str(-flirtROTDegrees), str(flirtROTDegrees),
				'-searchrz', str(-flirtROTDegrees), str(flirtROTDegrees),
				'-omat', NIIFileARTOutputAffineMat, '-cost', flirtCost, "-interp", flirtInterp]

				if flirtOutputFile != None:
					#CommandString = CommandString + " -out " + NIIFileForARTOutput 
					CMD.extend(["-out", NIIFileForARTOutput])
				print((" ".join(CMD)))
				subprocess.call(CMD)
				
				shutil.copyfile(NIIFileForARTOutput + ".gz", outputBase + "_native_cropped_to_template.nii.gz")

				#flirtMAT = open(NIIFileARTOutputAffineMat, 'r')
				flirtMAT = numpy.loadtxt(NIIFileARTOutputAffineMat)
				shutil.copyfile(NIIFileARTOutputAffineMat, outputBase + "_native_cropped_to_template.mat")

				# find out whether the output file is a nifti or compressed nifti
				#print NIIFileForARTOutput
				NII = nibabel.load(NIIFileForARTOutput + ".gz")
			
				NIIData = NII.get_data()
				NIIData = numpy.rot90(NIIData, 1)
				
				T = int(math.floor(NIIData.shape[1] / 2))
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

	if MSPMethod == 'long':
		FID = h5py.File(outputMAT, 'w')
		FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
		FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
		FID.create_dataset("MSPMethod", data=MSPMethod)
		FID.create_dataset("flirtMAT", data=flirtMAT)
		FID.create_dataset("flirtTemplateFile", data=flirtTemplateFile)
		FID.create_dataset("originalNativeFile", data=NIFTIFileName)

		FID.close()
	elif MSPMethod == 'symmetric':
		FID = h5py.File(outputMAT, 'w')

		FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
		FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
		FID.create_dataset("MSPMethod", data=MSPMethod)
		FID.create_dataset("skullCrop", data=skullCrop)
		FID.create_dataset("NIIOrnt", data=NIIOrnt)
		FID.create_dataset("NIIAffine", data=NIIAffine)
		FID.create_dataset("NIIShape", data=NIIShape)
		FID.create_dataset("axialNIIAffine", data=axialNIIAffine)
		FID.create_dataset("axialNIIShape", data=axialNIIShape)
		FID.create_dataset("axialCroppedNIIShape", data=axialCroppedNIIShape)
		FID.create_dataset("midSlice", data=midSlice)
		FID.create_dataset("IMGxx", data=IMGxx)
		FID.create_dataset("IMGyy", data=IMGyy)
		FID.create_dataset("IMGzz", data=IMGzz)
		FID.create_dataset("finalRotY", data=curRotY)
		FID.create_dataset("finalRotZ", data=curRotZ)
		FID.create_dataset("finalTransX", data=curTransX)
		
		#FID.create_dataset("parasagittalSlices", data=parasagittalSlices)
		#FID.create_dataset("parasagittalFX", data=parasagittalFX)
		#FID.create_dataset("parasagittalFY", data=parasagittalFY)
		#FID.create_dataset("parasagittalFZ", data=parasagittalFZ)
		FID.close()

	else:
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
		#print NII.get_header().get_zooms()
		#print NIIPixdims
		#ylab.imshow(midSagAVW)
		#pylab.set_cmap(pylab.cm.gray)
		#pylab.show()

	#if len(NIIShape) == 2:
	
	#if doGraphics:
	#	T = numpy.double(midSagAVW)
	#	T[numpy.isnan(T)] = 0
#
#		T = (T - numpy.min(T)) / (numpy.max(T) - numpy.min(T))
#		T = numpy.uint8(numpy.round(T * 255))
#		outputPNG = os.path.join(PNGDirectory, subjectID + "_midsag.png")
#			
#		scipy.misc.imsave(outputPNG, T)
#		del T

#def midsagExtract(inputBase, outputBase, MSPMethod):


