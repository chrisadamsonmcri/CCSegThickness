#!/usr/bin/python

import tempfile

import os
import scipy
import cv2

import math
import numpy
import nibabel
import pylab
import h5py
import shutil

import CCSegUtils


def segCCToNativeOneFile(segIMG, outputBase, midSagMATFile, segMATFile, dilateToParaSag = False, flirtInterpType = 'trilinear'):
	
	FID = h5py.File(segMATFile, 'r')
	
	resampledAVWShape = numpy.array(FID['resampledAVWShape'])
	LKCropRows = numpy.array(FID['LKCropRows'])
	LKCropCols = numpy.array(FID['LKCropCols'])
	midSagAVWxx = numpy.array(FID['midSagAVWxx'])
	midSagAVWyy = numpy.array(FID['midSagAVWyy'])
	resamplexx = numpy.array(FID['resamplexx'])
	resampleyy = numpy.array(FID['resampleyy'])

	#print seg['templatePixdims']
	FID.close()
	#print seg['IMG'].dtype
	
	FID = h5py.File(midSagMATFile, 'r')
	
	MSPMethod = str(numpy.array(FID['MSPMethod']))
	
	FID.close()

	midSagAVW = numpy.array(FID['midSagAVW'])
	
	skullCrop = numpy.array(FID["skullCrop"]) # the initial cropping indices of the background
	originalOrientationString = str(numpy.array(FID['originalOrientationString']))
	originalNativeFile = str(numpy.array(FID['originalNativeFile']))
	originalNativeCroppedFile = str(numpy.array(FID['originalNativeCroppedFile']))
	flirtTemplateFile = str(numpy.array(FID['flirtTemplateFile']))
	
	if flirtTemplateFile == "***":
		flirtTemplateFile = CCSegUtils.MNI152FLIRTTemplate()
	
	flirtMAT = numpy.array(FID["flirtMAT"]) # the transformation between  originalNativeCroppedFile -> flirtTemplateFile
	flirtCropZerosRows = numpy.array(FID["flirtCropZerosRows"])
	flirtCropZerosCols = numpy.array(FID["flirtCropZerosCols"])

	
	#midSagAVW = midSagAVW[:, ::-1]
	#midSagAVW = numpy.rot90(midSagAVW, -1)
	FID.close()
	
	if not numpy.array_equal(flirtMAT.shape, (4, 4)):
		print "the flirt matrix should be a 4x4 array, it is not"
		return

	del FID


	#pylab.subplot(1, 4, 1); CCSegUtils.showIMG(seg['initialSeg'])
	#pylab.subplot(1, 4, 2); CCSegUtils.showIMG(seg['finalSeg'])
	#pylab.subplot(1, 4, 3); CCSegUtils.showIMG(seg['finalSegNoArtefacts'])
	#pylab.subplot(1, 4, 4); CCSegUtils.showIMG(seg['IMG'])
	
	finalSegResampledAVWSpace = numpy.zeros(resampledAVWShape, dtype = numpy.uint8)

	#finalSegResampledAVWSpace[LKOutput['cropRows'], LKOutput['cropCols']] = finalSeg
	finalSegResampledAVWSpace[LKCropRows[0]:(LKCropRows[1] + 1), LKCropCols[0]:(LKCropCols[1] + 1)] = segIMG
	midSagAVWX, midSagAVWY = numpy.meshgrid(midSagAVWxx, midSagAVWyy)
	
	finalSegMidSagAVWSpace = CCSegUtils.interp2q(resamplexx, resampleyy, numpy.double(finalSegResampledAVWSpace), midSagAVWX, midSagAVWY, interpmethod = 'nearest', extrapval = 0)
	
	flirtTemplateNII = nibabel.load(flirtTemplateFile)
	
	#print finalSegMidSagAVWSpace.shape
	finalSegTemplateSpace = numpy.zeros((flirtTemplateNII.shape[2], flirtTemplateNII.shape[1]), dtype = numpy.uint8)
	#print finalSegTemplateSpace.shape
	finalSegTemplateSpace[flirtCropZerosRows[0]:(flirtCropZerosRows[-1] + 1), flirtCropZerosCols[0]:(flirtCropZerosCols[-1] + 1)] = finalSegMidSagAVWSpace
	
	finalSegTemplateSpaceNIIData = numpy.zeros((flirtTemplateNII.shape[1], flirtTemplateNII.shape[0], flirtTemplateNII.shape[2]), dtype = numpy.uint8)

	#NIIData = flirtTemplateNII.get_data()
	#NIIData = numpy.rot90(NIIData, 1)
	
	#print "flirt template: " + flirtTemplateFile

	T = math.floor(flirtTemplateNII.shape[0] / 2)
	
	#print flirtTemplateNII
	#print finalSegTemplateSpaceNIIData.shape
	#print finalSegTemplateSpace.shape
	finalSegTemplateSpaceNIIData[:, T] = numpy.rot90(finalSegTemplateSpace[:, ::-1], -1)
	if dilateToParaSag == True:
		finalSegTemplateSpaceNIIData[:, T - 1] = numpy.rot90(finalSegTemplateSpace[:, ::-1], -1)
		finalSegTemplateSpaceNIIData[:, T + 1] = numpy.rot90(finalSegTemplateSpace[:, ::-1], -1)

	finalSegTemplateSpaceNIIData = numpy.uint8(numpy.rot90(finalSegTemplateSpaceNIIData, -1))
	
	newNII = nibabel.Nifti1Image(finalSegTemplateSpaceNIIData, flirtTemplateNII.get_affine())
	
	NIITempDir = tempfile.mkdtemp()
	nibabel.save(newNII, os.path.join(NIITempDir, 'test.nii.gz'))
	
	invFlirtMAT = numpy.linalg.inv(flirtMAT)
	numpy.savetxt(os.path.join(NIITempDir, 'test.mat'), invFlirtMAT, fmt = '%.18f')

	# run flirt to project the segmentation to native space
	CommandString = os.environ['FSLDIR'] + '/bin/flirt -in ' + os.path.join(NIITempDir, 'test') + ' -ref ' + originalNativeCroppedFile + ' -applyxfm -init ' + os.path.join(NIITempDir, 'test.mat') + ' -interp ' + flirtInterpType + ' -out ' + (outputBase + "_native_axial_cropped.nii.gz")
	#print CommandString
	os.system(CommandString)
	
	segNativeCroppedNII = nibabel.load((outputBase + "_native_axial_cropped.nii.gz"))

	originalNativeNII = nibabel.load(originalNativeFile)
	finalSegNativeNIIData = numpy.zeros(originalNativeNII.shape, dtype = segNativeCroppedNII.get_data().dtype)
	finalSegNativeNIIData = numpy.rot90(finalSegNativeNIIData, 1)
	finalSegNativeNIIData[skullCrop[0, 0]:(skullCrop[0, 1] + 1), skullCrop[1, 0]:(skullCrop[1, 1] + 1), skullCrop[2, 0]:(skullCrop[2, 1] + 1)] = numpy.rot90(segNativeCroppedNII.get_data(), 1)
	finalSegNativeNIIData = numpy.rot90(finalSegNativeNIIData, -1)
	
	#print skullCrop
	finalSegNativeSpaceNII = nibabel.Nifti1Image(finalSegNativeNIIData, originalNativeNII.get_affine())
	
	nibabel.save(finalSegNativeSpaceNII, (outputBase + "_native_axial.nii.gz"))
		
	os.system(os.environ['FSLDIR'] + '/bin/fslswapdim ' + (outputBase + "_native_axial.nii.gz") + ' ' + originalOrientationString + ' ' + (outputBase + "_native.nii.gz"))

	shutil.rmtree(NIITempDir)

	del newNII
	del T

def segCCToNative(outputBase, doGraphics = False, dilateToParaSag = False):
	midSagMATFile = outputBase + "_midsag.hdf5"
	
	assert(os.path.isfile(midSagMATFile)),"midsag hdf5 file not found"

	segMATFile = outputBase + "_seg.hdf5"
	segManeditMATFile = outputBase + "_seg_manedit.hdf5"
	
	#finalSeg = numpy.array(seg['finalSeg'])
	
	assert(os.path.isfile(segMATFile)),"Seg mat file not found: " + segMATFile
	
	FID = h5py.File(segMATFile, 'r')
	
	autoFinalSeg = numpy.array(FID['finalSegNoArtefacts'])	
	templatePixdims = numpy.array(FID['templatePixdims'])
	#print seg['templatePixdims']
	FID.close()
	#print seg['IMG'].dtype
	

	#pylab.subplot(1, 4, 1); CCSegUtils.showIMG(seg['initialSeg'])
	#pylab.subplot(1, 4, 2); CCSegUtils.showIMG(seg['finalSeg'])
	#pylab.subplot(1, 4, 3); CCSegUtils.showIMG(seg['finalSegNoArtefacts'])
	#pylab.subplot(1, 4, 4); CCSegUtils.showIMG(seg['IMG'])
	
	
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()	
	if os.path.isfile(segManeditMATFile):
		FID = h5py.File(segManeditMATFile, 'r')
		finalSeg = numpy.array(FID['finalSegManEdit']) > 0
		FID.close()	
	else:
		finalSeg = numpy.array(autoFinalSeg)
	
	segCCToNativeOneFile(numpy.uint8(finalSeg), outputBase + "_seg", midSagMATFile, segMATFile, dilateToParaSag = dilateToParaSag)
	#print finalSeg.shape
	thicknessMATFile = outputBase + "_thickness.hdf5"
	
	if os.path.isfile(thicknessMATFile):
		
		finalSegxx = numpy.arange(finalSeg.shape[1])
		finalSegyy = numpy.arange(finalSeg.shape[0])
		
		finalSegX, finalSegY = numpy.meshgrid(finalSegxx, finalSegyy)

		FID = h5py.File(thicknessMATFile, 'r')

		labelNames = ['witelson', 'hoferFrahm', 'emsell']

		for curLabels in labelNames:
			curLabelIMG = numpy.array(FID[curLabels + 'Labels'])
			# for the parcellations, we need to get the parcellation images in the same space as the segmentation images
			# so we resample the parcellation image using the grid of the segmentation image, we need to reverse a scaling that was performed in the thickness part:
			curLabelIMG = CCSegUtils.interp2q(numpy.array(FID['xx']) / templatePixdims[0], numpy.array(FID['yy']) / templatePixdims[1], numpy.double(curLabelIMG), finalSegX, finalSegY, interpmethod = 'nearest', extrapval = 0)
			
			#SR = 1
			#SC = 3
			#pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(finalSeg)
			#pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(curLabelIMG)
			#pylab.subplot(SR, SC, 3); CCSegUtils.showIMG(T)
			#pylab.gcf().set_size_inches((20, 10), forward = True)
			#pylab.show()
			#quit()	
			segCCToNativeOneFile(numpy.uint8(curLabelIMG), outputBase + "_parcellation_" + curLabels.lower(), midSagMATFile, segMATFile, flirtInterpType = 'nearestneighbour', dilateToParaSag = dilateToParaSag)
			
		FID.close()

