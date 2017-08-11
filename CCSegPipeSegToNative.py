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
import CCSegPipeMidSagSymmetric
import subprocess

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
	
	if not MSPMethod == 'symmetric':
		raise

	midSagAVW = numpy.array(FID['midSagAVW'])
	skullCrop = numpy.array(FID["skullCrop"]) # the initial cropping indices of the background
	NIIOrnt = numpy.int16(numpy.array(FID["NIIOrnt"]))
	NIIAffine = numpy.array(FID["NIIAffine"])
	NIIShape = numpy.int16(numpy.array(FID["NIIShape"]))
	axialNIIShape = numpy.int16(numpy.array(FID["axialNIIShape"]))
	axialNIIAffine = numpy.array(FID["axialNIIAffine"])
	axialCroppedNIIShape = numpy.int16(numpy.array(FID["axialCroppedNIIShape"]))
	midSlice = numpy.int16(numpy.array(FID["midSlice"]))
	IMGxx = numpy.array(FID["IMGxx"])
	IMGyy = numpy.array(FID["IMGyy"])
	IMGzz = numpy.array(FID["IMGzz"])
	finalRotY = numpy.array(FID["finalRotY"])
	finalRotZ = numpy.array(FID["finalRotZ"])
	finalTransX = numpy.array(FID["finalTransX"])

	#midSagAVW = midSagAVW[:, ::-1]
	#midSagAVW = numpy.rot90(midSagAVW, -1)
	FID.close()
	
	#pylab.subplot(1, 4, 1); CCSegUtils.showIMG(seg['initialSeg'])
	#pylab.subplot(1, 4, 2); CCSegUtils.showIMG(seg['finalSeg'])
	#pylab.subplot(1, 4, 3); CCSegUtils.showIMG(seg['finalSegNoArtefacts'])
	#pylab.subplot(1, 4, 4); CCSegUtils.showIMG(seg['IMG'])
	
	# this is going to be the segmentation in original resampled shape
	finalSegResampledAVWSpace = numpy.zeros(resampledAVWShape, dtype = numpy.uint8)
	
	# place the segmentation in the bounding box specified by the Lucas Kanade bounding box
	finalSegResampledAVWSpace[LKCropRows[0]:(LKCropRows[1] + 1), LKCropCols[0]:(LKCropCols[1] + 1)] = segIMG
	
	midSagAVWX, midSagAVWY = numpy.meshgrid(midSagAVWxx, midSagAVWyy)
	# resample the midsagittal slice in resampled space to the original space
	finalSegMidSagAVWSpace = CCSegUtils.interp2q(resamplexx, resampleyy, numpy.double(finalSegResampledAVWSpace), midSagAVWX, midSagAVWY, interpmethod = 'nearest', extrapval = 0)
	#print finalSegMidSagAVWSpace.shape	
	#print axialCroppedNIIShape
	#print axialNIIShape
	#print skullCrop
	#print axialCroppedNIIShape
	finalSegNativeAxialCroppedSpace = numpy.zeros(axialCroppedNIIShape, dtype = numpy.int8)
	#print axialNIIShape
	#print finalSegNativeAxialCroppedSpace.shape
	#print finalSegMidSagAVWSpace.shape
	#quit()
	finalSegNativeAxialCroppedSpace[midSlice] = numpy.rot90(finalSegMidSagAVWSpace, -1)
	#if dilateToParaSag:
	finalSegNativeAxialCroppedSpace[midSlice - 1] = numpy.rot90(finalSegMidSagAVWSpace, -1)
	finalSegNativeAxialCroppedSpace[midSlice + 1] = numpy.rot90(finalSegMidSagAVWSpace, -1)
	
	#NIISaving = nibabel.Nifti1Image(numpy.uint8(finalSegNativeAxialCroppedSpace), axialNIIAffine)
	#nibabel.save(NIISaving, outputBase + "_native_cropped_beforetransform.nii.gz")

	#NIISaving = nibabel.Nifti1Image(finalSegNativeAxialCroppedSpace, NIIAffine)
	
	#nibabel.save(NIISaving, outputBase + "_native_test_axial.nii.gz")
	#print skullCrop
	
	#print finalSegNativeAxialCroppedSpace.shape
	#print IMGxx.shape
	#print IMGyy.shape
	#print IMGzz.shape
	finalSegNativeAxialCroppedSpace = CCSegPipeMidSagSymmetric.transformIMG(finalSegNativeAxialCroppedSpace, IMGxx, IMGyy, IMGzz, -finalRotY, -finalRotZ, -finalTransX, interpmethod = 'nearest') 
	#NIISaving = nibabel.Nifti1Image(numpy.uint8(finalSegNativeAxialCroppedSpace), axialNIIAffine)
	#nibabel.save(NIISaving, outputBase + "_native_cropped_aftertransform.nii.gz")
	
	
	#finalSegNativeAxialCroppedSpace = numpy.rot90(finalSegNativeAxialCroppedSpace, 1)
	
	#print skullCrop
	#print finalSegNativeAxialCroppedSpace.shape
	#print axialNIIShape
	#print NIIShape
	# uncrop using skull crop
	#finalSegNativeAxialSpace = numpy.zeros((NIIShape[1], NIIShape[0], NIIShape[2]), dtype = numpy.int8)
	finalSegNativeAxialSpace = numpy.zeros(axialNIIShape, dtype = numpy.int8)
	finalSegNativeAxialSpace[skullCrop[0, 0]:(skullCrop[0, 1] + 1), skullCrop[1, 0]:(skullCrop[1, 1] + 1), skullCrop[2, 0]:(skullCrop[2, 1] + 1)] = numpy.array(finalSegNativeAxialCroppedSpace)
	
	#finalSegNativeAxialSpace = numpy.rot90(finalSegNativeAxialSpace, -1);
	# transform to original orientation

	MNITemplate = CCSegUtils.MNI152FLIRTTemplate(skullStripped = False)
	MNITemplateNII = nibabel.load(MNITemplate)
	MNIAXCodes = nibabel.aff2axcodes(MNITemplateNII.get_affine())
	MNIOrnt = nibabel.orientations.axcodes2ornt(MNIAXCodes)
	
	MNIToNIITransformOrnt = nibabel.orientations.ornt_transform(MNIOrnt, NIIOrnt)
	
	nativeNIIIMG = nibabel.orientations.apply_orientation(finalSegNativeAxialSpace, MNIToNIITransformOrnt)
	NIISaving = nibabel.Nifti1Image(nativeNIIIMG, NIIAffine)
	nibabel.save(NIISaving, outputBase + "_native.nii.gz")

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

