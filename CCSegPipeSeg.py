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
import Otsu

# My Lucas-Kanade tracker
import LKTracker

from matplotlib.font_manager import FontProperties

import errno

# implements imfill(BW, 'holes')
def bwFillHoles(BW):
	
	assert(isinstance(BW, numpy.ndarray) and BW.dtype == numpy.bool)
	
	mask = numpy.pad(numpy.logical_not(BW), 1, mode='constant', constant_values = 1)
	marker = numpy.pad(numpy.zeros(BW.shape, dtype=numpy.bool), 1, mode='constant', constant_values = 1)
	
	numDilations = 0
	while True:
		oldMarker = numpy.array(marker)
		marker = numpy.logical_and(mask, scipy.ndimage.morphology.binary_dilation(marker, structure = numpy.ones((3, 3), dtype=numpy.bool)))

		if numpy.array_equal(oldMarker, marker):
			break;
		numDilations = numDilations + 1
		del oldMarker	
	print "bwFillHoles finished with " + str(numDilations) + " dilations"
	marker = numpy.take(marker, numpy.arange(1, marker.shape[0] - 1), axis = 0)
	marker = numpy.take(marker, numpy.arange(1, marker.shape[1] - 1), axis = 1)
	
	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(marker);
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(BW);

	return numpy.logical_not(marker)
	#return numpy.logical_or(BW, numpy.logical_not(marker))


# bwselect using masks
# label the Mask
# extract all values in the labeled Mask under the marker
# retain only labels that overlap with the marker
def bwSelectWithMask(mask, marker):
	outIMG = numpy.zeros(marker.shape, dtype = numpy.bool)
	if not numpy.any(marker):
		return outIMG
	else:
		BWLabels, numLabels = scipy.ndimage.measurements.label(mask, structure = numpy.ones([3, 3]))
		labelsUnderMarker = BWLabels[marker]
		labelsUnderMarker = labelsUnderMarker[labelsUnderMarker.nonzero()]
		labelsUnderMarker = numpy.unique(labelsUnderMarker)
		#print "labelsUnderMarker: " + str(labelsUnderMarker)
		for z in range(numpy.size(labelsUnderMarker)):
			outIMG = numpy.logical_or(outIMG, BWLabels == labelsUnderMarker[z])
		return outIMG
	
def regionProps(labelImage, properties):
	
	numLabels = numpy.max(labelImage)
	regionProps = dict()
	
	pixelLists = list()

	if numLabels >= 1:
		for curLabel in range(numLabels):
			pixelLists.append(numpy.where(labelImage == (curLabel + 1)))

	for z in range(len(properties)):
		if properties[z].lower() == 'area':
			if numLabels >= 1:
				#regionProps['area'] = scipy.ndimage.measurements.labeled_comprehension(labelImage, labelImage, numpy.arange(1, numLabels + 1), numpy.size, numpy.uint32, 0)
				regionProps['area'] = numpy.zeros((numLabels))
				#print numpy.size(pixelLists[curLabel][0])
				#print regionProps['area'].shape
				for curLabel in range(numLabels):
					#print numpy.size(pixelLists[curLabel][0])
					regionProps['area'][curLabel] = numpy.size(pixelLists[curLabel][0])
			else:
				regionProps['area'] = None
		elif properties[z].lower() == 'pixellist':
			regionProps['pixelList'] = pixelLists[:]
			#if numLabels >= 1:
			#	for curLabel in range(numLabels):
			#		regionProps['pixelList'].append(numpy.where(labelImage == (curLabel + 1)))
		elif properties[z].lower() == 'mask':
			regionProps['mask'] = list()
			if numLabels >= 1:
				for curLabel in range(numLabels):
					T = numpy.zeros(labelImage.shape, dtype = numpy.bool)
					T[pixelLists[curLabel]] = 1
					regionProps['mask'].append(numpy.array(T))

	return regionProps

def dicesCoefficient(A, B):
	return 2.0 * numpy.sum(numpy.logical_and(A, B)) / (numpy.sum(A) + numpy.sum(B))

def bwAreaOpen(BW, areaThresh):
	outIMG = numpy.zeros(BW.shape, dtype=numpy.bool)
	BWLabels, numLabels = scipy.ndimage.measurements.label(BW, structure = numpy.ones([3, 3]))
	if numLabels >= 1:
		BWRegionProps = regionProps(BWLabels, ['Area'])
	
		I = numpy.where(BWRegionProps['area'] > areaThresh)[0]
		
		for z in range(numpy.size(I)):
			outIMG = numpy.logical_or(outIMG, BWLabels == I[z] + 1)

	return outIMG


def imMatchHist(inputIMG, targetIMG):
	
	# performs histogram matching so that the intensities of InputIMG are remapped to correspond to the histogram of TargetIMG
	#targetIMGBins, binSpacing = numpy.linspace(numpy.min(targetIMG), numpy.max(targetIMG), num = numBins - 2, retstep = True)

	numBins = 10000

	targetIMGHist, targetIMGBins = numpy.histogram(numpy.ravel(targetIMG), bins = numBins, density = False)
	
	targetIMGHist = numpy.double(targetIMGHist)
	# density doesnt normalize to 1, so do that
	targetIMGHist = targetIMGHist / numpy.sum(targetIMGHist)
	targetIMGHistCDF = numpy.cumsum(targetIMGHist)
	# make the bins centres, not edges
	targetIMGBins = targetIMGBins[:-1] + numpy.diff(targetIMGBins) / 2
	#print targetIMGBins

	# remove zero count bins
	I = numpy.nonzero(targetIMGHist)
	targetIMGBins = targetIMGBins[I]
	targetIMGHistCDF = targetIMGHistCDF[I]
	
	inputIMGHist, inputIMGBins = numpy.histogram(numpy.ravel(inputIMG), bins = numBins, density = False)
	
	inputIMGHist = numpy.double(inputIMGHist)
	inputIMGHist = inputIMGHist / numpy.sum(inputIMGHist)
	inputIMGHistCDF = numpy.cumsum(inputIMGHist)
	inputIMGBins = inputIMGBins[:-1] + numpy.diff(inputIMGBins) / 2
	# remove zero-count bins
	I = numpy.nonzero(inputIMGHist)
	inputIMGBins = inputIMGBins[I]
	inputIMGHistCDF = inputIMGHistCDF[I]
	
	#pylab.subplot(1, 2, 1)
	#pylab.plot(inputIMGBins, inputIMGHistCDF, 'b-', targetIMGBins, targetIMGHistCDF, 'g-')
	#pylab.show()

	# find the point on the CDF of the input image of each of the input pixels
	CDFOfInputIMG = numpy.interp(numpy.ravel(inputIMG), inputIMGBins, inputIMGHistCDF)
	# find the bins that are the same point on the CDF of the target image
	matchedIMG = numpy.interp(CDFOfInputIMG, targetIMGHistCDF, targetIMGBins)

	matchedIMG = numpy.reshape(matchedIMG, inputIMG.shape)
	return matchedIMG
	
	#pylab.plot(targetIMGBins[:-1] + numpy.diff(targetIMGBins) / 2, targetIMGHist)
	#pylab.show()

	#print targetIMGHist.shape
	#print targetIMGBins.shape

	#% get the CDFs of the histogram of the target image
#	[TargetIMGHist, TargetIMGBins] = hist(TargetIMG(:), TargetIMGBins);
#	TargetIMGHist = TargetIMGHist ./ sum(TargetIMGHist);
#	TargetIMGHistCDF = cumsum(TargetIMGHist);
#	% get the CDFs of the histogram of the matched filter
#
#	InputIMGBins = linspace(min(InputIMG(:)), max(InputIMG(:)), NumBins - 2);
#	S = InputIMGBins(2) - InputIMGBins(1);
#	InputIMGBins = [InputIMGBins(1) - S, InputIMGBins, InputIMGBins(end) + S];
#	clear S;
#	[InputIMGHist, InputIMGBins] = hist(InputIMG(:), InputIMGBins);
#	InputIMGHist = InputIMGHist ./ sum(InputIMGHist);
#	InputIMGHistCDF = cumsum(InputIMGHist);
#
#	% use unique to create a CDF with no repeated elements
#	[UInputIMGHistCDF, I] = unique(InputIMGHistCDF);
#	UInputIMGHistCDFBins = InputIMGBins(I);
#
#	[UTargetIMGHistCDF, I] = unique(TargetIMGHistCDF);
#	UTargetIMGHistCDFBins = TargetIMGBins(I);
#	clear I J;
#
#	CDFInputIMGPixels = interp1(UInputIMGHistCDFBins, UInputIMGHistCDF, InputIMG(:));
#	MatchedIMG = interp1(UTargetIMGHistCDF, UTargetIMGHistCDFBins, CDFInputIMGPixels);
#
#	MatchedIMG = reshape(MatchedIMG, size(InputIMG));

#function [InitialResampledAVW, ...
#	MatchedFilterRemapped, ...
#	LKParameters, LKCost, ...
#	TX, TY, InterpX, InterpY, ...
#	TemplateLKIMG, TemplateProbLKIMG, FornixProbLKIMG, ...
#	ResampledGroundCropped, ResampledAVWCropped, ...
#	TemplateProbLKIMGCropped, FornixProbLKIMGCropped, ...
#	OriginalOtsuMask, OtsuMask, ...
#	TemplateOverlap, OtsuSEG] = do_lk_and_otsu(ResampledAVW, MatchedFilter, MatchedFilterProb, MatchedFilterFornixProb, ResampledGroundAVW, real_total_offset, DoLK)
#pylab.show()

	#pylab.imshow(RW, origin = 'lower')
	#print numpy.array([ normXCorrCenterI, normXCorrCenterJ])
	#pylab.imshow(resampledAVW, origin = 'lower')
	#pylab.plot(normXCorrCenterJ, normXCorrCenterI, marker='o', markerfacecolor='blue', markersize=12)
	#pylab.plot(I[1], I[0], marker='o', markerfacecolor='red', markersize=3, linestyle='none')
	#pylab.plot([centreXCorrMaxJ, centreXCorrMaxJ + croppedTemplateAVW.shape[1]], [centreXCorrMaxI, centreXCorrMaxI])
	#pylab.plot([centreXCorrMaxJ, centreXCorrMaxJ + croppedTemplateAVW.shape[1]], [centreXCorrMaxI + croppedTemplateAVW.shape[0], centreXCorrMaxI + croppedTemplateAVW.shape[0]])
	#pylab.plot([centreXCorrMaxJ, centreXCorrMaxJ], [centreXCorrMaxI, centreXCorrMaxI + croppedTemplateAVW.shape[0]])
	#pylab.plot([centreXCorrMaxJ + croppedTemplateAVW.shape[1], centreXCorrMaxJ + croppedTemplateAVW.shape[1]], [centreXCorrMaxI, centreXCorrMaxI + croppedTemplateAVW.shape[0]])

def segCCLKandOtsu(IMG, templateIMG, templateCCProbIMG, templateFornixProbIMG, groundIMG, initialTranslation, DoLK = True):
	#print "Initial translation: " + str(initialTranslation)
	#print "IMG size: " + str(IMG.shape)
	#print "Template size: " + str(templateIMG.shape)
	#print numpy.arange(initialTranslation[1], initialTranslation[1] + templateIMG.shape[1])
	croppedIMG = numpy.take(IMG, numpy.arange(initialTranslation[0], initialTranslation[0] + templateIMG.shape[1]), axis = 1)
	#print croppedIMG.shape
	croppedIMG = numpy.take(croppedIMG, numpy.arange(initialTranslation[1], initialTranslation[1] + templateIMG.shape[0]), axis = 0)

	remappedTemplateIMG = imMatchHist(templateIMG, croppedIMG)

	#remappedTemplateIMGHist, remappedTemplateIMGBins = numpy.histogram(numpy.ravel(remappedTemplateIMG), bins = 100, density = True)
	#remappedTemplateIMGBins = remappedTemplateIMGBins[:-1] + numpy.diff(remappedTemplateIMGBins)
	
	#croppedIMGHist, croppedIMGBins = numpy.histogram(numpy.ravel(croppedIMG), bins = 100, density = True)
	#croppedIMGBins = croppedIMGBins[:-1] + numpy.diff(croppedIMGBins)
	
	#pylab.plot(croppedIMGBins, croppedIMGHist)
	#pylab.plot(remappedTemplateIMGBins, remappedTemplateIMGHist)
	#pylab.show()

	if DoLK == True:
		numLKIterations = 100
	else:
		numLKIterations = 0

	LKParameters, LKCost = LKTracker.weightedAffineInvComp(IMG, remappedTemplateIMG, templateIMG / numpy.max(templateIMG), numpy.array([0, 0, 0, 0, initialTranslation[0], initialTranslation[1]]), numLKIterations)

#	[TX, TY, InterpX, InterpY] = coords_template_lk_img(LKParameters, ResampledAVW, MatchedFilter);
	
	# find the coordinates of the template after the warp onto the target image

	TX, TY, interpX, interpY = LKTracker.coordsOfAffineWarpedTemplate(LKParameters, IMG, templateIMG)
	#pylab.subplot(2, 2, 1);	CCSegUtils.showIMG(TX);
	#pylab.subplot(2, 2, 2);	CCSegUtils.showIMG(TY);
	#pylab.subplot(2, 2, 3);	CCSegUtils.showIMG(InterpX);
	#pylab.subplot(2, 2, 4);	CCSegUtils.showIMG(InterpY);
	#pylab.show()
	#quit()

	# resample the template in the image space
	templateLKIMG = CCSegUtils.interp2q(numpy.arange(1, templateIMG.shape[1] + 1), numpy.arange(1, templateIMG.shape[0] + 1), templateIMG, interpX, interpY, extrapval = numpy.nan)
	templateCCProbLKIMG = CCSegUtils.interp2q(numpy.arange(1, templateIMG.shape[1] + 1), numpy.arange(1, templateIMG.shape[0] + 1), templateCCProbIMG, interpX, interpY, extrapval = 0)
	templateFornixProbLKIMG = CCSegUtils.interp2q(numpy.arange(1, templateIMG.shape[1] + 1), numpy.arange(1, templateIMG.shape[0] + 1), templateFornixProbIMG, interpX, interpY, extrapval = 0)
	#pylab.subplot(2, 2, 1); CCSegUtils.showIMG(templateLKIMG);
	#pylab.subplot(2, 2, 2); CCSegUtils.showIMG(templateCCProbLKIMG);
	#pylab.subplot(2, 2, 3); CCSegUtils.showIMG(templateFornixProbLKIMG);
	#pylab.subplot(2, 2, 4); CCSegUtils.showIMG(numpy.isnan(templateLKIMG));
	#pylab.show()
	#quit()
	
	# find the valid bounding box of the template in the image space
	I = numpy.where(numpy.logical_not(numpy.isnan(templateLKIMG)))

	cropRows = numpy.arange(numpy.min(I[0]), numpy.max(I[0]) + 1)
	cropCols = numpy.arange(numpy.min(I[1]), numpy.max(I[1]) + 1)
	
	del I
	# crop the warped images
	croppedIMG = IMG.take(cropRows, axis = 0).take(cropCols, axis = 1)
	croppedTemplateLKIMG = templateLKIMG.take(cropRows, axis = 0).take(cropCols, axis = 1)
	croppedTemplateCCProbLKIMG = templateCCProbLKIMG.take(cropRows, axis = 0).take(cropCols, axis = 1)
	croppedTemplateFornixProbLKIMG = templateFornixProbLKIMG.take(cropRows, axis = 0).take(cropCols, axis = 1)
	
	if not groundIMG == None:
		croppedGroundIMG = groundIMG.take(cropRows, axis = 0).take(cropCols, axis = 1)
	else:
		croppedGroundIMG = None
	# reset the NaN values in the warped template to zero

	croppedTemplateLKIMG[numpy.where(numpy.isnan(croppedTemplateLKIMG))] = 0

#	[THRESH, OtsuSEG] = robust_otsu2(ResampledAVWCropped, [0.05, 0.98]);
#	OtsuSEG(ResampledAVWCropped > THRESH(end)) = 3;
#	OtsuMask = (OtsuSEG == 3);
	otsuSeg = Otsu.robustOtsu(croppedIMG, [0.05, 0.98], NumberClasses=3, maskOutZeros = True)
	otsuSegWM = (otsuSeg == 3)
	origOtsuSegWM = numpy.array(otsuSegWM)
	
	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(origOtsuSegWM);
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(croppedIMG);
	
	#pylab.show()
	#quit()
	
	del otsuSeg

	templateFornixMask = (croppedTemplateFornixProbLKIMG > 0.01)
#	FornixMaskIdx = find((FornixProbLKIMGCropped * 50) > 0.5);

#	OtsuMaskCC = bwconncomp(OtsuMask);
#	OtsuMaskL = labelmatrix(OtsuMaskCC);
#	OtsuMaskR = regionprops(OtsuMaskCC, 'Area', 'PixelIdxList');
	
	otsuSegWMLabels, numLabels = scipy.ndimage.measurements.label(otsuSegWM, structure = numpy.ones([3, 3]))
	otsuSegWMRegionProps = regionProps(otsuSegWMLabels, ['Area', 'Mask'])
	

	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(otsuSegWMLabels);
	
	#pylab.show()
	#quit()
	

	#	TemplateOverlap = zeros(length(OtsuMaskR), 1);

#	% keep all regions that have high overlap with template or big area
#	for z = 1:length(OtsuMaskR)
#		IDX = setdiff(OtsuMaskR(z).PixelIdxList, FornixMaskIdx);
#		TemplateOverlap(z) = sum(TemplateProbLKIMGCropped(IDX)) ./ length(IDX);
#		clear IDX;
#	end
	templateOverlap = numpy.zeros((numLabels))
	
	J = numpy.argmax(otsuSegWMRegionProps['area'])
	
	for z in range(numLabels):
		M = numpy.logical_and(otsuSegWMRegionProps['mask'][z], numpy.logical_not(templateFornixMask))
		
#		if z == J:
#			#print J
#			#print otsuSegWMRegionProps['area']
#			SR = 1
#			SC = 3
#			pylab.subplot(SR, SC, 1)
#			CCSegUtils.showIMG(otsuSegWMRegionProps['mask'][z])
#			pylab.subplot(SR, SC, 2)
#			CCSegUtils.showIMG(templateFornixMask)
#			pylab.subplot(SR, SC, 3)
#			CCSegUtils.showIMG(croppedTemplateCCProbLKIMG)
#
#			pylab.gcf().set_size_inches((20, 10), forward = True)
#			pylab.show()
#			quit()
#
		if numpy.any(M):
			#print numpy.sum(M)
			templateOverlap[z] = numpy.sum(croppedTemplateCCProbLKIMG[M]) / numpy.sum(M)
		del M
	#rint templateOverlap
	largeTemplateOverlapIDX = numpy.where(templateOverlap > 0.01)[0]
	
	#print otsuSegWMRegionProps['area'][largeTemplateOverlapIDX]
	#print templateOverlap[largeTemplateOverlapIDX]
	#quit()

	#print templateOverlap
	#print largeTemplateOverlapIDX
	#print otsuSegWMRegionProps['area']
	#print otsuSegWMRegionProps['area'][largeTemplateOverlapIDX]
	otsuSegCC = numpy.zeros(otsuSegWM.shape, dtype=numpy.bool)
	if numpy.size(largeTemplateOverlapIDX) > 0:

		I = numpy.where(otsuSegWMRegionProps['area'][largeTemplateOverlapIDX] >= 200)[0]
		
		if numpy.size(I) == 0:
			# select the region with the highest overlap

			I = numpy.argmax(templateOverlap)
			
			otsuSegCC = numpy.array(otsuSegWMRegionProps['mask'][I])
			#CCSegUtils.showIMG(otsuSegCC)

#			SR = 2
#			SC = 2
#			pylab.subplot(SR, SC, 1)
#			CCSegUtils.showIMG(otsuSegWMLabels == largeTemplateOverlapIDX[J] + 1)
#			
#			pylab.subplot(SR, SC, 2)
#			CCSegUtils.showIMG(croppedTemplateCCProbLKIMG)
#			
#			for z in range(numpy.size(largeTemplateOverlapIDX)):
#				otsuSegCC = numpy.logical_or(otsuSegCC, otsuSegWMLabels == largeTemplateOverlapIDX[z] + 1)
#			
#			pylab.subplot(SR, SC, 3)
#			CCSegUtils.showIMG(otsuSegCC)
#			
#			pylab.subplot(SR, SC, 4)
#			
#			I = numpy.argmax(otsuSegWMRegionProps['area'])
#
#			CCSegUtils.showIMG(otsuSegWMLabels == I + 1)
#
			#pylab.gcf().set_size_inches((20, 10), forward = True)
			#pylab.show()
			#quit()
		else:
			#print I
			for z in range(numpy.size(I)):

				otsuSegCC = numpy.logical_or(otsuSegCC, otsuSegWMRegionProps['mask'][largeTemplateOverlapIDX[I[z]]])#otsuSegWMLabels == largeTemplateOverlapIDX[I[z]] + 1)
	

	return (TX, TY, interpX, interpY, croppedIMG, croppedTemplateLKIMG, croppedTemplateCCProbLKIMG,	croppedTemplateFornixProbLKIMG, croppedGroundIMG, otsuSegCC, cropRows, cropCols, LKCost)

	#pylab.subplot(2, 2, 1); CCSegUtils.showIMG(croppedTemplateLKIMG);
	#pylab.subplot(2, 2, 2); CCSegUtils.showIMG(croppedTemplateCCProbLKIMG);
	#pylab.subplot(2, 2, 3); CCSegUtils.showIMG(croppedTemplateFornixProbLKIMG);
	#pylab.subplot(2, 2, 4); CCSegUtils.showIMG(otsuSegCC);
	
	#pylab.show()
	#quit()
	
#	I = find(TemplateOverlap > 0.01);
#	if(~isempty(I))
#		OtsuMask = ismember(OtsuMaskL, I);
#		OtsuMask = bwareaopen(OtsuMask, 200);
#		clear Junk OtsuMaskL OtsuMaskR I;
#	else
#		disp('Current XCORR and LK failed');
#		OtsuMask = false(size(OtsuMaskL));
#	end

def nearestAnglesDistances(templateSegMask, estSegMask):
	
	if numpy.all(templateSegMask == False) or numpy.all(estSegMask == False):
		return (None, None)
	
	templateSegContours, hierarchy = cv2.findContours(numpy.uint8(templateSegMask), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

	#assert(len(templateSegContours) == 1),"Template CC seg has multiple contours, this should never happen"
	if(not (len(templateSegContours) == 1)):
		return (None, None)

	
	smoothingFilter = numpy.array([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])
	smoothingFilter = numpy.atleast_2d(smoothingFilter)
	tangentFilter = numpy.array([1, 0, -1])
	tangentFilter = numpy.atleast_2d(tangentFilter)
	
	templateSegContour = numpy.squeeze(templateSegContours[0]).T

	#print smoothingFilter.shape
	#print templateSegContour.shape
	smoothedTemplateSegContour = scipy.ndimage.filters.convolve(numpy.double(templateSegContour), smoothingFilter, mode = 'wrap')
	tangentTemplateSegContour = scipy.ndimage.filters.convolve(numpy.double(smoothedTemplateSegContour), tangentFilter, mode = 'wrap')

	# normalise the tangents to unit magnitude, numpy does the broadcasting for us so we just need to divide by the magnitudes
	#print tangentTemplateSegContour
	#print numpy.sqrt(numpy.sum(tangentTemplateSegContour * tangentTemplateSegContour, axis = 0))
	
	T = numpy.sqrt(numpy.sum(tangentTemplateSegContour * tangentTemplateSegContour, axis = 0))
	T[numpy.where(T == 0)] = 1

	tangentTemplateSegContour = tangentTemplateSegContour / T
	
	del T
	#print templateSegContour[0]
	#print smoothedTemplateSegContour[0]

	#print smoothedTemplateSegContour[0].shape
	#print numpy.concatenate([smoothedTemplateSegContour[0], numpy.array([smoothedTemplateSegContour[0, 0]])])
	#print smoothedTemplateSegContour.shape

	#CCSegUtils.showIMG(templateSegMask); CCSegUtils.plotContour(smoothedTemplateSegContour);
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()
	estSegContoursFirst, hierarchy = cv2.findContours(numpy.uint8(estSegMask), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
	del hierarchy
	
	estSegContours = list()

	for z in range(len(estSegContoursFirst)):
		estSegContours.append(numpy.atleast_2d(numpy.squeeze(estSegContoursFirst[z]).T))
	
	del estSegContoursFirst
	
	smoothedEstSegContours = list()
	tangentEstSegContours = list()
	
	#CCSegUtils.showIMG(allSegMasks)
	#CCSegUtils.plotContour(smoothedTemplateSegContour)
	for z in range(len(estSegContours)):
		if estSegContours[z].shape[0] == 2:
			smoothedEstSegContours.append(scipy.ndimage.filters.convolve(numpy.double(estSegContours[z]), smoothingFilter, mode = 'wrap'))
			tangentEstSegContours.append(scipy.ndimage.filters.convolve(numpy.double(smoothedEstSegContours[-1]), tangentFilter, mode = 'wrap'))
		#print smoothedEstSegContours[z].shape
		#CCSegUtils.plotContour(smoothedEstSegContours[z])
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()
	
	allSmoothedEstSegContours = numpy.concatenate(smoothedEstSegContours, axis = 1)
	allTangentEstSegContours = numpy.concatenate(tangentEstSegContours, axis = 1)
	
	#print smoothedEstSegContours[z]

	allTangentEstSegContoursMAG = numpy.sqrt(numpy.sum(allTangentEstSegContours * allTangentEstSegContours, axis = 0))

	I = numpy.where(allTangentEstSegContoursMAG == 0)[0]
	
	if numpy.size(I) > 0:
		
		IDXToDelete = numpy.concatenate((I, numpy.mod(I + 1, numpy.size(allTangentEstSegContoursMAG))))

		#CCSegUtils.showIMG(estSegMask)
		#CCSegUtils.plotContour(allSmoothedEstSegContours)
		#pylab.plot(allSmoothedEstSegContours[0, IDXToDelete], allSmoothedEstSegContours[1, IDXToDelete], 'g*')
		#pylab.gcf().set_size_inches((20, 10), forward = True)
		#pylab.show()
		#quit()

		allSmoothedEstSegContours = numpy.delete(allSmoothedEstSegContours, I, axis = 1)
		allTangentEstSegContours = scipy.ndimage.filters.convolve(allSmoothedEstSegContours, tangentFilter, mode = 'wrap')

		#allTangentEstSegContoursMAG = 

	#rint allSmoothedEstSegContours.shape
	del smoothedEstSegContours
	del tangentEstSegContours

	# find the indices of the closest points on the estimated contours for each of the template points
	# put the tempate coordinates along the X axis, estimated coordinates along the Y axis
	# so therefore to evaluate distances for each template coordinate use axis = 0, to process each column
	smoothedTemplateSegContourX, allSmoothedEstSegContoursX = numpy.meshgrid(numpy.atleast_2d(smoothedTemplateSegContour[0]), numpy.atleast_2d(allSmoothedEstSegContours[0]))
	smoothedTemplateSegContourY, allSmoothedEstSegContoursY = numpy.meshgrid(numpy.atleast_2d(smoothedTemplateSegContour[1]), numpy.atleast_2d(allSmoothedEstSegContours[1]))

	D = numpy.sqrt((smoothedTemplateSegContourX - allSmoothedEstSegContoursX) * (smoothedTemplateSegContourX - allSmoothedEstSegContoursX) + (smoothedTemplateSegContourY - allSmoothedEstSegContoursY) * (smoothedTemplateSegContourY - allSmoothedEstSegContoursY))
	
	closestIDX = numpy.argmin(D, axis = 0)
	closestDistances = D[(closestIDX, numpy.arange(D.shape[1]))]
	del D

	del smoothedTemplateSegContourX; del allSmoothedEstSegContoursX;
	del smoothedTemplateSegContourY; del allSmoothedEstSegContoursY;

	# normalise the tangents to unit magnitude, numpy does the broadcasting for us so we just need to divide by the magnitudes
	#print allTangentEstSegContours

	allTangentEstSegContours = allTangentEstSegContours / numpy.sqrt(numpy.sum(allTangentEstSegContours * allTangentEstSegContours, axis = 0))

	# compute the dot product between template tangent vectors and those in the closest indices in the EstSegContours
	dotProducts = numpy.minimum(numpy.abs(numpy.sum(numpy.take(allTangentEstSegContours, closestIDX, axis = 1) * tangentTemplateSegContour, axis = 0)), 1.0)
	# clip to [0, 1]
	#dotProducts = numpy.abs(dotProducts)
	#dotProducts = numpy.minimum(dotProducts, 1)
	return (closestDistances, dotProducts)
	#dotProduct = numpy.maximum(dotProduct, -1)
	
	#CCSegUtils.showIMG(D)
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#allSegMasks = numpy.zeros(estSegMask.shape)

	#allSegMasks[numpy.where(numpy.logical_and(numpy.logical_not(templateSegMask), estSegMask))] = 1
	#allSegMasks[numpy.where(numpy.logical_and(templateSegMask, numpy.logical_not(estSegMask)))] = 2
	#allSegMasks[numpy.where(numpy.logical_and(templateSegMask, estSegMask))] = 3


	#quit()

def bwJoinMSTBoundaries(BW):
	assert(isinstance(BW, numpy.ndarray) and BW.dtype == numpy.bool)
	
#	[B] = bwboundaries(BW, 8);
	contours, hierarchy = cv2.findContours(numpy.uint8(BW), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
	del hierarchy

	if len(contours) == 1:
		return numpy.zeros(BW.shape, dtype = numpy.bool)
	else:
		for curContour in range(len(contours)):
			contours[curContour] = numpy.array(numpy.squeeze(contours[curContour]))
		
		# contourDistances contains the minimum distances contourDistances[I, J] = minimum distance between contours A and B
		contourDistances = numpy.zeros((len(contours), len(contours)))
		# contourMinIDXI[I, J] and contourMinIDXJ[I, J] are the indices of the elements in contours I and J, respectively that form the distance contourDistances[I, J]
		contourMinIDXI = numpy.zeros((len(contours), len(contours)))
		contourMinIDXJ = numpy.zeros((len(contours), len(contours)))
			
		for contourI in range(len(contours) - 1):
			for contourJ in range(contourI + 1, len(contours)):

				XI, XJ = numpy.meshgrid(numpy.take(contours[contourI], [0], axis = 1), numpy.take(contours[contourJ], [0], axis = 1))
				YI, YJ = numpy.meshgrid(numpy.take(contours[contourI], [1], axis = 1), numpy.take(contours[contourJ], [1], axis = 1))
				D = numpy.double((XI - XJ) * (XI - XJ) + (YI - YJ) * (YI - YJ))
				I = numpy.argmin(numpy.sqrt(D))
				H = numpy.unravel_index(I, XI.shape)
				
				contourDistances[contourI, contourJ] = D[H[0], H[1]]
				contourMinIDXI[contourI, contourJ] = H[1]
				contourMinIDXJ[contourI, contourJ] = H[0]
				del H
				del XI
				del XJ
				del YI
				del YJ
				del I
				del D
			#	CCSegUtils.showIMG(BW)
			#	pylab.plot(numpy.take(contours[contourI], [0], axis = 1), numpy.take(contours[contourI], [1], axis = 1), 'r')
			#	pylab.plot(numpy.take(contours[contourJ], [0], axis = 1), numpy.take(contours[contourJ], [1], axis = 1), 'b')
			#	pylab.plot(contours[contourI][H[1], 0], contours[contourI][H[1], 1], 'm*') 
			#	pylab.plot(contours[contourJ][H[0], 0], contours[contourJ][H[0], 1], 'm*') 
			#	pylab.show()


				#print D[H[1], H[0]]
			#	quit()
			#for contourJ in range(contourI + 1, len(contours)):
		#for contourI in range(len(contours) - 1):
		#contourDistances = contourDistances + contourDistances.T
		#contourMinIDXI = numpy.uint32(contourMinIDXI + contourMinIDXI.T)
		#contourMinIDXJ = numpy.uint32(contourMinIDXJ + contourMinIDXJ.T)
		
		#CCSegUtils.showIMG(BW)
		
		#contourI = 0
		#contourJ = 1

		#pylab.plot(numpy.take(contours[contourI], [0], axis = 1), numpy.take(contours[contourI], [1], axis = 1), 'r')
		#pylab.plot(numpy.take(contours[contourJ], [0], axis = 1), numpy.take(contours[contourJ], [1], axis = 1), 'b')

	#	pylab.plot(AX, AY, 'm-')
	#	pylab.show()
		
		if len(contours) == 2:
			edges = numpy.atleast_2d(numpy.array([0, 1], dtype=numpy.uint32))
		else:
			# minimum spanning tree
			#print contourDistances
			C = scipy.sparse.csr_matrix(contourDistances, dtype = numpy.double)
			MSP = scipy.sparse.csgraph.minimum_spanning_tree(C)
			edges = MSP.nonzero()
			edges = numpy.concatenate((numpy.atleast_2d(edges[0]), numpy.atleast_2d(edges[1])), axis = 0).T
			#edges = numpy.atleast_2d(numpy.array([0, 1], dtype=numpy.uint32))
			#print edges
			del C; del MSP;
			#quit()
		
		joiningSegments = numpy.zeros(BW.shape, dtype=numpy.bool)
		for curEdge in range(edges.shape[0]):
			contourI = edges[curEdge, 0]
			contourJ = edges[curEdge, 1]

			AX = numpy.array([contours[contourI][contourMinIDXI[contourI, contourJ], 0], contours[contourJ][contourMinIDXJ[contourI, contourJ], 0]])
			AY = numpy.array([contours[contourI][contourMinIDXI[contourI, contourJ], 1], contours[contourJ][contourMinIDXJ[contourI, contourJ], 1]])
			
			arcLength = numpy.sqrt((AX[0] - AX[1]) * (AX[0] - AX[1]) + (AY[0] - AY[1]) * (AY[0] - AY[1]))
			
			n = numpy.ceil(arcLength * numpy.sqrt(2.0))
			IX = numpy.linspace(AX[0], AX[1], n)
			IY = numpy.linspace(AY[0], AY[1], n)

			IX = numpy.uint32(numpy.round(IX))
			IY = numpy.uint32(numpy.round(IY))
			T = numpy.zeros_like(joiningSegments)

			T[(IY, IX)] = True
			
			del IX;	del IY;
			del n
			del arcLength
			del contourI; del contourJ;

			Angle = numpy.arctan2(AY[1] - AY[0], AX[1] - AX[0]);
			Angle = numpy.mod(Angle + 2 * numpy.pi, 2 * numpy.pi);

			R = numpy.array([5, 5]) + 7 * numpy.abs(numpy.array([numpy.cos(Angle), numpy.sin(Angle)]))
			AngleWeighting = -0.9 * numpy.cos(2 * (Angle + 45 * numpy.pi / 180));
			SQ = numpy.sqrt(R[0] * R[1]) * AngleWeighting;
			SIGMA = numpy.array([[R[0], SQ], [SQ, R[1]]])
			
			F, FWHM = CCSegUtils.gaussianFWHM2D(SIGMA)

			joiningSegments = numpy.logical_or(joiningSegments, scipy.ndimage.morphology.binary_dilation(T, structure = (F > FWHM)))

			del F; del FWHM; del Angle; del R; del SQ; del SIGMA;
#			[F, HalfMaximum] = gaussian_fwhm2d(SIGMA);
			del AX;	del AY;
#
#			JoiningSegments = imdilate(T, double(F > HalfMaximum));
#
		#CCSegUtils.showIMG(joiningSegments)
		#pylab.show()
		return joiningSegments		

		#quit()


#	if(length(B) == 1)
#		JoiningSegments = false(size(BW));
#	else
#
#		Distances = zeros(length(B));
#		DistancesMinI = zeros(length(B));
#		DistancesMinJ = zeros(length(B));
#
#		for BoundaryI = 1:length(B) - 1
#			for BoundaryJ = BoundaryI + 1:length(B)
#				XC = bsxfun(@minus, B{BoundaryI}(:, 1), B{BoundaryJ}(:, 1)');
#				YC = bsxfun(@minus, B{BoundaryI}(:, 2), B{BoundaryJ}(:, 2)');
#				SZ = size(XC);
#				XC = XC(:);
#				YC = YC(:);
#				[Distances(BoundaryI, BoundaryJ), I] = min(sqrt(XC .* XC + YC .* YC));
#				[DistancesMinI(BoundaryI, BoundaryJ), DistancesMinJ(BoundaryI, BoundaryJ)] = ind2sub(SZ, I);
#			end
#		end
#		Distances = Distances + Distances';
#		DistancesMinI = DistancesMinI + DistancesMinI';
#		DistancesMinJ = DistancesMinJ + DistancesMinJ';
#		if(length(B) > 2)
#			edges = prim(Distances);
#			% this is because we calculate the distances based on increasing indices of edges
#			% if the edges are returned with a greater index on the left it
#			% will create errors downstream, so sort them
#			edges = sort(edges, 2); 
#		else
#			edges = [1, 2];
#		end
#
#		for z = 1:size(edges, 1)
#			StartX = B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 2);
#			EndX = B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 2);
#			StartY = B{edges(z, 1)}(DistancesMinI(edges(z, 1), edges(z, 2)), 1);
#			EndY = B{edges(z, 2)}(DistancesMinJ(edges(z, 1), edges(z, 2)), 1);
#
#			ArcLength = sqrt((StartX - EndX) * (StartX - EndX) + (StartY - EndY) * (StartY - EndY));
#			n = ceil(ArcLength * sqrt(2)); % make sure we cover all pixels
#
#			IX = linspace(StartX, EndX, n);
#			IY = linspace(StartY, EndY, n);
#			IX = round(IX);
#			IY = round(IY);
#			I = sub2ind(size(BW), IY, IX);
#			I = unique(I);
#			T = false(size(BW));
#			T(I) = 1;
#			Angle = atan2(EndY - StartY, EndX - StartX);
#			Angle = mod(Angle + 2 * pi, 2 * pi);
#
#			R = [5, 5] + 7 * abs([cos(Angle), sin(Angle)]);
#
#			AngleWeighting = -0.9 * cos(2 * (Angle + 45 * pi / 180));
#
#			SQ = sqrt(R(1) * R(2)) * AngleWeighting;
#			SIGMA = [R(1), SQ; SQ, R(2)];
#			[F, HalfMaximum] = gaussian_fwhm2d(SIGMA);
#
#			JoiningSegments = imdilate(T, double(F > HalfMaximum));
#
#			clear I IX IY n ArcLength EndY StartY EndX StartX F T HalfMaximum SQ Angle AngleWeighting R;
#		end
#	end
#

def radialStrel(R):
	X, Y = numpy.meshgrid(numpy.arange(-R, R + 1), numpy.arange(-R, R + 1))

	S = numpy.sqrt(X * X + Y * Y)

	return (S <= R)

def segCC(outputBase, groundTruthFile = None, doLK = True, doGraphics = False, segNoAuto = False):
	
	if doGraphics:
		(outputDir, subjectID) = os.path.split(outputBase)
		PNGDirectory = os.path.join(outputDir, "seg")
		try:
			os.makedirs(PNGDirectory)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(PNGDirectory):
				pass
			else:
				raise Exception


	# bwFillHoles test

	#A = numpy.zeros((100, 100), dtype=numpy.bool)
	
	#A[24:75, 24:75] = True
	#A[39:51, 39:51] = False
	
	#bwFillHoles(A)
	#CCSegUtils.showIMG(A);
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()
	# testing interp2

	#xx = numpy.arange(1, 5)
	#yy = numpy.arange(1, 6)
	#V = numpy.random.normal(0, 1, (yy.size, xx.size))

	#xi = numpy.array([0, 1, xx[-1]])
	#yi = numpy.array([0, 1, yy[-1]])

	#rint V.shape
	#rint xx.shape
	#rint yy.shape
	#F = interp2q(xx, yy, V, xi, yi)
	#print F
	#quit()

	midSagMATFile = outputBase + "_midsag.hdf5"
	
	assert(os.path.isfile(midSagMATFile)),"midsag hdf5 file not found"

	FID = h5py.File(midSagMATFile, 'r')
	
	NIIPixdims = numpy.array(FID['NIIPixdims'])
	midSagAVW = numpy.array(FID['midSagAVW'])
	#skullCrop = numpy.array(FID["skullCrop"]) # the initial cropping indices of the background
	#originalOrientationString = str(numpy.array(FID['originalOrientationString']))
	#originalNativeFile = str(numpy.array(FID['originalNativeFile']))
	#originalNativeCroppedFile = str(numpy.array(FID['originalNativeCroppedFile']))
	flirtTemplateFile = str(numpy.array(FID['flirtTemplateFile']))
	flirtMAT = numpy.array(FID["flirtMAT"]) # the transformation between  originalNativeCroppedFile -> flirtTemplateFile
	#flirtCropZerosRows = numpy.array(FID["flirtCropZerosRows"])
	#flirtCropZerosCols = numpy.array(FID["flirtCropZerosCols"])

	#MSPMethod = str(numpy.array(FID['MSPMethod']))
	
	#midSagAVW = midSagAVW[:, ::-1]
	#midSagAVW = numpy.rot90(midSagAVW, -1)
	FID.close()	
	del FID
	
	
#	FID = h5py.File(midSagMATFile, 'w')
#
#	FID.create_dataset("NIIPixdims", data=NIIPixdims, compression = 'gzip')
#	FID.create_dataset("midSagAVW", data=midSagAVW, compression = 'gzip')
#	FID.create_dataset("MSPMethod", data=MSPMethod)
#
#	FID.close()
#	return
	#rint "NIIPixdims = " + str(NIIPixdims)
	#CCSegUtils.showIMG(midSagAVW)

	#pylab.show()
	#quit()
	
	if not groundTruthFile == None:
		assert(os.path.isfile(groundTruthFile)),"Ground truth file was given, but it does not exist: " + groundTruthFile

		groundTruthNII = nibabel.load(groundTruthFile)
		groundTruthAVW = GroundTruthNII.get_data()
		groundTruthAVW = numpy.rot90(GroundTruthAVW, 1)
	#if not groundTruthFile == None:
	
	# load atlases
	
	scriptPath = os.path.realpath(__file__)
	(scriptDir, tail) = os.path.split(scriptPath)

	atlasDir = os.path.join(scriptDir, 'data')
	
	fornixProbNII = nibabel.load(os.path.join(atlasDir, 'all_fornix_prob.nii.gz'))
	fornixProbAVW = numpy.rot90(numpy.squeeze(fornixProbNII.get_data()), 1)
	CCProbNII = nibabel.load(os.path.join(atlasDir, 'all_cc_prob.nii.gz'))
	CCProbAVW = numpy.rot90(numpy.squeeze(CCProbNII.get_data()), 1)
	templateNII = nibabel.load(os.path.join(atlasDir, 'all_msp_mean.nii.gz'))
	templateAVW = numpy.rot90(numpy.squeeze(templateNII.get_data()), 1)
	
	del atlasDir
	templatePixdims = templateNII.get_header().get_zooms()[1:]
	#print "midSagAVW.shape = " + str(midSagAVW.shape)

	# resample to the space of the template
	midSagAVWxx = numpy.arange(0, midSagAVW.shape[1]) * NIIPixdims[0]
	midSagAVWxx = midSagAVWxx - numpy.mean(midSagAVWxx)
	
	midSagAVWyy = numpy.arange(0, midSagAVW.shape[0]) * NIIPixdims[1]
	midSagAVWyy = midSagAVWyy - numpy.mean(midSagAVWyy)

	midSagAVWX, midSagAVWY = numpy.meshgrid(midSagAVWxx, midSagAVWyy)

	templatexx = numpy.arange(0, templateAVW.shape[1]) * templatePixdims[0]
	templatexx = templatexx - numpy.mean(templatexx) 

	templateyy = numpy.arange(0, templateAVW.shape[0]) * templatePixdims[1]
	templateyy = templateyy - numpy.mean(templateyy)
	
	# compute overlap box of midSagAVWxx, midSagAVWyy, templatexx, templateyy
	
	lowX  = numpy.maximum(midSagAVWxx[ 0], templatexx[ 0])
	lowY  = numpy.maximum(midSagAVWyy[ 0], templateyy[ 0])
	highX = numpy.minimum(midSagAVWxx[-1], templatexx[-1])
	highY = numpy.minimum(midSagAVWyy[-1], templateyy[-1])
	
	resamplexx = numpy.arange(lowX, highX + templatePixdims[0], templatePixdims[0])
	resampleyy = numpy.arange(lowY, highY + templatePixdims[1], templatePixdims[1])
	#$del lowX; del lowY; del highX; del highY;
	#resampleX, resampleY = numpy.meshgrid(numpy.arange(lowX, highX + templatePixdims[0], templatePixdims[0]), numpy.arange(lowY, highY + templatePixdims[1], templatePixdims[1]))
	resampleX, resampleY = numpy.meshgrid(resamplexx, resampleyy)	
	# set up the interpolator as a function
	#midSagAVWInterpolator = scipy.interpolate.interp2d(midSagAVWxx, midSagAVWyy, midSagAVW, fill_value=0)
	#resampledAVW = midSagAVWInterpolator(resamplexx, resampleyy)
	
	resampledAVW = CCSegUtils.interp2q(midSagAVWxx, midSagAVWyy, midSagAVW, resampleX, resampleY, interpmethod = 'linear', extrapval = 0)
	#del midSagAVWInterpolator
	
	if not groundTruthFile == None:
		#groundTruthAVWInterpolator = scipy.interpolate.interp2d(midSagAVWxx, midSagAVWyy, groundTruthAVW, fill_value=0)
		#resampledGroundAVW = groundTruthAVWInterpolator(resamplexx, resampleyy) 
		resampledGroundAVW = CCSegUtils.interp2q(midSagAVWxx, midSagAVWyy, groundTruthAVW, resampleX, resampleY, interpmethod = 'linear', extrapval = 0)
	else:
		resampledGroundAVW = None

	del resampleX; del resampleY;
	#del resamplexx; del resampleyy;
	del lowX; del highX; del lowY; del highY;

	# set up imcrop
	spatial_rect = numpy.array([163.7248,  168.6860,  201.4264,  125.0233])
	pixelsPerVerticalUnit = 1
	pixelsPerHorizUnit = 1
	
	ymin = 1
	xmin = 1
	pixelHeight = spatial_rect[3] * pixelsPerVerticalUnit;
	pixelWidth = spatial_rect[2] * pixelsPerHorizUnit;
	r1 = (spatial_rect[1] - ymin) * pixelsPerVerticalUnit + 1;
	c1 = (spatial_rect[0] - xmin) * pixelsPerHorizUnit + 1;
	r2 = numpy.round(r1 + pixelHeight);
	c2 = numpy.round(c1 + pixelWidth);
	r1 = numpy.round(r1);
	c1 = numpy.round(c1);
	
	r1 = numpy.int32(r1) - 1
	r2 = numpy.int32(r2)
	c1 = numpy.int32(c1) - 1
	c2 = numpy.int32(c2)

	del pixelsPerVerticalUnit; del pixelsPerHorizUnit; del pixelHeight; del pixelWidth; del xmin; del ymin;
	
	# crop 2D images by using IMG[r1:r2, c1:c2]
	croppedTemplateAVW = templateAVW[r1:r2, c1:c2]
	croppedCCProbAVW = CCProbAVW[r1:r2, c1:c2]	
	croppedFornixProbAVW = fornixProbAVW[r1:r2, c1:c2]

	#rint str(r1) + "," + str(r2) + ":" + str(c1) + "," + str(c2)
	#print V
	
	# this is for cases where all automatic attempts fail and you just want to do a manual segmentation
	if segNoAuto == True:
		print "Not doing automatic segmentation"
		outputMAT = outputBase + "_seg.hdf5"

		FID = h5py.File(outputMAT, 'w')
		
		FID.create_dataset("initialSeg", data = numpy.zeros(midSagAVW.shape, dtype = numpy.uint8), compression = 'gzip')
		FID.create_dataset("finalSeg", data = numpy.zeros(midSagAVW.shape, dtype = numpy.uint8), compression = 'gzip')
		FID.create_dataset("finalSegNoArtefacts", data = numpy.zeros(midSagAVW.shape, dtype = numpy.uint8), compression = 'gzip')
		FID.create_dataset("IMG", data = numpy.array(midSagAVW), compression = 'gzip')
		FID.create_dataset("templatePixdims", data = templatePixdims)

		# for backprojection
		FID.create_dataset("resampledAVWShape", data = resampledAVW.shape)
		FID.create_dataset("LKCropRows", data = numpy.array([0, midSagAVW.shape[0] - 1]))
		FID.create_dataset("LKCropCols", data = numpy.array([0, midSagAVW.shape[1] - 1]))
		FID.create_dataset("midSagAVWxx", data = midSagAVWxx)
		FID.create_dataset("midSagAVWyy", data = midSagAVWyy)
		FID.create_dataset("resamplexx", data = resamplexx)
		FID.create_dataset("resampleyy", data = resampleyy)

		FID.close()
		return

	# prepare WMSeg for normxcorr, do 3-class otsu on cropped AVW
	S = Otsu.robustOtsu(resampledAVW, [0.05, 0.98], NumberClasses=3, maskOutZeros = True)
	WMSeg = (S == 3)
	del S;

	I = numpy.nonzero(WMSeg)
	minI = numpy.min(I[0])

	#rint numpy.round(minI + 75 / templatePixdims[1])
	#rint numpy.round(resampledAVW.shape[0] / 2)
	#rint numpy.round(resampledAVW.shape[1] / 2)

	# compute the index within the NormXCorr2 image that would put the template in the centre of the image
	#NormXCorrCenterI = numpy.round(numpy.round(resampledAVW.shape[0] / 2) + croppedTemplateAVW.shape[0] / 2)
	#NormXCorrCenterJ = numpy.round(numpy.round(resampledAVW.shape[1] / 2) + croppedTemplateAVW.shape[1] / 2)

	F = CCSegUtils.maxGaussian1D(2)
	#print numpy.atleast_2d(F).T.shape
	WMSegSmoothed = scipy.ndimage.convolve(numpy.double(WMSeg), numpy.atleast_2d(F), mode='nearest')
	WMSegSmoothed = scipy.ndimage.convolve(WMSegSmoothed, numpy.atleast_2d(F).T, mode='nearest')
	del F; del WMSeg;

	normXCorrAVW = cv2.matchTemplate(numpy.single(WMSegSmoothed), numpy.single(croppedCCProbAVW), cv2.TM_CCORR_NORMED)
	#print "WMSegSmoothed.shape = " + str(WMSegSmoothed.shape)
		
	#normXCorrCenterI = minI + 75 / templatePixdims[1]  - croppedTemplateAVW.shape[0] / 2 #$numpy.round(numpy.round(resampledAVW.shape[0] / 2) + croppedTemplateAVW.shape[0] / 2)
	normXCorrCentreI = (WMSegSmoothed.shape[0] - croppedTemplateAVW.shape[0]) / 2
	normXCorrCentreJ = (WMSegSmoothed.shape[1] - croppedTemplateAVW.shape[1]) / 2
	del I; del minI;

	#pylab.plot(F)
	#CCSegUtils.showIMG(WMSegSmoothed)
	#pylab.show()
	del WMSegSmoothed;
	
	# make coordinate matrices for normXCorrAVW so that the origin is the pixel that centres the filter over the image
	normXCorrAVWxx = (numpy.arange(normXCorrAVW.shape[1]) - normXCorrCentreJ) * templatePixdims[0] 
	normXCorrAVWyy = (numpy.arange(normXCorrAVW.shape[0]) - normXCorrCentreI) * templatePixdims[1]
	
	normXCorrAVWX, normXCorrAVWY = numpy.meshgrid(normXCorrAVWxx, normXCorrAVWyy)

	RFromCentre = numpy.sqrt(numpy.multiply(normXCorrAVWX, normXCorrAVWX) + numpy.multiply(normXCorrAVWY, normXCorrAVWY))

	RW = 0.5 * numpy.tanh(-(RFromCentre - 40.0) / 5.0) + 0.5
	normXCorrAVW = normXCorrAVW * RW
	#normXCorrAVW = normXCorrAVW * 0
	# approximation of imregionalmax, MAX = IMG == greyscale_dilation(IMG)
	# this should be fine given that we have hills and valleys rather than constant greyscale values
	normXCorrAVWMaxI = numpy.where(numpy.logical_and(scipy.ndimage.morphology.grey_dilation(normXCorrAVW, size=(3, 3)) == normXCorrAVW, normXCorrAVW > 0))
	#print normXCorrAVWMaxI
	#quit()
	del RW; del normXCorrAVWxx; del normXCorrAVWyy; del normXCorrAVWX; del normXCorrAVWY;

	if numpy.size(normXCorrAVWMaxI) == 0:
		print "Warning, no regional maxima of normxcorr found, using centering coordinates only"

		LKInitialOffsets = numpy.array([normXCorrCentreI, normXCorrCentreJ])
	else:
		minI = numpy.argmin(RFromCentre[normXCorrAVWMaxI])
		#centreXCorrMaxI = normXCorrAVWMaxI[0][minI]
		#centreXCorrMaxJ = normXCorrAVWMaxI[1][minI]

		LKInitialOffsets = numpy.array([normXCorrAVWMaxI[0][minI], normXCorrAVWMaxI[1][minI]])

		# sort the normxcorr values in DESCENDING order
		sortedNormXCorrAVWMaxI = numpy.argsort(normXCorrAVW[normXCorrAVWMaxI])
		sortedNormXCorrAVWMaxI = sortedNormXCorrAVWMaxI[::-1]
		sortedNormXCorrAVWMaxI = numpy.array([normXCorrAVWMaxI[0][sortedNormXCorrAVWMaxI], normXCorrAVWMaxI[1][sortedNormXCorrAVWMaxI]])
			
		centreRank = numpy.nonzero(numpy.logical_and(sortedNormXCorrAVWMaxI[0, :] == normXCorrAVWMaxI[0][minI], sortedNormXCorrAVWMaxI[1, :] == normXCorrAVWMaxI[1][minI]))[0][0]
		# remove the centre
		sortedNormXCorrAVWMaxI = numpy.delete(sortedNormXCorrAVWMaxI, centreRank, axis = 1)
		#print sortedNormXCorrAVWMaxI
		# dont need to check for invalid maximums because the version of normxcorr2 does not use padding
		numLKTries = 2
		#$print sortedNormXCorrAVWMaxI
		#print normXCorrAVWMaxI
		#CCSegUtils.showIMG(normXCorrAVW)
		LKInitialOffsets = numpy.concatenate((numpy.atleast_2d(LKInitialOffsets), numpy.take(sortedNormXCorrAVWMaxI, numpy.arange(0, numpy.minimum(numLKTries, sortedNormXCorrAVWMaxI.shape[1])), axis = 1).T), axis = 0)

		# swap so that the first column is X and the second column is Y
		LKInitialOffsets = numpy.fliplr(LKInitialOffsets)
		
		#print "matched template shape = " + str(croppedCCProbAVW.shape)
		#print "normXCorrAVW.shape = " + str(normXCorrAVW.shape)
		
		#print "LKInitialOffsets"
		#print LKInitialOffsets

		#CCSegUtils.showIMG(normXCorrAVW)
		#pylab.plot(LKInitialOffsets[:, 0], LKInitialOffsets[:, 1], 'b*')
		#pylab.show()
		del minI;
	
	segCCLKandOtsuOutput = dict()

	# initialise all the LK output variables to empty lists
	#segCCLKandOtsuOutput['TX'], segCCLKandOtsuOutput['TY'], segCCLKandOtsuOutput['interpX'], segCCLKandOtsuOutput['interpY'], segCCLKandOtsuOutput['croppedIMG'], segCCLKandOtsuOutput['croppedTemplateLKIMG'], segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'], segCCLKandOtsuOutput['croppedTemplateFornixProbLKIMG'], segCCLKandOtsuOutput['croppedGroundIMG'], segCCLKandOtsuOutput['otsuSegCC'] = (list(),) * 10
	segCCLKandOtsuOutput['TX'] = list()
	segCCLKandOtsuOutput['TY'] = list()
	segCCLKandOtsuOutput['interpX'] = list()
	segCCLKandOtsuOutput['interpY'] = list()
	segCCLKandOtsuOutput['croppedIMG'] = list()
	segCCLKandOtsuOutput['croppedTemplateLKIMG'] = list()
	segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'] = list()
	segCCLKandOtsuOutput['croppedTemplateFornixProbLKIMG'] = list()
	segCCLKandOtsuOutput['croppedGroundIMG'] = list()
	segCCLKandOtsuOutput['otsuSegCC'] = list()
	segCCLKandOtsuOutput['cropRows'] = list()
	segCCLKandOtsuOutput['cropCols'] = list()
	
	#print LKInitialOffsets.shape
	# used to perform the goodness of fit for the output segmentations from the LK method
	LKSegFits = dict()
	#LKSegFits['PCoverage'], LKSegFits['DiceCoefficients'], LKSegFits['Distances'], LKSegFits['Angles'], LKSegFits['SumDistances'], LKSegFits['SumAngles'], LKSegFits['TransformAngleX'], LKSegFits['TransformAngleY'], LKSegFits['LKCCSeg'] = (list(),) * 9
	LKSegFits['LKCost'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['PCoverage'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['DiceCoefficients'] = numpy.zeros((LKInitialOffsets.shape[0]))
	#LKSegFits['Distances'] = list()
	#LKSegFits['TangentDots'] = list()
	LKSegFits['SumDistances'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['SumTangentDots'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['TransformAngleX'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['TransformAngleY'] = numpy.zeros((LKInitialOffsets.shape[0]))
	LKSegFits['RigidFactor'] = numpy.zeros((LKInitialOffsets.shape[0]))
	
	LKCCSeg = list()

	for z in range(LKInitialOffsets.shape[0]):

		for curKey in segCCLKandOtsuOutput:
			segCCLKandOtsuOutput[curKey].append(None)
		for curKey in LKSegFits:
			if isinstance(LKSegFits[curKey], list):
				LKSegFits[curKey].append(None)

		LKCCSeg.append(None)

	for z in range(LKInitialOffsets.shape[0]):
		
		segCCLKandOtsuOutput['TX'][z], segCCLKandOtsuOutput['TY'][z], segCCLKandOtsuOutput['interpX'][z], segCCLKandOtsuOutput['interpY'][z], segCCLKandOtsuOutput['croppedIMG'][z], segCCLKandOtsuOutput['croppedTemplateLKIMG'][z], segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'][z], segCCLKandOtsuOutput['croppedTemplateFornixProbLKIMG'][z], segCCLKandOtsuOutput['croppedGroundIMG'][z], segCCLKandOtsuOutput['otsuSegCC'][z], segCCLKandOtsuOutput['cropRows'][z], segCCLKandOtsuOutput['cropCols'][z], LKSegFits['LKCost'][z] = segCCLKandOtsu(resampledAVW, croppedTemplateAVW, croppedCCProbAVW, croppedFornixProbAVW, resampledGroundAVW, LKInitialOffsets[z, :], DoLK = doLK)
		

		# compute the costs in LKSegFits
		#PCoverage(z) = sum(TemplateProbLKIMGCropped{z}(:) .* OtsuMask{z}(:)) ./ sum(TemplateProbLKIMGCropped{z}(:));
		LKSegFits['PCoverage'][z] = numpy.sum(segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'][z][numpy.nonzero(segCCLKandOtsuOutput['otsuSegCC'][z])]) / numpy.sum(segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'][z])
		
		# now do the transformation angles of the transformed coordinate matrices
		LKSegFits['TransformAngleX'][z] = numpy.arctan2(segCCLKandOtsuOutput['TY'][z][0, 1] - segCCLKandOtsuOutput['TY'][z][0, 0], segCCLKandOtsuOutput['TX'][z][0, 1] - segCCLKandOtsuOutput['TX'][z][0, 0])
		LKSegFits['TransformAngleY'][z] = numpy.arctan2(segCCLKandOtsuOutput['TY'][z][1, 0] - segCCLKandOtsuOutput['TY'][z][0, 0], segCCLKandOtsuOutput['TX'][z][1, 0] - segCCLKandOtsuOutput['TX'][z][0, 0]) - numpy.pi / 2.0
		# make a conservative mask of the template CC probability map
		templateBasedCCSeg = (segCCLKandOtsuOutput['croppedTemplateCCProbLKIMG'][z] > 0.3)
		
		# make a provisional CC segmentation
		LKCCSeg[z] = numpy.logical_and(segCCLKandOtsuOutput['otsuSegCC'][z], segCCLKandOtsuOutput['croppedTemplateFornixProbLKIMG'][z] <= 0.3)
		LKCCSeg[z] = bwAreaOpen(LKCCSeg[z], 50)
		if numpy.any(LKCCSeg[z]):
			# compute Dice's Coefficient between the 
			LKSegFits['DiceCoefficients'][z] = dicesCoefficient(templateBasedCCSeg, LKCCSeg[z])

		#	LKSegFits['Distances'][z], LKSegFits['TangentDots'][z] = nearestAnglesDistances(templateBasedCCSeg, LKCCSeg[z])
			#LKSegFits['SumDistances'][z] = numpy.sum(LKSegFits['Distances'][z])
			#LKSegFits['SumTangentDots'][z] = numpy.sum(LKSegFits['TangentDots'][z])
			Distances, TangentDots = nearestAnglesDistances(templateBasedCCSeg, LKCCSeg[z])
			
			if Distances != None:
				LKSegFits['SumDistances'][z] = numpy.sum(Distances)
				LKSegFits['SumTangentDots'][z] = numpy.sum(TangentDots)
			else:
				LKSegFits['SumDistances'][z] = numpy.inf
				LKSegFits['SumTangentDots'][z] = numpy.inf
			del Distances; del TangentDots; # we don't really need the individual distances/tangents, so discard them
		else:
			LKSegFits['DiceCoefficients'][z] = numpy.nan
			LKSegFits['SumDistances'][z] = numpy.inf
			LKSegFits['SumTangentDots'][z] = 0.0

	LKSegFits['RigidFactor'] = numpy.abs(LKSegFits['TransformAngleX'] - LKSegFits['TransformAngleY'])
	#print LKSegFits
	validTransform = numpy.logical_and(numpy.logical_and(numpy.logical_and(LKSegFits['TransformAngleX'] < numpy.pi / 3.0, LKSegFits['TransformAngleY'] < numpy.pi / 3.0), LKSegFits['RigidFactor'] < numpy.pi / 4.0), numpy.isfinite(LKSegFits['DiceCoefficients']))

	# testing
	#validTransform[1] = True

	validTransformIDX = numpy.where(validTransform)[0]
	#print validTransformIDX.shape

	if numpy.size(validTransformIDX) == 0:
		# problem, no valid registrations
		print "No valid registration results"
		return 1
	elif numpy.size(validTransformIDX) == 1:
		print "Only one valid registration"
		LKSegFitsScores = None
		bestLKIDX = validTransformIDX[0]
	else:
		# make a weight vector for each of the measures
		LKSegFitsWeights = dict()
		LKSegFitsWeights['LKCost'] = 3
		LKSegFitsWeights['PCoverage'] = 1
		LKSegFitsWeights['DiceCoefficients'] = 1
		LKSegFitsWeights['SumDistances'] = 1
		LKSegFitsWeights['SumTangentDots'] = 1
		LKSegFitsWeights['RigidFactor'] = 1
		
		# make a dict that tells me whether the lowest or highest value wins for each measure
		LKSegFitsWinDirection = dict()
		# if this is lowest then the index with the lowest value wins,
		LKSegFitsWinDirection['LKCost'] = 'lowest'
		LKSegFitsWinDirection['PCoverage'] = 'highest'
		LKSegFitsWinDirection['DiceCoefficients'] = 'highest'
		LKSegFitsWinDirection['SumDistances'] = 'lowest'
		LKSegFitsWinDirection['SumTangentDots'] = 'highest'
		LKSegFitsWinDirection['RigidFactor'] = 'lowest'
		
		LKSegFitsScores = numpy.zeros(numpy.size(validTransformIDX))
		
		#print LKSegFitsScores
		#print LKSegFitsScores.shape

		#print LKSegFitsWeights
		#print validTransformIDX[1]
		for curKey in LKSegFitsWeights:
			if LKSegFitsWinDirection[curKey] == 'lowest':
				I = numpy.argmin(LKSegFits[curKey][validTransformIDX])
			else:
				I = numpy.argmax(LKSegFits[curKey][validTransformIDX])
			
			#print I
			#print curKey
			#print LKSegFits[curKey][validTransformIDX]

			#print I
			#print validTransformIDX
			LKSegFitsScores[I] = LKSegFitsScores[I] + LKSegFitsWeights[curKey]
			#print LKSegFitsScores
			del I

		bestLKIDX = validTransformIDX[numpy.argmax(LKSegFitsScores)]
		del LKSegFitsWeights;
		del LKSegFitsWinDirection;
	
	#oLKGraphics = True
	if doGraphics:

		SR = 2
		SC = 3
		
		AX = CCSegUtils.empty2DList((SR, SC))
		for z in range(LKInitialOffsets.shape[0]):
			AX[0][z] = pylab.subplot(SR, SC, z + 1)
			#CCSegUtils.showIMG(segCCLKandOtsuOutput['croppedIMG'][z])
			CCSegUtils.showIMG(resampledAVW)

			pylab.plot([segCCLKandOtsuOutput['TX'][z][0, 0], segCCLKandOtsuOutput['TX'][z][-1, 0]], [segCCLKandOtsuOutput['TY'][z][0, 0], segCCLKandOtsuOutput['TY'][z][-1, 0]])
			pylab.plot([segCCLKandOtsuOutput['TX'][z][0, -1], segCCLKandOtsuOutput['TX'][z][-1, -1]], [segCCLKandOtsuOutput['TY'][z][0, -1], segCCLKandOtsuOutput['TY'][z][-1, -1]])
			pylab.plot([segCCLKandOtsuOutput['TX'][z][0, 0], segCCLKandOtsuOutput['TX'][z][0, -1]], [segCCLKandOtsuOutput['TY'][z][0, 0], segCCLKandOtsuOutput['TY'][z][0, -1]])
			pylab.plot([segCCLKandOtsuOutput['TX'][z][-1, 0], segCCLKandOtsuOutput['TX'][z][-1, -1]], [segCCLKandOtsuOutput['TY'][z][-1, 0], segCCLKandOtsuOutput['TY'][z][-1, -1]])
		
			AX[1][z] = pylab.subplot(SR, SC, z + 1 + SC)
			#CCSegUtils.showIMG(segCCLKandOtsuOutput['otsuSegCC'][z])
			CCSegUtils.showIMG(LKCCSeg[z])

		upNudge = 0.1	
		for curRow in range(SR):
			for curCol in range(SC):
				AXPos = numpy.array(AX[curRow][curCol].get_position().bounds)
				AXPos[1] = AXPos[1] + upNudge * (curRow + 1.0) / 2.0
				AX[curRow][curCol].set_position(AXPos)
				del AXPos
		for curCol in range(SC):
			AXText = list()
			for curKey in LKSegFits:
				AXText.append(curKey + ": " + "%.3f" % LKSegFits[curKey][curCol])

			FontProps = FontProperties()
			FontProps.set_size(10)
			if not curCol in validTransformIDX:
				AXText.append("Invalid")
				textColour = 'r'
			elif curCol == bestLKIDX:
				AXText.append("Best")
				FontProps.set_weight('bold')
				textColour = 'k'
			else:
				textColour = '0.25'

			pylab.text(0, -0.01, "\n".join(AXText), horizontalalignment='left', verticalalignment='top', transform=AX[1][curCol].transAxes, fontproperties = FontProps, color = textColour) # the transAxes is the same as 'Units', 'normalized'
			del FontProps;
			del textColour;
			del AXText;
		del upNudge
		
		pylab.gcf().set_size_inches((20, 10), forward = True)

		outputPNG = os.path.join(PNGDirectory, subjectID + "_lk.png")
		pylab.savefig(outputPNG)
		CCSegUtils.cropAutoWhitePNG(outputPNG)

		#pylab.show()
		del AX;
		del SR;
		del SC;
	
	# extract the best result
	
	#print bestLKIDX
	#print validTransformIDX
	#print len(segCCLKandOtsuOutput)
	#print len(LKCCSeg)
	CCSeg = numpy.array(LKCCSeg[bestLKIDX])
	
	LKCCSeg = numpy.array(CCSeg)
	#print CCSeg.shape	
	#print segCCLKandOtsuOutput.keys()
	#['interpY', 'interpX', 'croppedTemplateLKIMG', 'TX', 'TY', 'otsuSegCC', 'croppedTemplateFornixProbLKIMG', 'croppedGroundIMG', 'croppedTemplateCCProbLKIMG', 'croppedIMG']
	
	LKOutput = dict()

	for curKey in segCCLKandOtsuOutput:
		if segCCLKandOtsuOutput[curKey][bestLKIDX] != None:
			LKOutput[curKey] = numpy.array(segCCLKandOtsuOutput[curKey][bestLKIDX])
			#print curKey + str(LKOutput[curKey].shape)
		else:
			LKOutput[curKey] = None
	
	# delete working variables from the LK routine that are no longer needed
	del LKSegFitsScores; del segCCLKandOtsuOutput;

	# LKOutput contains: interpY interpX croppedTemplateLKIMG TX TY otsuSegCC croppedTemplateFornixProbLKIMG croppedGroundIMG croppedTemplateCCProbLKIMG croppedIMG 
	
	# mask CCSeg with a conservative threshold on the template
	CCSeg = numpy.logical_and(CCSeg, LKOutput['croppedTemplateCCProbLKIMG'] > 0.5)
	
	# make PenaltyImage
	
	# get the intensity properties under the current segmentation
	CCSegGaussianFit = dict()
	CCSegGaussianFit['mean'] = numpy.mean(LKOutput['croppedIMG'][CCSeg])
	CCSegGaussianFit['var'] = numpy.var(LKOutput['croppedIMG'][CCSeg])
	
	CCSegGaussianProb = CCSegUtils.normPDF(LKOutput['croppedIMG'], CCSegGaussianFit['mean'], numpy.sqrt(CCSegGaussianFit['var']))
	
	CCSegGaussianProb = CCSegGaussianProb / numpy.max(CCSegGaussianProb)
	CCSegGaussianProb[numpy.where(LKOutput['croppedIMG'] > CCSegGaussianFit['mean'])] = 1

	del CCSegGaussianFit

	# do edge strength on the gaussian prob image

	SIGMA = 0.25
	gaussianFilter = numpy.atleast_2d(CCSegUtils.maxGaussian1D(SIGMA))
	gaussianFilterDeriv = numpy.atleast_2d(CCSegUtils.maxGaussian1D(SIGMA, Derivative = 1))
	
	aSmooth = scipy.ndimage.convolve(CCSegGaussianProb, gaussianFilter, mode='nearest')
	aSmooth = scipy.ndimage.convolve(aSmooth, gaussianFilter.T, mode='nearest')
	
	ax = scipy.ndimage.convolve(aSmooth, gaussianFilterDeriv, mode='nearest')
	ay = scipy.ndimage.convolve(aSmooth, gaussianFilterDeriv.T, mode='nearest')
	
	#pylab.subplot(2, 2, 1); pylab.plot(gaussianFilter.squeeze())
	#pylab.subplot(2, 2, 2); pylab.plot(gaussianFilterDeriv.squeeze())
	#pylab.subplot(2, 2, 3); CCSegUtils.showIMG(ax)
	#pylab.subplot(2, 2, 4); CCSegUtils.showIMG(ay)

	CCSegGaussianProbGradMAG = numpy.sqrt(ax * ax + ay * ay)

	del aSmooth; del ax; del ay; del gaussianFilter; del gaussianFilterDeriv; del SIGMA;
	CCSegGaussianProbGradMAG = CCSegGaussianProbGradMAG / numpy.max(CCSegGaussianProbGradMAG)

	# make penalty images
	penaltyA = numpy.minimum(LKOutput['croppedTemplateFornixProbLKIMG'] * 50, 1) * CCSegGaussianProb
	penaltyB = numpy.array(CCSegGaussianProbGradMAG)
	penaltyC = numpy.exp(-2.0 * CCSegGaussianProb * (LKOutput['croppedTemplateFornixProbLKIMG'] + 2.0))
	
	# normalise penalty images
	penaltyA = penaltyA / numpy.max(penaltyA)
	penaltyB = penaltyB / numpy.max(penaltyB)
	penaltyC = penaltyC / numpy.max(penaltyC)
	
	penaltyImage = penaltyA + 4.0 * penaltyB + 2.0 * penaltyC

	if doGraphics:
		SR = 2
		SC = 2
		pylab.clf()
		pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(penaltyA); pylab.title('penaltyA')
		pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(penaltyB); pylab.title('penaltyB')
		pylab.subplot(SR, SC, 3); CCSegUtils.showIMG(penaltyC); pylab.title('penaltyC')
		pylab.subplot(SR, SC, 4); CCSegUtils.showIMG(penaltyImage); pylab.title('penaltyImage')

		pylab.gcf().set_size_inches((20, 10), forward = True)
		outputPNG = os.path.join(PNGDirectory, subjectID + "_penalty.png")
		pylab.savefig(outputPNG)
		CCSegUtils.cropAutoWhitePNG(outputPNG)

	
	CCSeg = bwSelectWithMask(penaltyImage < 0.8, CCSeg)
	#pylab.clf()
	#pylab.subplot(SR, SC, 1)
	#CCSegUtils.showIMG(CCSeg);
	#pylab.title('After select')
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#quit()

	CCSeg = bwAreaOpen(CCSeg, 50)
	#pylab.subplot(SR, SC, 2)

	#CCSegUtils.showIMG(CCSeg)
	#pylab.title('After bwareaopen')
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()

	#pylab.show()
	#quit()
	assert(numpy.any(CCSeg)),"CCSeg is empty after retaining 200 pixel regions"

	CCSeg = scipy.ndimage.morphology.binary_fill_holes(CCSeg)
	
	# join segments if there are multiple components
	CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(CCSeg, structure = numpy.ones([3, 3]))
	
	if numCCSegLabels > 1:
		CCJoiningSegments = bwJoinMSTBoundaries(CCSeg)
		CCSeg = numpy.logical_or(CCSeg, numpy.logical_and(CCJoiningSegments, CCSegGaussianProb > 0.3))
		del CCJoiningSegments
	
	CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(CCSeg, structure = numpy.ones([3, 3]))
	
	if numCCSegLabels > 1:
		CCJoiningSegments = bwJoinMSTBoundaries(CCSeg)
		CCSeg = numpy.logical_or(CCSeg, CCJoiningSegments)
		del CCJoiningSegments
	del CCSegLabels; del numCCSegLabels;	
	
	# close with radial STREL of radius 4
	CCSeg = scipy.ndimage.morphology.binary_closing(CCSeg, structure = radialStrel(4))
	
	# retain biggest connected component
	CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(CCSeg, structure = numpy.ones([3, 3]))
	
	CCSegRegionProps = regionProps(CCSegLabels, ['Area'])
	M = numpy.argmax(CCSegRegionProps['area'])
	CCSeg = (CCSegLabels == (M + 1))
	del M; del CCSegRegionProps; del CCSegLabels; del numCCSegLabels;
	
	# fill holes
	CCSeg = scipy.ndimage.morphology.binary_fill_holes(CCSeg)
	
	# reconstruct with the more liberal thresholding of penaltyImage
	BW = numpy.logical_or(penaltyImage < 1.0, CCSeg)
	BWS = bwSelectWithMask(BW, CCSeg)
	addedByReconstruct = numpy.logical_and(BWS, numpy.logical_not(CCSeg))
	CCSeg = numpy.logical_or(CCSeg, addedByReconstruct)
	del BW; del BWS;
	
	CCSeg = scipy.ndimage.morphology.binary_closing(CCSeg, structure = radialStrel(3))
	
	finalSeg = numpy.array(CCSeg)

	### ANTERIOR VESSEL REMOVAL ###
	SE = radialStrel(3)
	
	# compute the morphological gradient
	T = numpy.array(LKOutput['croppedIMG'])
	T[numpy.where(numpy.logical_not(CCSeg))] = -numpy.inf
	T = scipy.ndimage.morphology.grey_dilation(T, structure = SE)
	topHat = T - LKOutput['croppedIMG']

	T = numpy.array(LKOutput['croppedIMG'])
	T[numpy.where(numpy.logical_not(CCSeg))] = numpy.inf
	T = scipy.ndimage.morphology.grey_erosion(T, structure = SE)
	bottomHat = LKOutput['croppedIMG'] - T

	del SE;
	morphGrad = bottomHat + topHat
	morphGrad[numpy.where(numpy.logical_not(CCSeg))] = 0
	
	del T; del bottomHat; del topHat;

	morphGradMedian = numpy.median(morphGrad[numpy.where(CCSeg)])

	outliersStrong = morphGrad > (morphGradMedian * 3.0)
	outliersWeak = morphGrad > (morphGradMedian * 2.0)

	CCSegUtils.showIMG(outliersWeak)

	if numpy.any(outliersStrong):
		outliers = bwSelectWithMask(outliersWeak, outliersStrong)
	else:
		outliers = numpy.zeros_like(CCSeg)
		# test, create outliers everywhere
		#outliers = numpy.random.uniform(size = CCSeg.shape)
		#outliers = numpy.logical_and(outliers > 0.975, CCSeg)
		#outliers = scipy.ndimage.morphology.binary_dilation(outliers, structure = radialStrel(2))
		#outliers = numpy.logical_and(outliers, CCSeg)
	
	# cull outliers in the left 2/3 of the CCSeg

	# get the column indices of finalSeg
	I = numpy.nonzero(finalSeg)[1]

	leftColumn = numpy.uint32(numpy.min(I) + (numpy.max(I) - numpy.min(I)) * 2.0 / 3.0)
	outliers[:, 0:leftColumn] = False
	del I; del leftColumn;
	
	# produce finalSegNoVeins
	if not numpy.any(outliers):
		finalSegNoVeins = numpy.array(finalSeg)
	else:
		SE = radialStrel(3)
		outliersDilated = scipy.ndimage.morphology.binary_dilation(outliers, structure = SE)
		finalSegNoOutliers = numpy.logical_and(finalSeg, numpy.logical_not(outliersDilated))
		T = scipy.ndimage.morphology.binary_dilation(outliersDilated, structure = numpy.zeros((3, 3), dtype = numpy.bool))
		finalSegNoOutliersBorders = numpy.logical_and(T, finalSegNoOutliers)
		del T;

		finalSegOpened = numpy.logical_or(scipy.ndimage.morphology.binary_dilation(finalSegNoOutliersBorders, SE), finalSegNoOutliers)
		finalSegOpened = numpy.logical_and(finalSegOpened, finalSeg)

		CCSegLabels, numCCSegLabels = scipy.ndimage.measurements.label(finalSeg, structure = numpy.ones([3, 3]))
		if numCCSegLabels > 1:
			finalSegOpened = numpy.logical_or(finalSegOpened, bwJoinMSTBoundaries(finalSegOpened))
		del CCSegLabels; del numCCSegLabels;
		
		finalSegOpened = scipy.ndimage.morphology.binary_closing(finalSegOpened, structure = SE)
		finalSegNoVeins = scipy.ndimage.morphology.binary_fill_holes(finalSegOpened)
		del finalSegOpened; del finalSegNoOutliers; del finalSegNoOutliersBorders; del SE;
	
	# thicken thin regions
	SE = numpy.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype = numpy.bool)

	# firstly open to get rid of 
	finalSegNoVeinsOpened = scipy.ndimage.morphology.binary_opening(finalSegNoVeins, structure = SE)
	finalSegNoVeinsOpenedLost = numpy.logical_and(numpy.logical_not(finalSegNoVeinsOpened), finalSegNoVeins)

	finalSegNoVeinsOpenedLostLabels, finalSegNoVeinsOpenedLostNumLabels = scipy.ndimage.measurements.label(finalSegNoVeinsOpenedLost, structure = numpy.ones([3, 3]))
	
	finalSegNoVeinsOpenedThickened = numpy.array(finalSegNoVeinsOpened)
	# for each of the labels in the regions that were lost by opening (so thin bits)
	# find out if there was a disconnection, if so dilate the lost bit and add it to the mask
	for z in range(finalSegNoVeinsOpenedLostNumLabels):
		T = numpy.logical_and(finalSegNoVeinsOpened, numpy.logical_not(finalSegNoVeinsOpenedLostLabels == z + 1))
		TLabels, TNumLabels = scipy.ndimage.measurements.label(T, structure = numpy.ones([3, 3]))

		if TNumLabels > 1:
			finalSegNoVeinsOpenedThickened = numpy.logical_or(finalSegNoVeinsOpenedThickened, scipy.ndimage.morphology.binary_dilation(finalSegNoVeinsOpenedLostLabels == z + 1, structure = SE))
		del T; del TLabels; del TNumLabels;
	
	finalSegNoArtefacts = scipy.ndimage.morphology.binary_dilation(finalSegNoVeinsOpenedThickened, structure = SE)
	finalSegNoArtefacts = scipy.ndimage.morphology.binary_fill_holes(finalSegNoArtefacts)

	finalSeg = scipy.ndimage.morphology.binary_dilation(finalSeg, structure = SE)
	finalSeg = scipy.ndimage.morphology.binary_fill_holes(finalSeg, structure = SE)

	#T = scipy.ndimage.morphology.binary_opening(CCSeg, structure = R)
	#T[10:15, 50:55] = True

		
	#TJoiningSegments = bwJoinMSTBoundaries(T)
	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(finalSeg)
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(finalSegNoArtefacts)
	
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()

	if doGraphics:
		pylab.clf()
		pylab.subplot(1, 2, 1)
		CCSegUtils.showIMG(LKOutput['croppedIMG'])
		
		finalSegContours, hierarchy = cv2.findContours(numpy.uint8(finalSeg), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
		lineProps = {'color': 'b', 'linewidth': 2}
		for z in range(len(finalSegContours)):
			CCSegUtils.plotContour(numpy.squeeze(finalSegContours[z]).T, lineProps = lineProps)
		pylab.title('Final seg')
		
		pylab.subplot(1, 2, 2)
		CCSegUtils.showIMG(LKOutput['croppedIMG'])
		finalSegNoArtefactsContours, hierarchy = cv2.findContours(numpy.uint8(finalSegNoArtefacts), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
		lineProps = {'color': 'r', 'linewidth': 2}
		for z in range(len(finalSegNoArtefactsContours)):
			CCSegUtils.plotContour(numpy.squeeze(finalSegNoArtefactsContours[z]).T, lineProps = lineProps)
		
		pylab.title('Final seg + vein removal')
		pylab.gcf().set_size_inches((20, 10), forward = True)
		outputPNG = os.path.join(PNGDirectory, subjectID + "_seg.png")
		pylab.savefig(outputPNG)
		CCSegUtils.cropAutoWhitePNG(outputPNG)
	
	outputMAT = outputBase + "_seg.hdf5"

	FID = h5py.File(outputMAT, 'w')

#OutputS.InitialSeg = CCSeg;
#OutputS.FinalSeg = FinalSeg;
#OutputS.FinalSegNoArtefacts = FinalSegOpenedMThickened;
#OutputS.IMG = ResampledAVWCropped;
#OutputS.TemplatePixdims = TemplatePixdims;
	
	# transform finalSeg back to native space
	# the first is to reverse a cropping LKOutput['otsuSeg'] -> resampledAVW

	FID.create_dataset("initialSeg", data = numpy.uint8(LKCCSeg), compression = 'gzip')
	FID.create_dataset("finalSeg", data = numpy.uint8(finalSeg), compression = 'gzip')
	FID.create_dataset("finalSegNoArtefacts", data = numpy.uint8(finalSegNoArtefacts), compression = 'gzip')
	FID.create_dataset("IMG", data = LKOutput['croppedIMG'], compression = 'gzip')
	FID.create_dataset("templatePixdims", data = templatePixdims)

	# for backprojection
	FID.create_dataset("resampledAVWShape", data = resampledAVW.shape)
	FID.create_dataset("LKCropRows", data = numpy.array([numpy.min(LKOutput['cropRows']), numpy.max(LKOutput['cropRows'])]))
	FID.create_dataset("LKCropCols", data = numpy.array([numpy.min(LKOutput['cropCols']), numpy.max(LKOutput['cropCols'])]))
	FID.create_dataset("midSagAVWxx", data = midSagAVWxx)
	FID.create_dataset("midSagAVWyy", data = midSagAVWyy)
	FID.create_dataset("resamplexx", data = resamplexx)
	FID.create_dataset("resampleyy", data = resampleyy)

	FID.close()
#	pylab.clf()
#	
#	lineProps = {'color': 'r', 'linewidth': 2}
#	
#	SR = 2
#	SC = 2
#	pylab.subplot(SR, SC, 1)
#	CCSegUtils.showIMG(LKOutput['croppedIMG'])
#	T, hierarchy = cv2.findContours(numpy.uint8(finalSegNoArtefacts), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
#	for z in range(len(T)):
#		CCSegUtils.plotContour(numpy.squeeze(T[z]).T, lineProps = lineProps)
#
#	pylab.subplot(SR, SC, 2)
#	CCSegUtils.showIMG(resampledAVW)
#	T, hierarchy = cv2.findContours(numpy.uint8(finalSegResampledAVWSpace), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
#	for z in range(len(T)):
#		CCSegUtils.plotContour(numpy.squeeze(T[z]).T, lineProps = lineProps)
#
#	#resamplexx = numpy.arange(lowX, highX + templatePixdims[0], templatePixdims[0])
#	#resampleyy = numpy.arange(lowY, highY + templatePixdims[1], templatePixdims[1])
#	#resampleX, resampleY = numpy.meshgrid(resamplexx, resampleyy)	
#	pylab.subplot(SR, SC, 3)
#	CCSegUtils.showIMG(midSagAVW)
#	T, hierarchy = cv2.findContours(numpy.uint8(finalSegMidSagAVWSpace), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
#	for z in range(len(T)):
#		CCSegUtils.plotContour(numpy.squeeze(T[z]).T, lineProps = lineProps)
#	
#	#pylab.subplot(SR, SC, 4)
#	#CCSegUtils.showIMG(midSagAVW2)
#	#T, hierarchy = cv2.findContours(numpy.uint8(finalSegTemplateSpace), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
#	#for z in range(len(T)):
##		CCSegUtils.plotContour(numpy.squeeze(T[z]).T, lineProps = lineProps)
#	
#
#	pylab.gcf().set_size_inches((20, 10), forward = True)
#	pylab.show()
#	
#		SR = 3
#		SC = 3
#		
#		#SZ = size(MatchedFilter);
#
#		pylab.subplot(SR, SC, 1);
#		CCSegUtils.showIMG(normXCorrAVW);
#
#		pylab.plot(LKInitialOffsets[0, 0], LKInitialOffsets[0, 1], 'r*')
#		pylab.plot(LKInitialOffsets[1:, 0], LKInitialOffsets[1:, 1], 'b*')
#		pylab.title('LK initial offsets')
#		
#		pylab.subplot(SR, SC, 2);
#		CCSegUtils.showIMG(resampledAVW)
#		
#		lineProps = {'color': 'r', 'linewidth': 2}
#
#		plotRectangle([LKInitialOffsets[bestLKIDX, 0], LKInitialOffsets[bestLKIDX, 1], croppedTemplateAVW.shape[1], croppedTemplateAVW.shape[0]], lineProps)
#		
#		lineProps = {'color': 'b', 'linewidth': 2}
#		pylab.plot([LKOutput['TX'][ 0,  0], LKOutput['TX'][-1,  0]], [LKOutput['TY'][ 0,  0], LKOutput['TY'][-1,  0]], **lineProps)
#		pylab.plot([LKOutput['TX'][ 0, -1], LKOutput['TX'][-1, -1]], [LKOutput['TY'][ 0, -1], LKOutput['TY'][-1, -1]], **lineProps)
#		pylab.plot([LKOutput['TX'][ 0,  0], LKOutput['TX'][ 0, -1]], [LKOutput['TY'][ 0,  0], LKOutput['TY'][ 0, -1]], **lineProps)
#		pylab.plot([LKOutput['TX'][-1,  0], LKOutput['TX'][-1, -1]], [LKOutput['TY'][-1,  0], LKOutput['TY'][-1, -1]], **lineProps)
#		pylab.title('Initial XCORR and LK alignment');
#
#		pylab.subplot(SR, SC, 3)
#		CCSegUtils.showIMG(LKCCSeg)
#		pylab.title("Initial segmentation")
#
		#pylab.subplot(SR, SC, 4);
		#CCSegUtils.showIMG(LKOutput['croppedIMG']);
#		hold on;
#		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#		set(CC, 'Color', 'r');
#		
#		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#		set(CC, 'Color', 'g');
#		title('Original Image');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
#		
#		subplot(SR, SC, 5);
#		imshow(ResampledAVWCropped, []);
#		[I, J] = find(AddedByReconstruct);
#		if(~isempty(I))
#			hold on;
#			plot(J, I, '*');
#		end
#		%	[~, CC] = contour(AddedByReconstruct, [0.5, 0.5]);
#		%set(CC, 'Color', 'm');
#		
#		%[~, CC] = contour(OriginalPSI, [0, 0]);
#		%set(CC, 'Color', 'b');
#		%[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#		%set(CC, 'Color', 'g');
#	% 	if(DoReconstruct)
#	% 		title('Original Image, added by reconstruct');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
#	% 	else
#	% 		title('Original Image, no reconstruct');%, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']
#	% 	end
#
#		subplot(SR, SC, 6);
#		imshow(GaussianProb, []);
#		hold on;
#		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#		set(CC, 'Color', 'r');
#		%[~, CC] = contour(OriginalPSI, [0, 0]);
#		%set(CC, 'Color', 'b');
#		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#		set(CC, 'Color', 'g');
#		title('Gaussian probability of pixels');
#
#		subplot(SR, SC, 7);
#		imshow(TemplateProbLKIMGCropped, []);
#		hold on;
#		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#		set(CC, 'Color', 'r');
#	% 	[~, CC] = contour(OriginalPSI, [0, 0]);
#	% 	set(CC, 'Color', 'b');
#		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#		set(CC, 'Color', 'g');
#		title('Probability of CC from template');
#		text(0.5, -0.01, {'In all these images', 'red: result', 'green: artefacts removed'}, 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
#	% 	subplot(SR, SC, 7);
#	% 	imshow(TemplateProbLKIMGCropped, []);
#	% 	hold on;
#	% 	[~, CC] = contour(OriginalOtsuMask, [0.5, 0.5]);
#	% 	set(CC, 'Color', 'r');
#	% 	%title(['Probability of CC from template with original otsu borders, ' num2str(TemplateOverlapOtsu * 100, '%.2f') '% overlap']);
#	% 	title(['Probability of CC from template with original otsu borders']);
#
#	% 	subplot(SR, SC, 8);
#	% 	imshow(OriginalOtsuSeg, []);
#	% 	hold on;
#	% 	[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#	% 	set(CC, 'Color', 'r');
#	% 	[~, CC] = contour(OriginalPSI, [0, 0]);
#	% 	set(CC, 'Color', 'b');
#	% 	title(['OtsuMask ' num2str(OriginalOtsuSegCounts(1:3)', '%d ')]);
#
#		subplot(SR, SC, 8);
#		imshow(PenaltyImage, []);
#		hold on;
#		[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#		set(CC, 'Color', 'r');
#		%[~, CC] = contour(OriginalPSI, [0, 0]);
#		%set(CC, 'Color', 'b');
#		[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#		set(CC, 'Color', 'g');
#		title('Exclusion mask');
#
#		if(GroundTruthGiven)
#			subplot(SR, SC, 9);
#			imshow(ResampledAVWCropped, []);	
#			hold on;
#			[~, CC] = contour(ResampledGroundCropped, [0.5, 0.5]);
#			set(CC, 'Color', 'b');
#			[~, CC] = contour(FinalSeg, [0.5, 0.5]);
#			set(CC, 'Color', 'r');
#			[~, CC] = contour(FinalSegOpenedMThickened, [0.5, 0.5]);
#			set(CC, 'Color', 'g');
#			title({'Ground truth (blue)', 'Final (red)', 'Artefacts Removed (green)'});		
#
#		end
#		OutputPNG = fullfile(OutputDir, 'seg', [OutputPrefix '_seg.png']);
#		%keyboard;
#		[~, ~, ~] = mkdir(fullfile(OutputDir, 'seg'));
#		
#		FigPos = fullscreen_fig_pos;
#		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
#		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
#		
#		IMG = imread(OutputPNG);
#		IMG = imautocropwhite(IMG);
#		imwrite(IMG, OutputPNG);
#		%delete(gcf);
#	end

#	pylab.subplot(1, 3, 1); CCSegUtils.showIMG(penaltyImage < 0.8);
#	pylab.subplot(1, 3, 2); CCSegUtils.showIMG(CCSeg);
#	pylab.subplot(1, 3, 3); CCSegUtils.showIMG(T);
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()


