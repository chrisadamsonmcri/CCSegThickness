import shutil
import os
import tempfile
import errno

#import scipy.interpolate

import numpy
import nibabel
import pylab
import h5py

# for connected components functions
#import pymorph
import scipy.ndimage

import math

import cv2

import CCSegUtils

from matplotlib.font_manager import FontProperties

import LaplaceThicknessMethod

	#SigmaB = 
#	MuK = numpy.cumsum(T * hist)
#	MuT = numpy.sum(MuK)
#	
#	OmegaZero = OmegaK
#	OmegaOne = 1 - OmegaK
#	
#	MuZero = numpy.zeros(OmegaZero.shape)
#	MuOne = numpy.zeros(OmegaZero.shape)
#	
#	OmegaKGreaterZero = (OmegaK > 0)
#	OmegaKLessOne = (OmegaK < 1)
#	MuZero[OmegaKGreaterZero] = MuK[OmegaKGreaterZero] / OmegaK[OmegaKGreaterZero]
#	MuOne[OmegaKLessOne] = (MuT - MuK[OmegaKLessOne]) / (1 - OmegaK[OmegaKLessOne])
#	
#	#SigmaZero = numpy.zeros(OmegaZero.shape)
#	#SigmaOne = numpy.zeros(OmegaZero.shape)
#	
#	#for i in range(256):
#	#	if OmegaOne[i] > 0:
#	#		IndicesAfterI = numpy.arange(i + 1, 255, 1)
#	#		SigmaOne[i] = numpy.sum(hist[IndicesAfterI] * (IndicesAfterI - MuOne[i]) * (IndicesAfterI - MuOne[i]) / OmegaOne[i])
#	
#	SigmaB3 = OmegaZero * OmegaOne * (MuOne - MuZero) * (MuOne - MuZero)
#	
#	SigmaB2 = numpy.zeros(OmegaZero.shape)
#
#	SigmaB2Mask = numpy.logical_and(OmegaKGreaterZero, OmegaKLessOne)
#	#SigmaB2[numpy.where(SigmaB2Mask)] = ((MuT * OmegaK[numpy.where(SigmaB2Mask)] - MuK[numpy.where(SigmaB2Mask)]) * (MuT * OmegaK[numpy.where(SigmaB2Mask)] - MuK[numpy.where(SigmaB2Mask)])) / (OmegaK[numpy.where(SigmaB2Mask)] * (1 - OmegaK[numpy.where(SigmaB2Mask)]))
#	SigmaB2[SigmaB2Mask] = ((MuT * OmegaK[SigmaB2Mask] - MuK[SigmaB2Mask]) * (MuT * OmegaK[SigmaB2Mask] - MuK[SigmaB2Mask])) / (OmegaK[SigmaB2Mask] * (1 - OmegaK[SigmaB2Mask]))
#	#print SigmaB
#	#print SigmaB2
#
#	print numpy.argmax(SigmaB2) + 1
#	print numpy.argmax(SigmaB) + 1
#

#def otsu2(IMG):
#	
#	hist, bin_edges = numpy.histogram(IMG, bins=range(257), range=None, normed=False, weights=None, density=True)
#	CumSumHist = numpy.cumsum(hist)
#	
#	#hist = numpy.matrix(numpy.atleast_2d(hist)).T
#	#print hist.shape
#
#	#XMatrix = numpy.arange(256)
#
#	MaxSigmaB = 0
#	for I in range(256):
#		for J in range(I + 1, 256):
#			OmegaZero = CumSumHist[I]
#			OmegaOne = CumSumHist[J] - CumSumHist[I]
#			OmegaTwo = CumSumHist[-1] - CumSumHist[J]
#
#			if OmegaZero > 0 and OmegaOne > 0 and OmegaTwo > 0:
#
#				T = numpy.arange(I)
#				#T = XMatrix[0:I]
#				#print T
#				#T = numpy.arange(0, I)
#				#MuZero = (1 / OmegaZero) * numpy.sum(T * hist[0:I])
#				MuZero = (1 / OmegaZero) * numpy.sum(T * hist[0:I])
#				XC = T - MuZero
#				SigmaZero = (1 / OmegaZero) * numpy.sum(XC * XC * hist[0:I])
#
#				#SigmaZero = (1 / OmegaZero) * ( numpy.sum(XC * XC * hist[0:I])
#				
#				T = numpy.arange(I, J)
#				MuOne = (1 / OmegaOne) * numpy.sum(T * hist[I:J])
#				XC = T - MuOne
#				SigmaOne = (1 / OmegaOne) * numpy.sum(XC * XC * hist[I:J])
#				
#				T = numpy.arange(J, 256)
#				MuTwo = (1 / OmegaTwo) * numpy.sum(T * hist[J:])
#				XC = T - MuTwo
#				SigmaTwo = (1 / OmegaTwo) * numpy.sum(XC * XC * hist[J:])
#				
#				SigmaB = OmegaZero * OmegaOne * OmegaTwo * (MuZero - MuOne - MuTwo) * (MuZero - MuOne - MuTwo)
#				
#				if SigmaB > MaxSigmaB:
#					MaxSigmaB = SigmaB
#					T = numpy.array([I, J])
#	return T
#

#RegisteredProfiles = zeros(size(P));
#
#NumNodes = size(P, 1);
#
#%MeanProfile = nanmean(P, 2);
#
#SE = strel('arbitrary', ones(floor(NumNodes / 5), 1));
#
#%Peaks = find(MeanProfile == imdilate(MeanProfile, SE));
#
#%LeftPeak = Peaks(1);
#%RightPeak = Peaks(end);
#
#% hardcoded
#LeftPeak = round(NumNodes * 0.12);
#RightPeak = round(NumNodes * 0.90);
#MeanPeaks = [LeftPeak, RightPeak];
#
#
#%[~, MeanProfileFirstPeak] = max(MeanProfile(1:floor(NumNodes / 2)));
#%[~, MeanProfileSecondPeak] = max(MeanProfile(ceil(NumNodes / 2):end));
#
#%plot(MeanProfile);
#
#%hold on;
#%plot(MeanProfileFirstPeak, MeanProfile(MeanProfileFirstPeak), '*');
#%plot(MeanProfileSecondPeak, MeanProfile(MeanProfileSecondPeak), '*');
#%plot(Peaks, MeanProfile(Peaks), '*');
#
#AllPeaks = zeros(size(P, 2), 2);
#AllReparamX = zeros(size(P));
#for z = 1:size(P, 2)
#	T = P(:, z);
#	T(isnan(T)) = 0;
#	CurPeaks = find(T == imdilate(T, SE));
#	
#	
#	%plot(T);
#	%hold on;
#	%keyboard;
#	if(~isempty(CurPeaks))
#		CurLeftPeaks = CurPeaks(CurPeaks < NumNodes / 2);
#		CurRightPeaks = CurPeaks(CurPeaks > NumNodes / 2);
#
#		[~, I] = max(T(CurLeftPeaks));
#		CurLeftPeak = CurLeftPeaks(I);
#		[~, I] = max(T(CurRightPeaks));
#		CurRightPeak = CurRightPeaks(I);
#		%plot(CurPeaks, T(CurPeaks), '*');
#		%CurLeftPeak = CurPeaks(1);
#		%CurRightPeak = CurPeaks(end);
#		AllPeaks(z, 1) = CurLeftPeak;
#		AllPeaks(z, 2) = CurRightPeak;
#		%if(CurLeftPeak ~= LeftPeak && CurRightPeak ~= RightPeak)
#
#			ReparamX = zeros(NumNodes, 1);
#			ReparamX(1:CurLeftPeak) = linspace(1, LeftPeak, CurLeftPeak);
#			ReparamX(CurLeftPeak:CurRightPeak) = linspace(LeftPeak, RightPeak, CurRightPeak - CurLeftPeak + 1);
#			ReparamX(CurRightPeak:end) = linspace(RightPeak, NumNodes, NumNodes - CurRightPeak + 1);
#			RegisteredProfiles(:, z) = interp1(ReparamX, T, 1:NumNodes);
#			AllReparamX(:, z) = ReparamX;
#		%end
#		%keyboard;
#	end
#	%plot
#	%keyboard;
#end

def registerProfile(inputProfile):
	
	templateLeftPeak = numpy.round(inputProfile.size * 0.12) - 1
	templateRightPeak = numpy.round(inputProfile.size * 0.90) - 1
	#MeanPeaks = [LeftPeak, RightPeak];
	
	inputProfileNoNaNs = numpy.array(inputProfile)
	
	inputProfileNoNaNs[numpy.logical_not(numpy.isfinite(inputProfile))] = 0
	
	SE = numpy.ones((1, numpy.round(inputProfile.size / 5.0)))
	D = scipy.ndimage.morphology.grey_dilation(numpy.atleast_2d(inputProfileNoNaNs), footprint = SE)
	D = numpy.ravel(D)
#	pylab.subplot(1, 2, 1)
#	pylab.plot(inputProfile)
#	pylab.subplot(1, 2, 2)
#	pylab.plot(D)
#	print D
#	pylab.show()
	#thicknessProfile[numpy.logical_not(validStreamlines)] = numpy.nan
	
	peaks = numpy.where(inputProfileNoNaNs == D)[0]
	#print peaks
	leftPeaks = peaks[numpy.where(peaks < inputProfile.size / 2 - 1)]
	rightPeaks = peaks[numpy.where(peaks > inputProfile.size / 2 - 1)]
	
	#print leftPeaks
	#print rightPeaks
	
	if numpy.size(leftPeaks) > 0:
		I = numpy.argmax(inputProfileNoNaNs[leftPeaks])
		leftPeak = leftPeaks[I]
	else:
		leftPeak = templateLeftPeak
	
	if numpy.size(rightPeaks) > 0:
		I = numpy.argmax(inputProfileNoNaNs[rightPeaks])
		rightPeak = rightPeaks[I]
	else:
		rightPeak = templateRightPeak

	reparamX = numpy.zeros((inputProfile.size))
	reparamX[:leftPeak + 1] = numpy.linspace(0, templateLeftPeak, leftPeak + 1)
	reparamX[leftPeak:rightPeak + 1] = numpy.linspace(templateLeftPeak, templateRightPeak, rightPeak - leftPeak + 1)
	reparamX[rightPeak:] = numpy.linspace(templateRightPeak, inputProfile.size - 1, inputProfile.size - rightPeak)
	
	resampledProfile = numpy.interp(numpy.arange(inputProfile.size), reparamX, inputProfileNoNaNs, left=None, right=None)

	#pylab.plot(inputProfileNoNaNs)
	#pylab.plot(resampledProfile)

	#pylab.show()

	#quit()
	return resampledProfile
	#scipy.ndimage.morphology.grey_dilation(numpy.atleast_2d(T), footprint = numpy.ones((1, dilateSize)
	#pylab.plot(scaledContours['outer'][0], scaledContours['outer'][1], **lineProps)


def plotRectangle(Position, lineProps):
	
	assert(isinstance(lineProps, dict)),"LineProps must be a dict"

	pylab.plot([Position[0], Position[0]], [Position[1], Position[1] + Position[3]], **lineProps)
	pylab.plot([Position[0], Position[0] + Position[2]], [Position[1], Position[1]], **lineProps)
	pylab.plot([Position[0] + Position[2], Position[0] + Position[2]], [Position[1], Position[1] + Position[3]], **lineProps)
	pylab.plot([Position[0], Position[0] + Position[2]], [Position[1] + Position[3], Position[1] + Position[3]], **lineProps)

#def fn(val):
#	
#	print("fn says: %s : " % (val))
#	return numpy.size(val)

def thicknessCC(outputBase, groundTruthFile = None, numThicknessNodes = 100, doGraphics = False):
	
	assert(isinstance(numThicknessNodes, int) and numThicknessNodes > 0),"number of thickness nodes must be non-zero, positive and of type int"

	segMATFile = outputBase + "_seg.hdf5"
	segManeditMATFile = outputBase + "_seg_manedit.hdf5"
	
	assert(os.path.isfile(segMATFile)),"Seg mat file not found: " + segMATFile
	
	(head, tail) = os.path.split(outputBase)
	subjectID = tail

	FID = h5py.File(segMATFile, 'r')
	
	seg = dict()

	seg['IMG'] = numpy.array(FID['IMG']);
	#seg['initialSeg'] = numpy.array(FID['initialSeg']) > 0
	seg['finalSeg'] = numpy.array(FID['finalSeg']) > 0
	seg['finalSegNoArtefacts'] = numpy.array(FID['finalSegNoArtefacts']) > 0 # the greater than 0 is to convert it to a boolean array
	seg['templatePixdims'] = numpy.array(FID['templatePixdims']) 
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
		print "using man edit seg"
		FID = h5py.File(segManeditMATFile, 'r')
		finalSeg = numpy.array(FID['finalSegManEdit']) > 0
		FID.close()	
	else:
		finalSeg = numpy.array(seg['finalSeg'])
	
	if doGraphics:
		outputPNG = os.path.join(head, 'endpoint', subjectID + "_endpoints.png")
	else:
		outputPNG = None
	try:
		finalContours = LaplaceThicknessMethod.endpointsFind(finalSeg, outputPNG = outputPNG)
	except Exception as e:
		raise(e)
	
	#return

	xx, yy, scaledContours, streamlines, validStreamlines, solvedImage = LaplaceThicknessMethod.laplaceEquation2DContours(finalContours['xi'], finalContours['yi'], finalContours['xo'][::-1], finalContours['yo'][::-1], seg['templatePixdims'][0], seg['templatePixdims'][1], 1.0 / numpy.min(seg['templatePixdims']), 100, redoFactorisation = True)
	
	thicknessProfile = numpy.zeros((numThicknessNodes))
	for z in range(numThicknessNodes):
		D = numpy.diff(streamlines[z], axis = 1)

		thicknessProfile[z] = numpy.sum(numpy.sqrt(numpy.sum(D * D, axis = 0)))
	
	thicknessProfile[numpy.logical_not(validStreamlines)] = 0

	registeredThicknessProfile = registerProfile(thicknessProfile)

	if doGraphics:
		PNGDirectory = os.path.join(head, "thickness")
		try:
			os.makedirs(PNGDirectory)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(PNGDirectory):
				pass
			else:
				raise Exception

		outputPNG = os.path.join(PNGDirectory, subjectID + "_thickness.png")

		pylab.clf()
		
		pylab.subplot(2, 1, 1)
		for z in range(len(streamlines)):
			if validStreamlines[z] == True:
				lineProps = {'color': 'b', 'linewidth': 2}
			else:
				lineProps = {'color': 'm', 'linewidth': 2}
			CCSegUtils.plotContour(streamlines[z], closed = False, lineProps = lineProps)
			
			if numpy.mod(z, 5) == 0:
				pylab.text(streamlines[z][0, 0], streamlines[z][1, 0], str(z))
			#CCSegUtils.plotContour(streamlinesOuter[z], closed = False)
			#CCSegUtils.plotContour(streamlinesInner[z], closed = False)
		lineProps = {'color': 'g', 'linewidth': 2}
		pylab.plot(scaledContours['inner'][0], scaledContours['inner'][1], **lineProps)
		lineProps = {'color': 'r', 'linewidth': 2}
		pylab.plot(scaledContours['outer'][0], scaledContours['outer'][1], **lineProps)
		
		pylab.gca().invert_yaxis()
		pylab.gca().get_xaxis().set_ticks([])
		pylab.gca().get_yaxis().set_ticks([])
		pylab.title('Streamlines')

		pylab.subplot(2, 1, 2)
		thicknessProfilePlot, = pylab.plot(thicknessProfile)
		registeredThicknessProfilePlot, = pylab.plot(registeredThicknessProfile)
		pylab.legend([thicknessProfilePlot, registeredThicknessProfilePlot], ['Original thickness profile', 'Registered thickness profile'], loc = 9)
		pylab.title('Thickness profile')
		pylab.xlabel('Node')
		pylab.ylabel('Thickness (mm)')

		#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(FY, extent = [xx[0], xx[-1], yy[0], yy[-1]], ticks = True);
		#for z in range(numpy.size(sxi)):
		#	CCSegUtils.plotContour(streamlinesOuter[z], closed = False)
		#	CCSegUtils.plotContour(streamlinesInner[z], closed = False)
		#pylab.plot(scaledxi[0], scaledyi[0])
		#pylab.plot(scaledxo[0], scaledyo[0])
		
		#print scaledxi
		#print scaledyi
#	CCSegUtils.showIMG(FX)
		
		#print outputPNG
		pylab.savefig(outputPNG)
		CCSegUtils.cropAutoWhitePNG(outputPNG)
		
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()
	
	# perform Witelson and "Hofer and Frahm" parcellations on the CCs, output parcellations of the streamlines
	#../witelson_hoferfrahm_parcellation_schemes.png

	# first, get the leftmost and rightmost pixels

	maskClosed = numpy.isfinite(solvedImage)
	I = numpy.where(maskClosed)

	minCol = numpy.min(I[1])
	maxCol = numpy.max(I[1])

	minColRow = numpy.mean(I[0][numpy.where(I[1] == minCol)])  
	maxColRow = numpy.mean(I[0][numpy.where(I[1] == maxCol)])  
	pylab.clf()

	minColX = xx[minCol]
	minColY = yy[minColRow]
	maxColX = xx[maxCol]
	maxColY = yy[maxColRow]
	
	# get the slope of the line, rise over run
	cropLineM = (maxColY - minColY) / (maxColX - minColX)
	# get the Y intercept
	cropLineC = minColY - (minColX * cropLineM)

	#print str(cropLineM) + " * X + " + str(cropLineC)

	X, Y = numpy.meshgrid(xx, yy)
	
	midLineCross = Y - (cropLineM * X + cropLineC)
	#midLineCross = numpy.sign(midLineCross)
	
	witelsonLabels = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)
	hoferFrahmLabels = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)
	emsellLabels = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)
	
	I = numpy.where(maskClosed)
	midLineCrossInMask = midLineCross[I]
	
	midLineCrossInMaskIMG = numpy.array(midLineCross)
	midLineCrossInMaskIMG[numpy.logical_not(maskClosed)] = numpy.nan

	validX = X[I]
	validY = Y[I]
	
	# according to the witelson and HoferFrahm parcellations, the lines of division are perpendicular to the line that cuts the CC
	# find the scalar projection of the vector go)ing from (minColX, minColY) to (validX, validY) on (minColX, minColY) -> (maxColX, maxColY)
	# construct vectors A and B, the scalar projection is A . B / |B|^2
	AX = validX - minColX
	AY = validY - minColY
	
	BX = maxColX - minColX
	BY = maxColY - minColY
	
	antPostProportion = (AX * BX + AY * BY) / (BX * BX + BY * BY)
	del AX
	del AY
	del BX
	del BY
	
	# because our images are posterior to anterior, we need to invert the direction, i.e.:
	antPostProportion = 1.0 - antPostProportion

	witelsonI = numpy.zeros(validX.shape, dtype = numpy.uint8)
	
	witelsonI[numpy.where(antPostProportion < (1.0 / 3.0))] = 1	
	witelsonI[numpy.where(numpy.logical_and(antPostProportion < (1.0 / 2.0), witelsonI == 0))] = 2	
	witelsonI[numpy.where(numpy.logical_and(antPostProportion < (2.0 / 3.0), witelsonI == 0))] = 3	
	witelsonI[numpy.where(numpy.logical_and(antPostProportion < (4.0 / 5.0), witelsonI == 0))] = 4	
	witelsonI[numpy.where(witelsonI == 0)] = 5	
	
	witelsonI[midLineCrossInMask > 0] = 0
	
	#SIMG = numpy.zeros(maskClosed.shape)
	#SIMG[I] = witelsonI
	witelsonLabels[I] = witelsonI
	
	hoferFrahmI = numpy.zeros(validX.shape, dtype = numpy.uint8)
	
	hoferFrahmI[numpy.where(antPostProportion < (1.0 / 6.0))] = 1	
	hoferFrahmI[numpy.where(numpy.logical_and(antPostProportion < (1.0 / 2.0), hoferFrahmI == 0))] = 2	
	hoferFrahmI[numpy.where(numpy.logical_and(antPostProportion < (2.0 / 3.0), hoferFrahmI == 0))] = 3	
	hoferFrahmI[numpy.where(numpy.logical_and(antPostProportion < (3.0 / 4.0), hoferFrahmI == 0))] = 4	
	hoferFrahmI[numpy.where(hoferFrahmI == 0)] = 5
	hoferFrahmI[midLineCrossInMask > 0] = 0
	
	hoferFrahmLabels[I] = hoferFrahmI

	emsellI = numpy.zeros(validX.shape, dtype = numpy.uint8)
	
	emsellI[numpy.where(antPostProportion < (1.0 / 6.0))] = 2
	emsellI[numpy.where(numpy.logical_and(antPostProportion < (1.0 / 3.0), numpy.logical_and(emsellI == 0, midLineCrossInMask <= 0)))] = 3
	emsellI[numpy.where(numpy.logical_and(antPostProportion < (1.0 / 2.0), numpy.logical_and(emsellI == 0, midLineCrossInMask <= 0)))] = 4
	emsellI[numpy.where(numpy.logical_and(antPostProportion < (2.0 / 3.0), numpy.logical_and(emsellI == 0, midLineCrossInMask <= 0)))] = 5
	emsellI[numpy.where(numpy.logical_and(antPostProportion < (4.0 / 5.0), numpy.logical_and(emsellI == 0, midLineCrossInMask <= 0)))] = 6
	emsellI[numpy.where(numpy.logical_and(emsellI == 0, midLineCrossInMask <= 0))] = 7
	emsellI[numpy.where(numpy.logical_and(emsellI == 0, antPostProportion > (1.0 / 2.0)))] = 8
	
	
	# we need to find the apex of the genu
	
	# it is the pixel closest to the rightmost element of the inferior boundary
	# firstly find out which contour is the inferior contour
	# it is the contour whose minimum y coordinate is greatest
	if numpy.min(scaledContours['inner'][1]) > numpy.min(scaledContours['outer'][1]):
		inferiorContourXY = numpy.array(scaledContours['inner'])
	else:
		inferiorContourXY = numpy.array(scaledContours['outer'])
	
	#print scaledContours
	# find the rightmost point
	rightMostInnerXY = numpy.take(inferiorContourXY, [numpy.argmax(inferiorContourXY[0])], axis = 1)
	
	# get the anterior-posterior proportion that is closest to this point
	
	rightMostInnerXYClosestIDX = numpy.argmin(numpy.abs(rightMostInnerXY[0] - validX) + numpy.abs(rightMostInnerXY[1] - validY))
	
	antPostProportionClosest = antPostProportion[rightMostInnerXYClosestIDX]
	midLineCrossInMaskClosest = midLineCrossInMask[rightMostInnerXYClosestIDX]

	emsellI[numpy.where(numpy.logical_and(numpy.logical_and(antPostProportion > antPostProportionClosest, antPostProportion < (1.0 / 2.0)), midLineCrossInMask > midLineCrossInMaskClosest))] = 1


	################### THIS DOESNT ALWAYS WORK ######################## NEW METHOD ABOVE
	# it is the pixel with the highest S value in the right half on the 0 contour of midLineCrossInMask
	
	# if pixels on 0 contour of midLineCrossInMask
	# make midLineCrossInMask image
	#midLineCrossInMaskIMG = numpy.zeros(maskClosed.shape)
	#midLineCrossInMaskIMG.fill(numpy.inf)
	# make SIMG
	#IMG = numpy.zeros(maskClosed.shape)
	#IMG.fill(numpy.nan)

	#midLineCrossInMaskIMG[I] = midLineCrossInMask
	#IMG[I] = S
	
	# the roll shifts elements down one row, so this says pixelsAtSignChangeIMG[I, J] = midLineCrossInMaskIMG[I, J] == 1 && midLineCrossInMaskIMG[I - 1, J] == -1

	#pixelsAtSignChangeIMG = numpy.logical_and(midLineCrossInMaskIMG == 1, numpy.roll(midLineCrossInMaskIMG, 1, axis = 0) < 1)
	
	#pixelsAtSignChange = pixelsAtSignChangeIMG[I]

	# find highest S in the right half 
	#SAtPixelsAtSignChange = S[pixelsAtSignChange]

	#R = 1
	#SC = 1
	#pylab.subplot(SR, SC, 1)
	#CCSegUtils.showIMG(FY, extent = [xx[0], xx[-1], yy[0], yy[-1]], ticks = True);
	#CCSegUtils.showIMG(midLineCrossInMaskIMG, extent = [xx[0], xx[-1], yy[0], yy[-1]], ticks = True)
	#pylab.plot(inferiorContourXY[0, rightMostIDX], inferiorContourXY[1, rightMostIDX], 'r*')
#	pylab.colorbar()
#	pylab.subplot(SR, SC, 2)
#	CCSegUtils.showIMG(midLineCrossInMaskIMG == 1)
#	
#	pylab.subplot(SR, SC, 3)
#	CCSegUtils.showIMG(pixelsAtSignChangeIMG)
#	
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#quit()
#	
	
#	try:
#		SAtApex = numpy.max(SAtPixelsAtSignChange[SAtPixelsAtSignChange < (1.0 / 2.0)])
#		print "Proportion at rostrum apex: " + str(SAtApex)
#		emsellI[numpy.where(numpy.logical_and(numpy.logical_and(S > SAtApex, S < (1.0 / 2.0)), midLineCrossInMask > 0))] = 1
#	except Exception:
#		print "Rostrum apex not found, region 1 will be empty"
#
	emsellLabels[I] = emsellI
	
	# make RGB images
	witelsonRGB = numpy.tile(numpy.atleast_3d(numpy.double(maskClosed)), (1, 1, 3))
	hoferFrahmRGB = numpy.tile(numpy.atleast_3d(numpy.double(maskClosed)), (1, 1, 3))
	emsellRGB = numpy.tile(numpy.atleast_3d(numpy.double(maskClosed)), (1, 1, 3))
	
	#witelsonIMG = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)
	#hoferFrahmIMG = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)
	#emsellIMG = numpy.zeros(maskClosed.shape, dtype = numpy.uint8)

	#witelsonIMG[I] = witelsonLabels
	#hoferFrahmIMG[I] = hoferFrahmLabels
	#emsellIMG[I] = emsellLabels

	CMAP = numpy.array([[1, 1, 1], [1, 0, 0], [0, 0.5, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0.5, 0, 0], [0, 1, 1], [1, 0, 0.5]])
	
	#emsellLabels = numpy.mod(emsellLabels, CMAP.shape[0])
	

	for z in range(3):
		T = witelsonRGB[:, :, z]
		T[I] = CMAP[witelsonLabels[I], z]
		witelsonRGB[:, :, z] = T
		T = hoferFrahmRGB[:, :, z]
		T[I] = CMAP[hoferFrahmLabels[I], z]
		hoferFrahmRGB[:, :, z] = T
		T = emsellRGB[:, :, z]
		T[I] = CMAP[emsellLabels[I], z]
		emsellRGB[:, :, z] = T
		del T
	
	del I
	del witelsonI
	del hoferFrahmI
	del emsellI

	#CCSegUtils.showIMG(maskClosed, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	#pylab.subplot(2, 1, 1)
	#CCSegUtils.showIMG(witelsonRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	#pylab.title('Witelson parcellation')
	
	#pylab.subplot(2, 1, 2)
	#CCSegUtils.showIMG(hoferFrahmRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	#pylab.title('Hofer and Frahm parcellation')
	#pylab.plot(xx, cropLineM * xx + cropLineC, 'b-')

	#pylab.plot(minColX, minColY, 'g*')
	#pylab.plot(maxColX, maxColY, 'r*')

	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	#pylab.savefig(outputPNG)
		

	startV = numpy.zeros((2, numThicknessNodes))

	# make start points for each streamline
	for z in range(numThicknessNodes):
		D = numpy.diff(streamlines[z], axis = 1)
		cumArcLength = numpy.cumsum(numpy.sqrt(numpy.sum(D * D, axis = 0)))
		cumArcLength = cumArcLength / cumArcLength[-1]
		cumArcLength = numpy.concatenate((numpy.array([0]), cumArcLength))

		startV[0, z] = numpy.interp(0.5, cumArcLength, streamlines[z][0])
		startV[1, z] = numpy.interp(0.5, cumArcLength, streamlines[z][1])
	
	witelsonNodeLabels = CCSegUtils.interp2q(xx, yy, numpy.double(witelsonLabels), startV[0], startV[1], interpmethod = 'nearest')
	witelsonNodeLabels = numpy.uint8(witelsonNodeLabels)
	
	hoferFrahmNodeLabels = CCSegUtils.interp2q(xx, yy, numpy.double(hoferFrahmLabels), startV[0], startV[1], interpmethod = 'nearest')
	hoferFrahmNodeLabels = numpy.uint8(hoferFrahmNodeLabels)
	
	emsellNodeLabels = CCSegUtils.interp2q(xx, yy, numpy.double(emsellLabels), startV[0], startV[1], interpmethod = 'nearest')
	emsellNodeLabels = numpy.uint8(emsellNodeLabels)
	
	if doGraphics:
		pylab.clf()
		
		lineProps = {'color': 'y', 'linewidth': 2}
		
		SR = 3
		SC = 1
		pylab.subplot(SR, SC, 1)
		CCSegUtils.showIMG(witelsonRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
		for z in range(numThicknessNodes):
			pylab.text(startV[0, z], startV[1, z], str(witelsonNodeLabels[z]), horizontalalignment='center', verticalalignment='center')
		CCSegUtils.plotStreamlines(streamlines, lineProps = lineProps)
		pylab.title('Witelson parcellation')
		
		pylab.subplot(SR, SC, 2)
		CCSegUtils.showIMG(hoferFrahmRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
		for z in range(numThicknessNodes):
			pylab.text(startV[0, z], startV[1, z], str(hoferFrahmNodeLabels[z]), horizontalalignment='center', verticalalignment='center')
		CCSegUtils.plotStreamlines(streamlines, lineProps = lineProps)
		pylab.title('Hofer and Frahm parcellation')

		pylab.subplot(SR, SC, 3)
		CCSegUtils.showIMG(emsellRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
		for z in range(numThicknessNodes):
			pylab.text(startV[0, z], startV[1, z], str(emsellNodeLabels[z]), horizontalalignment='center', verticalalignment='center')
		CCSegUtils.plotStreamlines(streamlines, lineProps = lineProps)
		pylab.title('Emsell parcellation')

		parcellationPNGDirectory = os.path.join(head, "parcellations")
		try:
			os.makedirs(parcellationPNGDirectory)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(parcellationPNGDirectory):
				pass
			else:
				raise Exception

	
		outputPNG = os.path.join(parcellationPNGDirectory, subjectID + "_parcellations.png")
		
		pylab.gcf().set_size_inches((20, 10), forward = True)
		#pylab.show()
		#quit()

		pylab.savefig(outputPNG)
		CCSegUtils.cropAutoWhitePNG(outputPNG)
		
	pixelArea = (xx[1] - xx[0]) * (yy[1] - yy[0])
	# now do the stats
	witelsonStats = CCSegUtils.parcellationStats(witelsonLabels, pixelArea, thicknessProfile, witelsonNodeLabels)
	hoferFrahmStats = CCSegUtils.parcellationStats(hoferFrahmLabels, pixelArea, thicknessProfile, hoferFrahmNodeLabels)
	emsellStats = CCSegUtils.parcellationStats(emsellLabels, pixelArea, thicknessProfile, emsellNodeLabels)
	
	#print witelsonStats
	#print hoferFrahmStats

	#quit()
	outputMAT = outputBase + "_thickness.hdf5"

	FID = h5py.File(outputMAT, 'w')
	
	FID.create_dataset("xx", data = xx, compression = 'gzip')
	FID.create_dataset("yy", data = yy, compression = 'gzip')
	finalContoursGroup = FID.create_group("finalContours")
	finalContoursGroup.create_dataset("inner", data = scaledContours['inner'], compression = 'gzip')
	finalContoursGroup.create_dataset("outer", data = scaledContours['outer'], compression = 'gzip')
	
	#FID.create_dataset("startV", data = startV, compression = 'gzip')
	FID.create_dataset("validStreamlines", data = numpy.uint8(validStreamlines), compression = 'gzip')
	FID.create_dataset("thicknessProfile", data = thicknessProfile, compression = 'gzip')
	FID.create_dataset("registeredThicknessProfile", data = registeredThicknessProfile, compression = 'gzip')
	
	streamlinesGroup = FID.create_group("streamlines")
	for z in range(numThicknessNodes):
		streamlinesGroup.create_dataset(str(z), data = streamlines[z], compression = 'gzip')

	FID.create_dataset("solvedImage", data = solvedImage, compression = 'gzip')
	
	FID.create_dataset("witelsonLabels", data = witelsonLabels, compression = 'gzip')
	FID.create_dataset("witelsonNodeLabels", data = witelsonNodeLabels, compression = 'gzip')
	
	witelsonGroup = FID.create_group("witelsonStats")
	for curKey in witelsonStats.iterkeys():
		witelsonGroup.create_dataset(curKey, data = witelsonStats[curKey], compression = 'gzip')
	
	FID.create_dataset("hoferFrahmLabels", data = hoferFrahmLabels, compression = 'gzip')
	FID.create_dataset("hoferFrahmNodeLabels", data = hoferFrahmNodeLabels, compression = 'gzip')
	
	hoferFrahmGroup = FID.create_group("hoferFrahmStats")
	for curKey in hoferFrahmStats.iterkeys():
		hoferFrahmGroup.create_dataset(curKey, data = hoferFrahmStats[curKey], compression = 'gzip')
	
	FID.create_dataset("emsellLabels", data = emsellLabels, compression = 'gzip')
	FID.create_dataset("emsellNodeLabels", data = emsellNodeLabels, compression = 'gzip')
	
	emsellGroup = FID.create_group("emsellStats")
	for curKey in emsellStats.iterkeys():
		emsellGroup.create_dataset(curKey, data = emsellStats[curKey], compression = 'gzip')
	
	
	FID.create_dataset("startV", data = startV, compression = 'gzip')

	#FID.create_dataset("maskClosed", data = maskClosed, compression = 'gzip')

	FID.close()

	#Seg.FinalSeg = HDFReadFunc(SegMatFileName, '/finalSeg'); Seg.FinalSeg = logical(Seg.FinalSeg');
	#Seg.FinalSegArtefactsRemoved = hdf5read(SegMatFileName, '/finalSegNoArtefacts'); Seg.FinalSegArtefactsRemoved = logical(Seg.FinalSegArtefactsRemoved');
	#Seg.InitialSeg = hdf5read(SegMatFileName, '/initialSeg'); Seg.InitialSeg = logical(Seg.InitialSeg');
	#Seg.TemplatePixdims = hdf5read(SegMatFileName, '/templatePixdims'); Seg.TemplatePixdims = Seg.TemplatePixdims';
	
	#NIIPixdims = numpy.array(FID['NIIPixdims'])
	#midSagAVW = numpy.array(FID['midSagAVW'])
	#MSPMethod = str(numpy.array(FID['MSPMethod']))
	
	#midSagAVW = midSagAVW[:, ::-1]
	#midSagAVW = numpy.rot90(midSagAVW, -1)
