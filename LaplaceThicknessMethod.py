import numpy
import numpy.linalg

import os
import errno

import pylab
import scipy.ndimage
import scipy.sparse.linalg

import CCSegUtils

import cv2

from inspect import stack

import skimage.measure

try:
	from scikits.sparse.cholmod import cholesky
	usingCholesky = True
	CholeskyOfA = None
	LUOfA = None
except Exception:
	#print "Cholesky method not found, using LU decomposition, install scikits.sparse.cholmod"
	usingCholesky = False
	CholeskyOfA = None
	LUOfA = None


# extrema from MATLAB
def regionPropsExtrema(BW):
	
	pixelList = BW.nonzero()

	r = pixelList[0]
	c = pixelList[1]

	minR = numpy.min(r);
	maxR = numpy.max(r);
	minC = numpy.min(c);
	maxC = numpy.max(c);

	minRSet = (r == minR);
	maxRSet = (r == maxR);
	minCSet = (c == minC);
	maxCSet = (c == maxC);

	#% Points 1 and 2 are on the top row.
	r1 = minR;
	r2 = minR;
	#% Find the minimum and maximum column coordinates for
	#% top-row pixels.
	tmp = c[minRSet];
	c1 = numpy.min(tmp);
	c2 = numpy.max(tmp);

	#% Points 3 and 4 are on the right column.
	#% Find the minimum and maximum row coordinates for
	#% right-column pixels.
	tmp = r[maxCSet];
	r3 = numpy.min(tmp);
	r4 = numpy.max(tmp);
	c3 = maxC;
	c4 = maxC;

	#% Points 5 and 6 are on the bottom row.
	r5 = maxR;
	r6 = maxR;
	#% Find the minimum and maximum column coordinates for
	#% bottom-row pixels.
	tmp = c[maxRSet];
	c5 = numpy.max(tmp);
	c6 = numpy.min(tmp);

	#% Points 7 and 8 are on the left column.
	#% Find the minimum and maximum row coordinates for
	#% left-column pixels.
	tmp = r[minCSet];
	r7 = numpy.max(tmp);
	r8 = numpy.min(tmp);
	c7 = minC;
	c8 = minC;
	
	return numpy.array([[c1, r1], [c2, r2], [c3, r3], [c4, r4], [c5, r5], [c6, r6], [c7, r7], [c8, r8]])

def contourResampleArcLength(X, n):
	D = numpy.diff(X, axis = 1)

	cumArcLength = numpy.cumsum(numpy.sqrt(numpy.sum(D * D, axis = 0)))
	cumArcLength = cumArcLength / cumArcLength[-1]
	cumArcLength = numpy.concatenate((numpy.array([0]), cumArcLength))

	#print cumArcLength.shape

	outContour = numpy.zeros((X.shape[0], int(n)))
	
	for z in range(X.shape[0]):
		outContour[z] = numpy.interp(numpy.linspace(0, 1, n), cumArcLength, X[z])
	return outContour

# returns the arc length of the contour in X, dimensions are the rows
def arcLength(X):

	if X.shape[1] <= 1:
		return 0
	else:
		D = numpy.diff(X, axis = 1)

		return numpy.sum(numpy.sqrt(numpy.sum(D * D, axis = 0)))

def numericGradient2D(IMG, xSpacing = 1.0, ySpacing = 1.0, mask = None):
	
	if mask is None:
		T = numpy.gradient(IMG)
		return (T[1], T[0])
	else:
		assert(numpy.array_equal(IMG.shape, mask.shape)),"IMG and mask must be the same shape"

		maskIDX = numpy.where(mask)
		
		northIDX = (numpy.maximum(maskIDX[0] - 1,                 0), maskIDX[1])
		southIDX = (numpy.minimum(maskIDX[0] + 1, mask.shape[0] - 1), maskIDX[1])
		westIDX = (maskIDX[0], numpy.maximum(maskIDX[1] - 1,                 0))
		eastIDX = (maskIDX[0], numpy.minimum(maskIDX[1] + 1, mask.shape[1] - 1))
		
		northMask = mask[northIDX]
		notNorthMask = numpy.logical_not(northMask)
		southMask = mask[southIDX]
		notSouthMask = numpy.logical_not(southMask)
		westMask = mask[westIDX]
		notWestMask = numpy.logical_not(westMask)
		eastMask = mask[eastIDX]
		notEastMask = numpy.logical_not(eastMask)
		
		northSouthMask = numpy.logical_and(northMask, southMask)

		T = numpy.ones((numpy.count_nonzero(mask)))
		T[numpy.where(northSouthMask)] = 2.0
		northSouthSpacing = numpy.ones(mask.shape)
		northSouthSpacing[maskIDX] = T

		eastWestMask = numpy.logical_and(eastMask, westMask)
		T.fill(1)
		T[numpy.where(eastWestMask)] = 2.0
		eastWestSpacing = numpy.ones(mask.shape)
		eastWestSpacing[maskIDX] = T
		del T;

		# values that are north
		northValues = numpy.zeros(mask.shape)
		# the northValues that are in the mask should be taken from the north indices
		northValues[(maskIDX[0][northMask], maskIDX[1][northMask])] = IMG[(northIDX[0][northMask], northIDX[1][northMask])]
		# the northValues that are NOT in the mask should be taken from the current indices
		northValues[(maskIDX[0][notNorthMask], maskIDX[1][notNorthMask])] = IMG[(maskIDX[0][notNorthMask], maskIDX[1][notNorthMask])]
		
		# values that are south
		southValues = numpy.zeros(mask.shape)
		# the southValues that are in the mask should be taken from the south indices
		southValues[(maskIDX[0][southMask], maskIDX[1][southMask])] = IMG[(southIDX[0][southMask], southIDX[1][southMask])]
		# the southValues that are NOT in the mask should be taken from the current indices
		southValues[(maskIDX[0][notSouthMask], maskIDX[1][notSouthMask])] = IMG[(maskIDX[0][notSouthMask], maskIDX[1][notSouthMask])]
		# values that are west
		westValues = numpy.zeros(mask.shape)
		# the westValues that are in the mask should be taken from the west indices
		westValues[(maskIDX[0][westMask], maskIDX[1][westMask])] = IMG[(westIDX[0][westMask], westIDX[1][westMask])]
		# the westValues that are NOT in the mask should be taken from the current indices
		westValues[(maskIDX[0][notWestMask], maskIDX[1][notWestMask])] = IMG[(maskIDX[0][notWestMask], maskIDX[1][notWestMask])]
		
		# values that are east
		eastValues = numpy.zeros(mask.shape)
		# the eastValues that are in the mask should be taken from the east indices
		eastValues[(maskIDX[0][eastMask], maskIDX[1][eastMask])] = IMG[(eastIDX[0][eastMask], eastIDX[1][eastMask])]
		# the eastValues that are NOT in the mask should be taken from the current indices
		eastValues[(maskIDX[0][notEastMask], maskIDX[1][notEastMask])] = IMG[(maskIDX[0][notEastMask], maskIDX[1][notEastMask])]
		
#		pylab.subplot(2, 3, 1); CCSegUtils.showIMG(northValues)
#		pylab.subplot(2, 3, 2); CCSegUtils.showIMG(southValues)
#		pylab.subplot(2, 3, 3); CCSegUtils.showIMG((southValues - northValues) / northSouthSpacing)
#		pylab.subplot(2, 3, 4); CCSegUtils.showIMG(westValues)
#		pylab.subplot(2, 3, 5); CCSegUtils.showIMG(eastValues)
#		pylab.subplot(2, 3, 6); CCSegUtils.showIMG((eastValues - westValues) / eastWestSpacing)
#
#		pylab.gcf().set_size_inches((20, 10), forward = True)
#		pylab.show()
#		quit()
		return ((eastValues - westValues) / eastWestSpacing, (southValues - northValues) / northSouthSpacing)

#function [MaskInnerBoundary, ...
#	MaskOuterBoundary, ...
#	MaskFree, ...
#	SolvedImage, ...
#	XY, ...
#	X, ...
#	Y, ...
#	NormFX, ...
#	NormFY, ...
#	ValidStreamlines, ...
#	StartV, ...
#	MaskClosed] = laplace_get_points_2d_auto_mw(xi, yi, xo, yo, XScale, YScale, Delta, NumStreamlines, MaskClosedPrecomputed)

def laplaceEquationResetFactors():
	global CholeskyOfA
	global LUOfA
	CholeskyOfA = None
	LUOfA = None
	
# python port of laplace_get_points_2d_auto_mw.m
def laplaceEquation2DContours(xi, yi, xo, yo, xScale, yScale, delta, numStreamlines, onlyNeedStartV = False, redoFactorisation = False):
	displayMessages = False

	contoursBoundingBox = dict()
	
	contoursBoundingBox['minX'] = numpy.min(numpy.concatenate((xi, xo), axis = 1))
	contoursBoundingBox['maxX'] = numpy.max(numpy.concatenate((xi, xo), axis = 1))
	contoursBoundingBox['minY'] = numpy.min(numpy.concatenate((yi, yo), axis = 1))
	contoursBoundingBox['maxY'] = numpy.max(numpy.concatenate((yi, yo), axis = 1))
	
	xx = numpy.arange(numpy.floor(contoursBoundingBox['minX']) - 1, numpy.ceil(contoursBoundingBox['maxX']) + 1 + 1, 1.0 / delta) * xScale
	yy = numpy.arange(numpy.floor(contoursBoundingBox['minY']) - 1, numpy.ceil(contoursBoundingBox['maxY']) + 1 + 1, 1.0 / delta) * yScale
	
	#print xx
	#print yy

	X, Y = numpy.meshgrid(xx, yy)

	scaledxi = xi * xScale
	scaledyi = yi * yScale
	scaledxo = xo * xScale
	scaledyo = yo * yScale

	xSpacing = xx[1] - xx[0]
	ySpacing = yy[1] - yy[0]

	#if not onlyNeedStartV:
	#	print contoursBoundingBox
	#	print xx[0]
	#	print xx[-1]
	#	print yy[0]
	#	print yy[-1]
	# resample the contours veensely to rasterise them
	arcLengthI = arcLength(numpy.concatenate((scaledxi, scaledyi), axis = 0))
	arcLengthO = arcLength(numpy.concatenate((scaledxo, scaledyo), axis = 0))

	contourIResampled = contourResampleArcLength(numpy.concatenate((scaledxi, scaledyi), axis = 0), numpy.ceil(arcLengthI) * 10)
	contourOResampled = contourResampleArcLength(numpy.concatenate((scaledxo, scaledyo), axis = 0), numpy.ceil(arcLengthO) * 10)
	
	del arcLengthI; del arcLengthO;
	
	XI = numpy.atleast_2d(contourIResampled[0])
	YI = numpy.atleast_2d(contourIResampled[1])
	XO = numpy.atleast_2d(contourOResampled[0])
	YO = numpy.atleast_2d(contourOResampled[1])
	
	# find the nearest pixel of each of the XI, YI, XO, YO
	XIColumns = numpy.int32(numpy.round((XI - xx[0]) / xSpacing))
	XOColumns = numpy.int32(numpy.round((XO - xx[0]) / xSpacing))
	YIRows = numpy.int32(numpy.round((YI - yy[0]) / ySpacing))
	YORows = numpy.int32(numpy.round((YO - yy[0]) / ySpacing))
	
	#print X.shape
	maskInnerBoundary = numpy.zeros(X.shape, dtype = numpy.bool)
	maskOuterBoundary = numpy.zeros(X.shape, dtype = numpy.bool)

	maskInnerBoundary[(YIRows, XIColumns)] = True
	maskOuterBoundary[(YORows, XOColumns)] = True
	del YIRows; del YORows; del XIColumns; del XOColumns;
	del contourOResampled; del contourIResampled;

	arcLengthI = arcLength(numpy.concatenate((scaledxi, scaledyi), axis = 0))
	arcLengthO = arcLength(numpy.concatenate((scaledxo, scaledyo), axis = 0))
	#rint numpy.array([[scaledxi[0, 0]], [scaledyi[0, 0]]]).shape
	closedContour = numpy.concatenate((numpy.concatenate((scaledxi, scaledyi), axis = 0), numpy.concatenate((scaledxo, scaledyo), axis = 0), numpy.array([[scaledxi[0, 0]], [scaledyi[0, 0]]])), axis = 1)
	closedContourResampled = contourResampleArcLength(closedContour, numpy.ceil(arcLength(closedContour)) * 10)
	#print closedContour.shape
	#print closedContourResampled.shape
	
	closedColumns = numpy.int32(numpy.round((closedContourResampled[0] - xx[0]) / xSpacing))
	closedRows = numpy.int32(numpy.round((closedContourResampled[1] - yy[0]) / ySpacing))
	
	maskClosed = numpy.zeros(X.shape, dtype = numpy.bool)
	maskClosed[(closedRows, closedColumns)] = True
	
	del closedContour
	del closedContourResampled
	del closedColumns
	del closedRows

	#if not onlyNeedStartV:
	#	CCSegUtils.showIMG(maskClosed)
	#	CSegUtils.plotContour(closedContour)
	#	pylab.gcf().set_size_inches((20, 10), forward = True)
	#	pylab.show()
	#	quit()
	
	#maskClosed = scipy.ndimage.morphology.binary_fill_holes(numpy.logical_or(maskInnerBoundary, maskOuterBoundary), structure = numpy.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype = numpy.bool))
	maskClosed = scipy.ndimage.morphology.binary_fill_holes(maskClosed, structure = numpy.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype = numpy.bool))
	
	maskFree = numpy.logical_and(maskClosed, numpy.logical_not(numpy.logical_or(maskInnerBoundary, maskOuterBoundary)))
#	if not onlyNeedStartV:
#
#		pylab.subplot(2, 2, 1); CCSegUtils.showIMG(maskInnerBoundary, extent = [xx[0], xx[-1], yy[0], yy[-1]])#pylab.imshow(maskInnerBoundary, extent = [xx[0], xx[-1], yy[0], yy[-1]], origin = 'lower') #CCSegUtils.showIMG(maskInnerBoundary)
#		lineProps = {'color': 'r', 'linewidth': 2}
#		#CCSegUtils.plotContour(numpy.concatenate((XI, YI), axis = 0), lineProps = lineProps, closed = False);
#		CCSegUtils.plotContour(numpy.concatenate((scaledxi, scaledyi), axis = 0), lineProps = lineProps, closed = False);
#		lineProps = {'color': 'g', 'linewidth': 2}
#		CCSegUtils.plotContour(numpy.concatenate((scaledxo, scaledyo), axis = 0), lineProps = lineProps, closed = False);
#		pylab.subplot(2, 2, 2); CCSegUtils.showIMG(maskOuterBoundary, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#		lineProps = {'color': 'b', 'linewidth': 2}
#		CCSegUtils.plotContour(numpy.concatenate((XO, YO), axis = 0), lineProps = lineProps, closed = False);
#		pylab.subplot(2, 2, 3);	CCSegUtils.showIMG(numpy.logical_or(maskInnerBoundary, maskOuterBoundary), extent = [xx[0], xx[-1], yy[0], yy[-1]])
#		pylab.subplot(2, 2, 4);	CCSegUtils.showIMG(maskClosed, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#
#		pylab.gcf().set_size_inches((20, 10), forward = True)
#		pylab.show()
#		quit()
#	pylab.subplot(2, 2, 3); CCSegUtils.showIMG(numpy.logical_or(maskInnerBoundary, maskOuterBoundary), extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	
#	pylab.subplot(2, 2, 3); CCSegUtils.showIMG(maskClosed, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	pylab.subplot(2, 2, 4); CCSegUtils.showIMG(maskFree, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	
	del XI; del YI; del XO; del YO;
	
	# we now have maskFree (to be solved), maskInnerBoundary, maskOuterBoundary 
	################ SET UP LAPLACE EQUATIONS
	innerPotential = 0
	outerPotential = 1
	numVariables = numpy.sum(maskFree)
	variableIDXIMG = numpy.zeros(X.shape, dtype = numpy.int32)
	variableIDXIMG.fill(-1)
	
	# try to get the matlab ordering for comparison
	#variableIDXIMG = variableIDXIMG.T
	n = -1
	variableIDXIMG = numpy.rot90(variableIDXIMG, n)
	maskFree = numpy.rot90(maskFree, n)
	variableIDXIMG[numpy.where(maskFree)] = numpy.arange(numVariables)
	variableIDXIMG = numpy.rot90(variableIDXIMG, -n)
	maskFree = numpy.rot90(maskFree, -n)
	del n

	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(maskFree)
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(variableIDXIMG)
	
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()

	#quit()
	
	laplaceBVector = numpy.zeros((numVariables))
	
	#numNeighbours = numpy.zeros_like(laplaceBVector)

	#CCSegUtils.showIMG(variableIDXIMG, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	
	maskFreeIDX = numpy.where(maskFree)
	#T = numpy.concatenate((numpy.zeros_like(maskFreeIDX[0]), numpy.zeros_like(maskFreeIDX[0]) + 1, numpy.zeros_like(maskFreeIDX[0]) + 2))
	#maskFreeIDXForRGB = (numpy.tile(maskFreeIDX[0], (3)), numpy.tile(maskFreeIDX[1], (3)), T)
	#del T

	# try to do this in python
	northNeighbours = (maskFreeIDX[0] - 1, maskFreeIDX[1])
	southNeighbours = (maskFreeIDX[0] + 1, maskFreeIDX[1])
	westNeighbours = (maskFreeIDX[0], maskFreeIDX[1] - 1)
	eastNeighbours = (maskFreeIDX[0], maskFreeIDX[1] + 1)
	
	northNeighboursVariables = variableIDXIMG[northNeighbours]
	southNeighboursVariables = variableIDXIMG[southNeighbours]
	westNeighboursVariables = variableIDXIMG[westNeighbours]
	eastNeighboursVariables = variableIDXIMG[eastNeighbours]
	global CholeskyOfA
	global LUOfA
	global usingCholesky
	#print redoFactorisation
	if (CholeskyOfA == None and LUOfA == None) or redoFactorisation == True:
		#print "factorising"
		# we havent got anything precomputed so make the A matrix
		# make the laplaceAMatrix
		# we use a lil_matrix because we are changing the sparsity structure, the lil_matrix is efficient at this
		laplaceAMatrix = scipy.sparse.lil_matrix((numVariables, numVariables), dtype=numpy.double)

		northNeighboursOtherVariables = numpy.where(northNeighboursVariables != -1)
		southNeighboursOtherVariables = numpy.where(southNeighboursVariables != -1)
		westNeighboursOtherVariables = numpy.where(westNeighboursVariables != -1)
		eastNeighboursOtherVariables = numpy.where(eastNeighboursVariables != -1)
	
		# for the north/south neighbours of I
		#  if we see another variable J, set laplaceAMatrix[I, J] = -xSpacingSquared
		#  if we see one of the inner boundary pixels, add xSpacingSquared * innerPotential to laplaceBVector[I]
		#  if we see one of the outer boundary pixels, add xSpacingSquared * outerPotential to laplaceBVector[I]
		# neighbours are other variables
		# the I is the collection of variables whose neighbours have variable indices
		# the J is the indices of the other variables
		I = variableIDXIMG[(maskFreeIDX[0][northNeighboursOtherVariables], maskFreeIDX[1][northNeighboursOtherVariables])]
		J = northNeighboursVariables[northNeighboursOtherVariables]
		laplaceAMatrix[(I, J)] = -(xSpacing * xSpacing)
		del I; del J;

		# neighbours are other variables
		# the I is the collection of variables whose neighbours have variable indices
		# the J is the indices of the other variables
		I = variableIDXIMG[(maskFreeIDX[0][southNeighboursOtherVariables], maskFreeIDX[1][southNeighboursOtherVariables])]
		J = southNeighboursVariables[southNeighboursOtherVariables]
		laplaceAMatrix[(I, J)] = -(xSpacing * xSpacing)
		del I; del J;

		# for the east/west neighbours of I
		#  if we see another variable J, set laplaceAMatrix[I, J] = -ySpacingSquared
		#  if we see one of the inner boundary pixels, add ySpacingSquared * innerPotential to laplaceBVector[I]
		#  if we see one of the outer boundary pixels, add ySpacingSquared * outerPotential to laplaceBVector[I]

		# neighbours are other variables
		# the I is the collection of variables whose neighbours have variable indices
		# the J is the indices of the other variables
		I = variableIDXIMG[(maskFreeIDX[0][westNeighboursOtherVariables], maskFreeIDX[1][westNeighboursOtherVariables])]
		J = westNeighboursVariables[westNeighboursOtherVariables]
		laplaceAMatrix[(I, J)] = -(ySpacing * ySpacing)
		del I; del J;

		# neighbours are other variables
		# the I is the collection of variables whose neighbours have variable indices
		# the J is the indices of the other variables
		I = variableIDXIMG[(maskFreeIDX[0][eastNeighboursOtherVariables], maskFreeIDX[1][eastNeighboursOtherVariables])]
		J = eastNeighboursVariables[eastNeighboursOtherVariables]
		laplaceAMatrix[(I, J)] = -(ySpacing * ySpacing)
		del I; del J;

		# do the diagonal, it equals 2 * xSpacingSquared + 2 * ySpacingSquared
		laplaceAMatrix[(numpy.arange(numVariables), numpy.arange(numVariables))] = 2.0 * (xSpacing * xSpacing + ySpacing * ySpacing)
		
		if usingCholesky:
			CholeskyOfA = cholesky(laplaceAMatrix.tocsc())
			LUOfA = None
			#Z = L(laplaceBVector)
			#del L
		else:#	except Exception:
			LUOfA = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(laplaceAMatrix))
			CholeskyOfA = None
			#Z = LUofA.solve(laplaceBVector)
			#del LUofA
		
		del laplaceAMatrix
		del northNeighboursOtherVariables
		del southNeighboursOtherVariables
		del westNeighboursOtherVariables
		del eastNeighboursOtherVariables
	#elif laplaceCholeskyDecompPrecomputed != None and laplaceLUDecompPrecomputed == None:
	#	CholeskyOfA = laplaceCholeskyDecompPrecomputed
	#	LUOfA = None
	#elif laplaceCholeskyDecompPrecomputed == None and laplaceLUDecompPrecomputed != None:
	#	LUOfA = laplaceLUDecompPrecomputed
	#	CholeskyOfA = None

	northNeighboursInnerBoundary = numpy.where(maskInnerBoundary[northNeighbours])
	northNeighboursOuterBoundary = numpy.where(maskOuterBoundary[northNeighbours])
	
	# neighbours are on inner boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * innerPotential
	I = variableIDXIMG[(maskFreeIDX[0][northNeighboursInnerBoundary], maskFreeIDX[1][northNeighboursInnerBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (xSpacing * xSpacing) * innerPotential
	del I;

	# neighbours are on outer boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * outerPotential
	I = variableIDXIMG[(maskFreeIDX[0][northNeighboursOuterBoundary], maskFreeIDX[1][northNeighboursOuterBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (xSpacing * xSpacing) * outerPotential
	del I;
	del northNeighboursInnerBoundary; del northNeighboursOuterBoundary;

	# do the south neighbours
	southNeighboursInnerBoundary = numpy.where(maskInnerBoundary[southNeighbours])
	southNeighboursOuterBoundary = numpy.where(maskOuterBoundary[southNeighbours])

	# neighbours are on inner boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * innerPotential
	I = variableIDXIMG[(maskFreeIDX[0][southNeighboursInnerBoundary], maskFreeIDX[1][southNeighboursInnerBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (xSpacing * xSpacing) * innerPotential
	del I;

	# neighbours are on outer boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * outerPotential
	I = variableIDXIMG[(maskFreeIDX[0][southNeighboursOuterBoundary], maskFreeIDX[1][southNeighboursOuterBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (xSpacing * xSpacing) * outerPotential
	del I;
	del southNeighboursInnerBoundary; del southNeighboursOuterBoundary;
	
	# do the west neighbours
	westNeighboursInnerBoundary = numpy.where(maskInnerBoundary[westNeighbours])
	westNeighboursOuterBoundary = numpy.where(maskOuterBoundary[westNeighbours])

	# neighbours are on inner boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * innerPotential
	I = variableIDXIMG[(maskFreeIDX[0][westNeighboursInnerBoundary], maskFreeIDX[1][westNeighboursInnerBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (ySpacing * ySpacing) * innerPotential
	del I;

	# neighbours are on outer boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * outerPotential
	I = variableIDXIMG[(maskFreeIDX[0][westNeighboursOuterBoundary], maskFreeIDX[1][westNeighboursOuterBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (ySpacing * ySpacing) * outerPotential
	del I;
	del westNeighboursInnerBoundary; del westNeighboursOuterBoundary;
	
	# do the east neighbours
	eastNeighboursInnerBoundary = numpy.where(maskInnerBoundary[eastNeighbours])
	eastNeighboursOuterBoundary = numpy.where(maskOuterBoundary[eastNeighbours])

	# neighbours are on inner boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * innerPotential
	I = variableIDXIMG[(maskFreeIDX[0][eastNeighboursInnerBoundary], maskFreeIDX[1][eastNeighboursInnerBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (ySpacing * ySpacing) * innerPotential
	del I;

	# neighbours are on outer boundary
	# the I is the collection of free indices that are on the boundary
	# increment the B vector xSpacingSquared * outerPotential
	I = variableIDXIMG[(maskFreeIDX[0][eastNeighboursOuterBoundary], maskFreeIDX[1][eastNeighboursOuterBoundary])]
	laplaceBVector[(I)] = laplaceBVector[(I)] + (ySpacing * ySpacing) * outerPotential
	del I;
	del eastNeighboursInnerBoundary; del eastNeighboursOuterBoundary;

	del northNeighbours; del southNeighbours; del eastNeighbours; del westNeighbours;
		
	#import time	
	#t = time.time()
	# do stuff
	#print laplaceAMatrix.shape
	#L = scipy.linalg.cholesky(laplaceAMatrix.todense())
	#LUofA = scipy.sparse.linalg.splu(scipy.sparse.csc_matrix(laplaceAMatrix))
	#print LS
	#print laplaceBVector
	#elapsed = time.time() - t
	#print "LU: " + str(elapsed)
	#print usingCholesky	

	#elapsed = time.time() - t
	#print "CHOL: " + str(elapsed)
	#import LinAlgHelpers

	# my implementation was too slow
	#print "Starting Chol"
	#t = time.time()
	#ZL = LinAlgHelpers.chol(laplaceAMatrix)
	#print "Finished Chol"

	#Z1 = LinAlgHelpers.substituteSolve(ZL, laplaceBVector)
	#Z = LinAlgHelpers.substituteSolve(ZL.T, Z1)
	#elapsed = time.time() - t
	#print "My CHOL: " + str(elapsed)
	
	
	#laplaceCholeskyDecompPrecomputed = None, laplaceLUDecompPrecomputed = None
	#quit()
	if usingCholesky:
		#if laplaceCholeskyDecompPrecomputed != None:
		#	CholeskyOfA

		Z = CholeskyOfA(laplaceBVector)
		#del L
	else:
		#if laplaceLUDecompPrecomputed != None:
		#	LUOfA = laplaceLUDecompPrecomputed
		#print LUOfA
		#print laplaceBVector
		#print laplaceLUDecompPrecomputed
		
		Z = LUOfA.solve(laplaceBVector)
		#del LUOfA
	#if not onlyNeedStartV:
	#	print Z	
	del laplaceBVector

	solvedImage = numpy.zeros(maskFree.shape)
	solvedImage[numpy.where(maskInnerBoundary)] = innerPotential
	solvedImage[numpy.where(maskOuterBoundary)] = outerPotential
	solvedImage[numpy.where(maskFree)] = Z[variableIDXIMG[numpy.where(maskFree)]]
	solvedImage[numpy.where(numpy.logical_not(maskClosed))] = numpy.nan
	#if not onlyNeedStartV:
	#	solvedImage[numpy.where(numpy.logical_not(maskFree))] = 0
	#print solvedImage[numpy.where(maskFree)]
	#if not onlyNeedStartV:
	#	CCSegUtils.showIMG(solvedImage, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	#	pylab.colorbar()
	#	pylab.gcf().set_size_inches((20, 10), forward = True)
	#	pylab.show()
	#	quit()

	# the gradients get normalised but this is to take into account possible different x to y sizes
	FX, FY = numericGradient2D(solvedImage, xSpacing = xSpacing, ySpacing = ySpacing, mask = maskClosed)
	
	# normalise the gradients
	gradMAG = numpy.sqrt(FX * FX + FY * FY)
	gradMAG[numpy.where(gradMAG == 0)] = 1

	FX = FX / gradMAG
	FY = FY / gradMAG
	del gradMAG
	
	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(FX, extent = [xx[0], xx[-1], yy[0], yy[-1]])
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(FY, extent = [xx[0], xx[-1], yy[0], yy[-1]])

	#if numStreamlines > 0:
	C = skimage.measure.find_contours(solvedImage, 0.5)
	outC = list()

	for z in range(len(C)):
		if C[z].shape[0] > 2:
			#print C[z].shape
			I = numpy.where(numpy.all(numpy.isfinite(C[z]), axis = 1))
		#	print I
			T = numpy.atleast_2d(numpy.flipud(numpy.squeeze(numpy.take(C[z], I, axis = 0)).T))
			if T.shape[1] > 2:
				outC.append(T)
			del T
			del I
	
	if len(outC) > 1:
		# choose the largest one

		#print "Warning more than one isocontour, picking largest one"
		maxSZ = 0
		maxIDX = 0

		for z in range(len(outC)):
			if outC[z].shape[1] > maxSZ:
				maxSZ = outC[z].shape[1]
				maxIDX = z
		
		startV = numpy.array(outC[maxIDX])
	else:
		startV = numpy.array(outC[0])

	#	print len(outC)
	#	for z in range(len(outC)):
	#		print outC[z].shape
	#	pylab.plot(outC[z][0], outC[z][1])
		#CCSegUtils.showIMG(solvedImage)
	#	pylab.gcf().set_size_inches((20, 10), forward = True)
	#	pylab.show()
	
	
	#assert(len(outC) == 1),"More than one midpoint isocontour"

	
	if onlyNeedStartV:
		return startV
	
	#rint startV
	del outC; del C;

	startV = contourResampleArcLength(startV, numStreamlines)
	
	startV[0] = numpy.interp(startV[0], numpy.arange(numpy.size(xx)), xx)
	startV[1] = numpy.interp(startV[1], numpy.arange(numpy.size(yy)), yy)
	
	#% fix to avoid oscillation at boundaries
	#% pad the gradient at the boundaries
	SE = numpy.ones((3, 3), dtype = numpy.bool)
	
	I = numpy.where(maskClosed)
	oldFX = numpy.array(FX)
	oldFY = numpy.array(FY)
	S = scipy.ndimage.morphology.grey_dilation(FX, footprint = SE)
	S[I] = 0;
	T = scipy.ndimage.morphology.grey_erosion(FX, footprint = SE)
	T[I] = 0;
	FX = FX + S + T;

	S = scipy.ndimage.morphology.grey_dilation(FY, footprint = SE);
	S[I] = 0;
	T = scipy.ndimage.morphology.grey_erosion(FY, footprint = SE)
	T[I] = 0;
	FY = FY + S + T;

#	pylab.subplot(2, 2, 1); CCSegUtils.showIMG(oldFX, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	lineProps = {'color': 'r', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XI, YI), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'b', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XO, YO), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'g', 'linewidth': 2}
#	CCSegUtils.plotContour(startV, lineProps = lineProps, closed = False);
#	pylab.subplot(2, 2, 2); CCSegUtils.showIMG(oldFY, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	lineProps = {'color': 'r', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XI, YI), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'b', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XO, YO), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'g', 'linewidth': 2}
#	CCSegUtils.plotContour(startV, lineProps = lineProps, closed = False);
#
#	
	#pylab.subplot(2, 2, 3); CCSegUtils.showIMG(FX, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	lineProps = {'color': 'r', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XI, YI), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'b', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XO, YO), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'g', 'linewidth': 2}
#	CCSegUtils.plotContour(startV, lineProps = lineProps, closed = False);
#	pylab.subplot(2, 2, 4); CCSegUtils.showIMG(FY, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#	lineProps = {'color': 'r', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XI, YI), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'b', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((XO, YO), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'g', 'linewidth': 2}
#	CCSegUtils.plotContour(startV, lineProps = lineProps, closed = False);
	#for z in range(len(outC)):
		#print outC[z].shape
	
	#CCSegUtils.showIMG(solvedImage)
	#CCSegUtils.plotContour(outC[0], closed = False)
	
	## STREAMLINE PART ###

	sxi = numpy.interp(startV[0], xx, numpy.arange(numpy.size(xx)))
	syi = numpy.interp(startV[1], yy, numpy.arange(numpy.size(yy)))
	
	import Streamline2DCython
	
	#pylab.subplot(2, 2, 1); CCSegUtils.showIMG(oldFX)
	#pylab.subplot(2, 2, 2); CCSegUtils.showIMG(oldFY)
	#pylab.subplot(2, 2, 3); CCSegUtils.showIMG(FX); pylab.plot(sxi, syi, 'g-');
	#pylab.subplot(2, 2, 4); CCSegUtils.showIMG(FY); pylab.plot(sxi, syi, 'g-');
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()

	#print xx.shape
	#print yy.shape
	#print FX.shape
	#print FY.shape
	
	#print xx[0]
	#print xx[-1]
	streamlinesOuter = list()	
	streamlinesInner = list()
	streamlinesOuterIntersected = list()	
	streamlinesInnerIntersected = list()	
	#streamlinesOuter = [None,] * 100	
	#streamlinesInner = [None,] * 100

	#z = 50
	#streamlinesInner[z] = Streamline2DCython.streamline2D(xx, yy, -FX, -FY, sxi[z], syi[z], 0.2, 10000, numpy.ravel(scaledxi), numpy.ravel(scaledyi))
	

	#streamline = [None,] * 100
	for z in range(numpy.size(sxi)):
	#for z in numpy.array([0.000, 14.000, 22.000, 78.000, 99.000], dtype = numpy.int):
		#print "Contour " + str(z)
		T, I = Streamline2DCython.streamline2D(xx, yy, -FX, -FY, sxi[z], syi[z], 0.2, 10000, numpy.ravel(scaledxi), numpy.ravel(scaledyi))
		streamlinesOuter.append(numpy.array(T))
		streamlinesOuterIntersected.append(I)
		del I
		del T

		T, I = Streamline2DCython.streamline2D(xx, yy, FX, FY, sxi[z], syi[z], 0.2, 10000, numpy.ravel(scaledxo), numpy.ravel(scaledyo))
		streamlinesInner.append(numpy.array(T))
		streamlinesInnerIntersected.append(I)
		del I
		del T

	streamlinesOuterIntersected = numpy.array(streamlinesOuterIntersected)
	streamlinesOuterIntersected = (streamlinesOuterIntersected == 1)
	streamlinesInnerIntersected = numpy.array(streamlinesInnerIntersected)
	streamlinesInnerIntersected = (streamlinesInnerIntersected == 1)
	
	validStreamlines = numpy.logical_and(streamlinesOuterIntersected, streamlinesInnerIntersected)
	
	streamlines = list()
	
	#print len(streamlinesInner)
	#print len(streamlinesOuter)

	for z in range(len(streamlinesInner)):
		streamlines.append(numpy.concatenate((numpy.fliplr(streamlinesInner[z][:, 1:]), streamlinesOuter[z]), axis = 1))
	
#	CCSegUtils.showIMG(FX, extent = [xx[0], xx[-1], yy[0], yy[-1]], ticks = True);
	#for z in range(len(streamlines)):
	#	lineProps = {'color': 'b', 'linewidth': 2}
	#	CCSegUtils.plotContour(streamlines[z], closed = False, lineProps = lineProps)
		#CCSegUtils.plotContour(streamlinesOuter[z], closed = False)
		#CCSegUtils.plotContour(streamlinesInner[z], closed = False)
	#lineProps = {'color': 'g', 'linewidth': 2}
	#pylab.plot(scaledxi[0], scaledyi[0], **lineProps)
	#lineProps = {'color': 'r', 'linewidth': 2}
	#pylab.plot(scaledxo[0], scaledyo[0], **lineProps)
#
#	pylab.subplot(1, 2, 2); CCSegUtils.showIMG(FY, extent = [xx[0], xx[-1], yy[0], yy[-1]], ticks = True);
#	for z in range(numpy.size(sxi)):
#		CCSegUtils.plotContour(streamlinesOuter[z], closed = False)
#		CCSegUtils.plotContour(streamlinesInner[z], closed = False)
#	pylab.plot(scaledxi[0], scaledyi[0])
#	pylab.plot(scaledxo[0], scaledyo[0])
	
	#print scaledxi
	#print scaledyi
#	CCSegUtils.showIMG(FX)
	#pylab.gca().invert_yaxis()
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	
	scaledContours = dict()
	scaledContours['inner'] = numpy.concatenate((scaledxi, scaledyi), axis = 0)
	scaledContours['outer'] = numpy.concatenate((scaledxo, scaledyo), axis = 0)
	
	return (xx, yy, scaledContours, streamlines, validStreamlines, solvedImage)
	#quit()

	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	

	#quit()

	# for the 

	
	#print I
	#maskRGB = numpy.double(numpy.concatenate((numpy.atleast_3d(maskFree), numpy.atleast_3d(maskInnerBoundary), numpy.atleast_3d(maskOuterBoundary)), axis = 2))
#
	#pylab.subplot(2, 2, 1); CCSegUtils.showIMG(maskRGB, extent = [xx[0], xx[-1], yy[0], yy[-1]])
#
	#T = numpy.array(maskRGB)
	#B = numpy.zeros_like(maskFree)
	#for z in range(len(I)):
	#	B = numpy.logical_or(B, variableIDXIMG == I[(z)])
	
	#B = B[numpy.where(maskFree)]

	#T[maskFreeIDXForRGB] = numpy.tile(B, (3))
	#pylab.subplot(2, 2, 2); CCSegUtils.showIMG(T, extent = [xx[0], xx[-1], yy[0], yy[-1]]); pylab.title('Other neighbours')
#	
#	T = numpy.array(maskRGB)
#	T[maskFreeIDXForRGB] = numpy.tile(innerBoundaryIDX, (3))
#	pylab.subplot(2, 2, 3); CCSegUtils.showIMG(T, extent = [xx[0], xx[-1], yy[0], yy[-1]]); pylab.title('Inner boundary')
#	
#	T = numpy.array(maskRGB)
#	T[maskFreeIDXForRGB] = numpy.tile(outerBoundaryIDX, (3))
#	pylab.subplot(2, 2, 4); CCSegUtils.showIMG(T, extent = [xx[0], xx[-1], yy[0], yy[-1]]); pylab.title('Inner boundary')
#
	
		


	# for the east/west neighbours of I
	#  if we see another variable J, set laplaceAMatrix[I, J] = -ySpacingSquared
	#  if we see one of the inner boundary pixels, add ySpacingSquared * innerPotential to laplaceBVector[I]
	#  if we see one of the outer boundary pixels, add ySpacingSquared * outerPotential to laplaceBVector[I]
	
	# for the north/south neighbours
	#  if we see another variable, set laplaceAMatrix[I, J] = -ySpacingSquared

	#  if we see for the north neighbours, if we see another variable, set laplaceAMatrix[I, J] = -xSpacingSquared
	# numNeighbours = 2 * xSpacingSquared + 2 * ySpacingSquared



	#closedHighResContour = numpy.concatenate
	#print XI.shape
	#print YI
	#del contourIResampled; del contourOResampled;
	#pylab.plot(BWContour[0, leftStart + curLeftOffset], BWContour[1, leftStart + curLeftOffset], 'b*');
	#print curContours['xi']
	#pylab.gca().invert_yaxis()
	
	#quit()


def endpointsFindObjectiveFunction(BWContour, leftOffset, rightOffset, rightStartIDX):
	# make the contour start from the current leftOffset
	
	BWContourLeftOffset	= numpy.roll(BWContour, -leftOffset, axis = 1)
	
	curContours = dict()
	
	#print BWContourLeftOffset.shape
	# the "inner" contour is from the start to the current left offset, we subtract the left offset because we have shifted the contour
	curContours['xi'] = numpy.atleast_2d(BWContourLeftOffset[0, 0:(rightStartIDX - leftOffset + rightOffset + 1)])
	curContours['yi'] = numpy.atleast_2d(BWContourLeftOffset[1, 0:(rightStartIDX - leftOffset + rightOffset + 1)])
	# the "outer" contour is from the end of the inner contour to the end
	curContours['xo'] = numpy.atleast_2d(BWContourLeftOffset[0, (rightStartIDX - leftOffset + rightOffset + 1):])
	curContours['yo'] = numpy.atleast_2d(BWContourLeftOffset[1, (rightStartIDX - leftOffset + rightOffset + 1):])
	
	#if rightOffset > 1:
	#	lineProps = {'color': 'r', 'linewidth': 2}
	#	CCSegUtils.plotContour(numpy.concatenate((curContours['xi'], curContours['yi']), axis = 0), lineProps = lineProps, closed = False);
	#	lineProps = {'color': 'b', 'linewidth': 2}
	#	CCSegUtils.plotContour(numpy.concatenate((curContours['xo'], curContours['yo']), axis = 0), lineProps = lineProps, closed = False);
		#pylab.plot(BWContour[0, leftStart + curLeftOffset], BWContour[1, leftStart + curLeftOffset], 'b*');
		#print curContours['xi']
		#pylab.gcf().set_size_inches((20, 10), forward = True)
	#	pylab.title(str(rightOffset))
	#	pylab.show()
		#quit()
	
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(maskClosed);
	
	startV = laplaceEquation2DContours(curContours['xi'], curContours['yi'], curContours['xo'][::-1], curContours['yo'][::-1], 1.0, 1.0, 1.0, 0, onlyNeedStartV = True)
	
	#print startV
	arcLength = numpy.sum(numpy.sqrt(numpy.sum(numpy.diff(startV, axis = 1) * numpy.diff(startV, axis = 1), axis = 0)))
	#print arcLength
	return arcLength

def contourClockwise(BW, BWContour):

	paddedBW = numpy.pad(BW, (1, 1), mode = 'constant', constant_values = 0)
	
	TBWContour = BWContour + 1

	T = numpy.roll(BWContour, -1, axis = 1) - BWContour
	#T = T + 1
	
	clockwiseArray = numpy.zeros((BWContour.shape[1]))

	for z in range(BWContour.shape[1]):
		# get the orientations
		if T[0, z] < 0 and T[1, z] < 0: # up-left
			if paddedBW[TBWContour[1, z] + 1, TBWContour[0, z] - 1]:
				# the left hand side of this vector is down-left
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] - 1, TBWContour[0, z] + 1]:
				# the right hand side of this vector is up-right
				clockwiseArray[z] = 1
		elif T[0, z] == 0 and T[1, z] < 0: # up
			if paddedBW[TBWContour[1, z], TBWContour[0, z] - 1]:
				# the left hand side of this vector is left
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z], TBWContour[0, z] + 1]:
				# the right hand side of this vector is right
				clockwiseArray[z] = 1
		elif T[0, z] > 0 and T[1, z] < 0: # up-right
			if paddedBW[TBWContour[1, z] - 1, TBWContour[0, z] - 1]:
				# the left hand side of this vector is up-left
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] + 1, TBWContour[0, z] + 1]:
				# the right hand side of this vector is down-right
				clockwiseArray[z] = 1
		elif T[0, z] < 0 and T[1, z] == 0: # left
			if paddedBW[TBWContour[1, z] + 1, TBWContour[0, z]]:
				# the left hand side of this vector is down
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] - 1, TBWContour[0, z]]:
				# the right hand side of this vector is up
				clockwiseArray[z] = 1
		elif T[0, z] > 0 and T[1, z] == 0: # right
			if paddedBW[TBWContour[1, z] - 1, TBWContour[0, z]]:
				# the left hand side of this vector is up
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] + 1, TBWContour[0, z]]:
				# the right hand side of this vector is down
				clockwiseArray[z] = 1
		elif T[0, z] < 0 and T[1, z] > 0: # down-left
			if paddedBW[TBWContour[1, z] + 1, TBWContour[0, z] + 1]:
				# the left hand side of this vector is down-right
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] - 1, TBWContour[0, z] - 1]:
				# the right hand side of this vector is up-left
				clockwiseArray[z] = 1
		elif T[0, z] == 0 and T[1, z] > 0: # down
			if paddedBW[TBWContour[1, z], TBWContour[0, z] + 1]:
				# the left hand side of this vector is right
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z], TBWContour[0, z] - 1]:
				# the right hand side of this vector is left
				clockwiseArray[z] = 1
		elif T[0, z] > 0 and T[1, z] > 0: # down-right
			if paddedBW[TBWContour[1, z] - 1, TBWContour[0, z] + 1]:
				# the left hand side of this vector is up-right
				clockwiseArray[z] = -1
			elif paddedBW[TBWContour[1, z] + 1, TBWContour[0, z] - 1]:
				# the right hand side of this vector is down-left
				clockwiseArray[z] = 1
	
	#print clockwiseArray
	#CCSegUtils.showIMG(BW)
	#pylab.plot(BWContour[0], BWContour[1], 'r-')
	
	#I = numpy.where(clockwiseArray == -1); pylab.quiver(BWContour[0, I], BWContour[1, I], T[0, I], T[1, I], color = 'r', angles = 'xy')
	#I = numpy.where(clockwiseArray == 1); pylab.quiver(BWContour[0, I], BWContour[1, I], T[0, I], T[1, I], color = 'g', angles = 'xy')
	#pylab.show()
	#quit()
	clockwiseVotes = numpy.size(numpy.where(clockwiseArray == 1))
	anticlockwiseVotes = numpy.size(numpy.where(clockwiseArray == -1))
	if clockwiseVotes > anticlockwiseVotes:
		return "clockwise"
	elif clockwiseVotes < anticlockwiseVotes:
		return "anticlockwise"
	else:
		return None

		
def endpointsFind(BW, outputPNG = None):
	assert(BW.dtype == numpy.bool and isinstance(BW, numpy.ndarray)),"BW must be a boolean numpy.ndarray"

	BWLabels, BWNumLabels = scipy.ndimage.measurements.label(BW, structure = numpy.ones([3, 3]))

	assert(BWNumLabels == 1),"Segmentation has multiple regions"

	del BWLabels; del BWNumLabels;
	
	BWCols = numpy.nonzero(BW)[1]
	midCol = int(numpy.floor((numpy.min(BWCols) + numpy.max(BWCols)) / 2.0))
	BWLeft = numpy.array(BW)
	BWRight = numpy.array(BW)
	
	BWLeft[:, midCol:] = False
	BWRight[:, 0:midCol] = False
	
	BWLeftExtrema = regionPropsExtrema(BWLeft)
	BWRightExtrema = regionPropsExtrema(BWRight)
	
	leftStart = numpy.squeeze(numpy.take(BWLeftExtrema, [4], axis = 0))
	rightStart = numpy.squeeze(numpy.take(BWRightExtrema, [5], axis = 0))
	
	#print BWLeftExtrema
	#print BWRightExtrema

	#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(BWLeft); pylab.plot(BWLeftExtrema[:, 0], BWLeftExtrema[:, 1], 'b*');
	#% the desirable extrema on the left half is bottom-right (5)
	#print leftStart
	#print leftStart.shape
	#pylab.plot(leftStart[0], leftStart[1], markersize=10, color='r', marker = '*')
	#pylab.subplot(1, 2, 2); CCSegUtils.showIMG(BWRight); pylab.plot(BWRightExtrema[:, 0], BWRightExtrema[:, 1], 'b*'); 
	#print rightStart
	#pylab.plot(rightStart[0], rightStart[1], markersize=10, color='r', marker = '*')
	#% the desirable extrema on the right half is bottom-left (6)
	
	del BWLeftExtrema; del BWRightExtrema;
	
	BWContours, hierarchy = cv2.findContours(numpy.uint8(BW), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
	del hierarchy

	BWContour = numpy.squeeze(BWContours[0]).T
	
	if contourClockwise(BW, BWContour) == "anticlockwise":
		BWContour = numpy.fliplr(BWContour)

	del BWContours
	
	# the matlab version worked out the closest point on the boundary contour to leftStart using the Euclidean distance
	# we will use the chessboard distance
	I = numpy.argmin(numpy.sum(numpy.abs(numpy.atleast_2d(leftStart).T - BWContour), axis = 0))
	
	# roll (circshift) the contour so leftStart is the first element
	BWContour = numpy.roll(BWContour, -I, axis = 1)
	
	# now find the contour element closest to the right start extrema
	rightStartIDX = numpy.argmin(numpy.sum(numpy.abs(numpy.atleast_2d(rightStart).T - BWContour), axis = 0))
	
	# compute the search range for the left point
	# indices increasing ->
	#               leftLimitPositive
	# -----------------|------------------------------
	# |                                              |
	# |                |< elementsFromMidColPoint >  |
	# |                                              |
	# + leftStart (contour index 0)                midCol
	# |                                              |
	# |                |< elementsFromMidColPoint >  |
	# |                                              |
	# -----------------|------------------------------
	#               leftLimitNegative
	# <- indices increasing (end of contour)
	elementsFromMidcolPoint = 20
	I = numpy.where(BWContour[0] >= midCol)[0]
	#print I
	
	# so the left limit positive is the first element that is greater than midcol
	leftLimitPositive = I[0] - elementsFromMidcolPoint - 1
	# so the left limit negative is the last element that is greater than midcol
	leftLimitNegative = -I[-1] + elementsFromMidcolPoint + 1
	
	del I

	BWContourRightStart = numpy.roll(BWContour, -rightStartIDX, axis = 1)

	# compute the search range for the right point
	#         indices increasing -> (end of contour)
	#                         rightLimitNegative
	# ------------------------------|---------------
	# |                                            |
	# |< elementsFromMidColPoint >  |              |
	# |                                            |
	# midCol                                       + rightStart (contour index 0)
	# |                                            |
	# |< elementsFromMidColPoint >  |              |
	# |                                            |
	# ------------------------------|---------------
	#               rightLimitPositive
	#                          <- indices increasing

	I = numpy.where(BWContourRightStart[0] <= midCol)[0]
	
	rightLimitPositive = I[0] - elementsFromMidcolPoint - 1
	rightLimitNegative = -I[-1] + elementsFromMidcolPoint + 1
	
	del I

	allLeftOffsets = numpy.arange(leftLimitNegative, leftLimitPositive + 1)
	allRightOffsets = numpy.arange(rightLimitNegative, rightLimitPositive + 1)
	
	curLeftOffset = 0
	curRightOffset = 0

	curLeftOffsetIDX = numpy.where(allLeftOffsets == curLeftOffset)[0][0]
	curRightOffsetIDX = numpy.where(allRightOffsets == curRightOffset)[0][0]

	#print curLeftOffsetIDX
	#print curRightOffsetIDX

	# output of the objective function values as they are computed
	arcLengthsLeft = numpy.zeros((numpy.size(allLeftOffsets)))
	arcLengthsRight = numpy.zeros((numpy.size(allRightOffsets)))
	#arcLengthsLeftStartV = [None,] * numpy.size(allLeftOffsets)
	#arcLengthsRightStartV = [None,] * numpy.size(allRightOffsets)
	searchSize = 10
	
	# use this to precompute the closed mask for the Laplace method, for speed
	#askClosed = None
	maskClosed = numpy.array(BW) # for testing

	#% algorithm
	#% for each direction perform gradient ascent until maxima is reached
	#% compute means for arc lengths on the negative and positive gradients within SearchSize
	#% am I at a maxima (current objective function greater than both means and greater than all arc lengths in the neighbourhood)
	#%	yes, terminate
	#%	no
	#%		compute gradients either side by looking at SearchSize points either side
	#%		move in the direction where the mean is greater
	#%keyboard;
	
	lastFourMoves = numpy.zeros((4), dtype = numpy.int8)
	
	#LUOfA = None
	#CholeskyOfA = None
	laplaceEquationResetFactors()
	# do the left optimisation
	while True:
		#rint "Left Iteration: " + str(curLeftOffsetIDX) + ", " + str(curLeftOffset)

		if arcLengthsLeft[curLeftOffsetIDX] == 0.0:
			#pylab.subplot(1, 2, 1); CCSegUtils.showIMG(BW); CCSegUtils.plotContour(BWContour); pylab.plot(BWContour[0, leftStart + curLeftOffset], BWContour[1, leftStart + curLeftOffset], 'b*');

			arcLengthsLeft[curLeftOffsetIDX] = endpointsFindObjectiveFunction(BWContour, curLeftOffset, curRightOffset, rightStartIDX)
		#CCSegUtils.plotContour(arcLengthsLeftStartV[curLeftOffsetIDX], closed = False); pylab.title(str(arcLengthsLeft[curLeftOffsetIDX]));
		#pylab.gcf().set_size_inches((20, 10), forward = True)
		#pylab.show()

		curObjectiveValue = arcLengthsLeft[curLeftOffsetIDX]
		
		# find the mean of the objective function if we move the left offset negative
		negObjectiveMean = 0.0
		if curLeftOffsetIDX > 0:
			negIDX = numpy.arange(curLeftOffsetIDX - 1, numpy.maximum(-1, curLeftOffsetIDX - searchSize - 1), -1)
			for z in negIDX:
				if arcLengthsLeft[z] == 0.0:
					arcLengthsLeft[z] = endpointsFindObjectiveFunction(BWContour, allLeftOffsets[z], curRightOffset, rightStartIDX)
			negObjectiveMean = numpy.mean(arcLengthsLeft[(negIDX)])
		
		# find the mean of the objective function if we move the left offset negative
		posObjectiveMean = 0.0
		if curLeftOffsetIDX < numpy.size(allLeftOffsets) - 1:
			posIDX = numpy.arange(curLeftOffsetIDX + 1, numpy.minimum(numpy.size(allLeftOffsets), curLeftOffsetIDX + searchSize + 1))
			for z in posIDX:
				if arcLengthsLeft[z] == 0.0:
					#rint "doing " + str(z)
					arcLengthsLeft[z] = endpointsFindObjectiveFunction(BWContour, allLeftOffsets[z], curRightOffset, rightStartIDX)
			posObjectiveMean = numpy.mean(arcLengthsLeft[(posIDX)])
		
		#print "cur:"
		#print curObjectiveValue
		#print "neg:"
		#print negIDX
		#print arcLengthsLeft[(negIDX)]
		#print "pos:"
		#print posIDX
		#print arcLengthsLeft[(posIDX)]

		# check to see if I am at a maxima
		# check to see if my current position is greater than all the arc lengths negative from here
		if curLeftOffsetIDX > 0:
			curGreaterThanAllNeg = numpy.all(curObjectiveValue >= arcLengthsLeft[(negIDX)])
		else:
			curGreaterThanAllNeg = True
		
		if curLeftOffsetIDX <  numpy.size(allLeftOffsets) - 1:
			curGreaterThanAllPos = numpy.all(curObjectiveValue >= arcLengthsLeft[(posIDX)])
		else:
			curGreaterThanAllPos = True

		if ((curObjectiveValue > posObjectiveMean and curObjectiveValue > negObjectiveMean) or (curGreaterThanAllNeg and curGreaterThanAllPos)):
			#print __file__ + "->" + stack()[0][3] + " (line " + str(stack()[0][2]) +  "): Maxima"
			# we are at a maxima
			break;
		else:
			# move the last three moves down
			lastFourMoves[0:3] = lastFourMoves[1:]
			
			if posObjectiveMean > negObjectiveMean:
				# move positive
				curLeftOffsetIDX = curLeftOffsetIDX + 1
				lastFourMoves[3] = 1
			else:
				# move negative
				curLeftOffsetIDX = curLeftOffsetIDX - 1
				lastFourMoves[3] = -1
			
			curLeftOffset = allLeftOffsets[curLeftOffsetIDX]

			OscCheckArray = numpy.array([-1, 1, -1, 1], dtype = numpy.int8)
			# check for oscillation
			if numpy.array_equal(lastFourMoves, OscCheckArray) or numpy.array_equal(lastFourMoves, -OscCheckArray):
				print __file__ + "->" + stack()[0][3] + " (line " + str(stack()[0][2]) +  "): Oscillating"
				break;
			#print "Last four moves: " + lastFourMoves
	#rint arcLengthsLeft
	#rint curLeftOffset
	#print curLeftOffsetIDX
	#print curLeftOffset
	
	# new method for the right side
	maxSoFar = 0.0

	for z in numpy.arange(numpy.maximum(0, curRightOffsetIDX - searchSize), numpy.size(allRightOffsets)):
		if arcLengthsRight[z] == 0.0:
			#rint "doing " + str(z)
			arcLengthsRight[z] = endpointsFindObjectiveFunction(BWContour, curLeftOffset, allRightOffsets[z], rightStartIDX)
			if arcLengthsRight[z] < maxSoFar / 2.0:
				break;
			else:
				maxSoFar = numpy.maximum(maxSoFar, arcLengthsRight[z])

	
	# find the last maxima before the arc lengths start to decrease a lot
	arcLengthsRightDone = numpy.nonzero(arcLengthsRight)[0]
	
	T = arcLengthsRight[(arcLengthsRightDone)]
	dilateSize = searchSize
	if (dilateSize % 2) == 0:
		dilateSize = dilateSize + 1
	
	dilationMinusT = numpy.ravel(scipy.ndimage.morphology.grey_dilation(numpy.atleast_2d(T), footprint = numpy.ones((1, dilateSize), dtype = numpy.bool), mode = 'nearest') - T)
	I = numpy.where(dilationMinusT == 0)[0]
	curRightOffset = allRightOffsets[(arcLengthsRightDone[(I[-1])])]
	curRightOffsetIDX = numpy.nonzero(allRightOffsets == curRightOffset)[0][0]
	del T
	#print curLeftOffset
	#print curRightOffset
	
	if curLeftOffset != 0:
		leftShiftedBWContour = numpy.roll(BWContour, -curLeftOffset, axis = 1)
	else:
		leftShiftedBWContour = numpy.array(BWContour)
	
	finalContours = dict()

	finalContours['xi'] = numpy.atleast_2d(leftShiftedBWContour[0][(numpy.arange(rightStartIDX - curLeftOffset + curRightOffset))])
	finalContours['yi'] = numpy.atleast_2d(leftShiftedBWContour[1][(numpy.arange(rightStartIDX - curLeftOffset + curRightOffset))])
	finalContours['xo'] = numpy.atleast_2d(leftShiftedBWContour[0][(numpy.arange(rightStartIDX - curLeftOffset + curRightOffset, leftShiftedBWContour.shape[1]))])
	finalContours['yo'] = numpy.atleast_2d(leftShiftedBWContour[1][(numpy.arange(rightStartIDX - curLeftOffset + curRightOffset, leftShiftedBWContour.shape[1]))])
	
	if outputPNG != None:
		(head, tail) = os.path.split(outputPNG)
		try:
			os.makedirs(head)
		except OSError as exc: # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(head):
				pass
			else:
				raise Exception
		T = numpy.array(BW)
		I = numpy.nonzero(T)

		BBoxI = numpy.arange(numpy.min(I[0]) - 1, numpy.max(I[0]) + 2)
		BBoxJ = numpy.arange(numpy.min(I[1]) - 1, numpy.max(I[1]) + 2)

		T = numpy.take(numpy.take(T, BBoxI, axis = 0), BBoxJ, axis = 1)
		pylab.clf()

		pylab.subplot2grid((2, 2), (0, 0), rowspan = 1, colspan = 2)
		CCSegUtils.showIMG(T, ticks = True)
		lineProps = {'color': 'r', 'linewidth': 2}
		
		starProps = {'markersize': 20}
		pylab.plot(leftStart[0] - BBoxJ[0], leftStart[1] - BBoxI[0], 'm*', **starProps)
		pylab.plot(rightStart[0] - BBoxJ[0], rightStart[1] - BBoxI[0], 'm*', **starProps)

		pylab.plot(leftShiftedBWContour[0, 0] - BBoxJ[0], leftShiftedBWContour[1, 0] - BBoxI[0], 'g*', **starProps)
		T = rightStartIDX - curLeftOffset + curRightOffset
		pylab.plot(leftShiftedBWContour[0, T] - BBoxJ[0], leftShiftedBWContour[1, T] - BBoxI[0], 'g*', **starProps)
		
		#print BWContour.shape
		I = numpy.arange(0, BWContour.shape[1], 5)
		#print I
		for z in range(numpy.size(I)):
			pylab.text(BWContour[0][I[z]] - BBoxJ[0], BWContour[1][I[z]] - BBoxI[0], str(I[z]), color = 'm')

		#S = numpy.roll(BWContour, 1, axis = 1)
		
		#print S[:, 0:10]
		#print BWContour[:, 0:10]
		T = numpy.roll(BWContour, -1, axis = 1) - BWContour
		#G = numpy.atleast_2d(CCSegUtils.maxGaussian1D(1, Derivative = 1))
		#G = numpy.atleast_2d(numpy.array([1, 0, -1]))
		
		#T = scipy.ndimage.convolve(numpy.double(BWContour), G, mode = 'wrap')
		#print G
		#print BWContour
		#print T
		#pylab.quiver(BWContour[0] - BBoxJ[0], BWContour[1] - BBoxI[0], -T[0], T[1], color = 'b')
		#print T[:, C]
		
		pylab.quiver(BWContour[0] - BBoxJ[0], BWContour[1] - BBoxI[0], T[0], T[1], color = 'b', angles='xy')
		
		pylab.plot(BWContour[0, 0] - BBoxJ[0], BWContour[1, 0] - BBoxI[0], 'r*', **starProps)
		#pylab.plot(BWContour[0] - BBoxJ[0], BWContour[1] - BBoxI[0], 'r-')

		pylab.subplot(2, 2, 3); pylab.plot(allLeftOffsets, arcLengthsLeft, '-', linewidth = 2); pylab.plot(curLeftOffset, arcLengthsLeft[curLeftOffsetIDX], 'm*'); pylab.title(str(curLeftOffset))
		
		pylab.subplot(2, 2, 4); pylab.plot(allRightOffsets, arcLengthsRight, '-', linewidth = 2); pylab.plot(curRightOffset, arcLengthsRight[curRightOffsetIDX], 'm*'); pylab.title(str(curRightOffset))
#	CCSegUtils.plotContour(numpy.concatenate((numpy.atleast_2d(finalContours['xi']), numpy.atleast_2d(finalContours['yi'])), axis = 0), lineProps = lineProps, closed = False);
#	lineProps = {'color': 'b', 'linewidth': 2}
#	CCSegUtils.plotContour(numpy.concatenate((numpy.atleast_2d(finalContours['xo']), numpy.atleast_2d(finalContours['yo'])), axis = 0), lineProps = lineProps, closed = False);
#	pylab.plot(BWContour[0, 0], BWContour[1, 0], 'g*', markersize = 10)
#	pylab.plot(BWContour[0, rightStartIDX], BWContour[1, rightStartIDX], 'm*', markersize = 10)
		pylab.gcf().set_size_inches((20, 15), forward = True)
		#ylab.show()
		#quit()
		pylab.savefig(outputPNG)
#
	#print finalContours
	#quit()

	return finalContours


	#print leftLimitPositive
	#print leftLimitNegative
		
	#
	#CCSegUtils.showIMG(BW)
	#CCSegUtils.plotContour(BWContour)
	
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	
	# calculate extrema for each half
