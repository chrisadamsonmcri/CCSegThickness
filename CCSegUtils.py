import numpy
import pylab

import os

import scipy



#function [F, HalfMaximum] = gaussian_fwhm2d(SIGMA)

# imitates the imglob utility in fsl
# takes a image filename and checks to see if the corresponding nifti file exists, returns the prefix

supportedExtensions = ['nii', 'nii.gz']


def findNIFTIFromPrefix(filePrefix):
	for curExtension in supportedExtensions:
		if os.path.isfile(filePrefix + '.' + curExtension):
			return filePrefix + '.' + curExtension
	return None

def imglob(fileName):

	if isinstance(fileName, list):
		T = list()
		for curFileName in fileName:
			T.append(imglob(curFileName))
		return T
	else:

		filePrefix = fileName
		# strip the extension
		for curExtension in supportedExtensions:
			if fileName.endswith("." + curExtension):
				filePrefix = filePrefix[0:-(len(curExtension) + 1)]
				break
		
		for curExtension in supportedExtensions:
			if os.path.isfile(filePrefix + '.' + curExtension):
				return filePrefix
		
		return None

def gaussianFWHM2D(SIGMA):

	assert(isinstance(SIGMA, numpy.ndarray)),"SIGMA must be an array"

	DetSIGMA = SIGMA[0, 0] * SIGMA[1, 1] - SIGMA[1, 0] * SIGMA[1, 0]
	
	assert(DetSIGMA > 0),"SIGMA must be positive-definite"

	precisionMatrix = numpy.array([[SIGMA[1, 1], -SIGMA[0, 1]], [-SIGMA[1, 0], SIGMA[0, 0]]]) / DetSIGMA

	XWidth = numpy.ceil(numpy.abs(SIGMA[0, 0]) / 3);
	YWidth = numpy.ceil(numpy.abs(SIGMA[1, 1]) / 3);

	xx = numpy.arange(-XWidth, XWidth + 1)
	yy = numpy.arange(-YWidth, YWidth + 1)

	X, Y = numpy.meshgrid(xx, yy)

	XY = numpy.concatenate((numpy.atleast_2d(numpy.ravel(X)), numpy.atleast_2d(numpy.ravel(Y))), axis = 0)

	maximum = numpy.sqrt(DetSIGMA) / (2.0 * numpy.pi)
	halfMaximum = maximum / 2.0

	quadForm = -0.5 * (numpy.matrix(XY.T) * numpy.matrix(precisionMatrix))

	F = numpy.sum(numpy.array(quadForm) * numpy.array(XY.T), axis = 1)
	F = numpy.reshape(F, (numpy.size(yy), numpy.size(xx)))
	F = numpy.exp(F) * maximum

	return (F, halfMaximum)
#YWidth = ceil(abs(SIGMA(2, 2)) / 3);
#xx = -XWidth:XWidth;
#yy = -YWidth:YWidth;
#
#[X, Y] = meshgrid(xx, yy);
#
#XY = [X(:)'; Y(:)'];
#clear X Y;
#
#Maximum = 1 ./ ((2 * pi) .* sqrt(DetSIGMA));
#HalfMaximum = Maximum / 2;
#QuadForm = -0.5 * (XY' * PrecisionMatrix);
#
#F = sum(QuadForm .* XY', 2);
#F = reshape(F, length(yy), length(xx));
#F = exp(F) .* Maximum;

def parcellationStats(labelImage, pixelArea, arcLengths, arcLengthLabels):
	
	assert(numpy.size(arcLengthLabels) == numpy.size(arcLengths)),"arcLengths and labels must be the same size"
	numLabels = numpy.max(arcLengthLabels)
	
	statsToCompute = ['min', 'max', 'mean', 'median', 'std', 'var', 'area']

	STATS = dict()

	for z in range(len(statsToCompute)):
		STATS[statsToCompute[z]] = numpy.zeros(numLabels)
	
	for z in range(1, numLabels + 1):
		I = numpy.where(arcLengthLabels == z)
		if numpy.size(I) > 0:
			STATS['min'][z - 1] = numpy.min(arcLengths[I])
			STATS['max'][z - 1] = numpy.max(arcLengths[I])
			STATS['mean'][z - 1] = numpy.mean(arcLengths[I])
			STATS['median'][z - 1] = numpy.median(arcLengths[I])
			STATS['std'][z - 1] = numpy.std(arcLengths[I])
			STATS['var'][z - 1] = numpy.var(arcLengths[I])
		STATS['area'][z - 1] = numpy.count_nonzero(labelImage == z) * pixelArea
	return STATS

def plotStreamlines(C, lineProps = None):
	for z in range(len(C)):
		plotContour(C[z], lineProps = lineProps, closed = False)

def plotContour(C, lineProps = None, closed = True):
	
	if closed == True:
		AX = numpy.concatenate([C[0], numpy.array([C[0, 0]])])
		AY = numpy.concatenate([C[1], numpy.array([C[1, 0]])])
	else:
		AX = numpy.array(C[0])
		AY = numpy.array(C[1])

	if lineProps != None:
		pylab.plot(AX, AY, **lineProps)
	else:
		pylab.plot(AX, AY)

def normPDF(X, MU, SIGMA):
	XC = X - MU

	return numpy.exp(-XC * XC / SIGMA / SIGMA / 2.0) / numpy.sqrt(2.0 * numpy.pi) / SIGMA

def empty2DList(SZ):
	#assert(isinstance(SZ, tuple)),"SZ must be a tuple"
	assert(len(SZ) == 2),"SZ must have 2 elements"

	I = list()
	for z in range(SZ[0]):
		I.append([None,] * SZ[1])
	return I

def cropAutoWhitePNG(IMGFileName, padding = (10, 10)):
	assert(os.path.isfile(IMGFileName)),"PNG file does not exist"

	IMG = pylab.imread(IMGFileName)
	
	if IMG.shape[2] == 4:
		T = numpy.squeeze(numpy.take(IMG, [0, 1, 2], axis = 2))
		T = numpy.any(T < 1, axis = 2)
	elif IMG.shape[2] == 3:
		T = numpy.any(T < 1, axis = 2)
	elif IMG.ndim < 3:
		T = (T < 1)

	I = numpy.where(T)
	
	croppedIMG = numpy.array(IMG)

	for z in range(2):
		croppedIMG = numpy.take(croppedIMG, numpy.arange(numpy.min(I[z]), numpy.max(I[z]) + 1), axis = z)
	
	cornerPadding = numpy.ones((padding[0], padding[1], croppedIMG.shape[2]))
	topBottomPadding = numpy.ones((padding[0], croppedIMG.shape[1], croppedIMG.shape[2]))
	leftRightPadding = numpy.ones((croppedIMG.shape[0], padding[1], croppedIMG.shape[2]))

	T = numpy.concatenate((
	numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1),
	numpy.concatenate((leftRightPadding, croppedIMG, leftRightPadding), axis = 1),
	numpy.concatenate((cornerPadding, topBottomPadding, cornerPadding), axis = 1)), axis = 0)

	scipy.misc.imsave(IMGFileName, T)

	#pylab.clf()

	#pylab.imshow(T)
	#pylab.gcf().set_size_inches((20, 10), forward = True)
	#pylab.show()
	
	#pass	

def showIMG(IMG, extent = None, ticks = False):
	pylab.imshow(IMG, origin = 'lower', extent = extent)

	pylab.set_cmap(pylab.cm.gray)
	
	if extent == None:
		pylab.ylim((0, IMG.shape[0] - 1))
		pylab.xlim((0, IMG.shape[1] - 1))
	else:
		pylab.ylim((extent[2], extent[3]))
		pylab.xlim((extent[0], extent[1]))
	
	if not ticks:
		pylab.gca().get_xaxis().set_ticks([])
		pylab.gca().get_yaxis().set_ticks([])
	pylab.gca().invert_yaxis()

def pylabShow():
	pylab.gcf().set_size_inches((20, 10), forward = True)
	pylab.show()

# if we are using the MNI template, 
# this is done by detecting three asterisks
# if there is no asterisks, just return the file name without modification

def MNI152FLIRTTemplate():
	scriptPath = os.path.realpath(__file__)
	(head, tail) = os.path.split(scriptPath)
	
	return os.path.join(head, 'data', 'MNI152_T1_1mm_centered.nii.gz')

def interp2q(xx, yy, V, xi, yi, interpmethod = 'linear', extrapval = numpy.nan):
	
	assert(numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0]),"numpy.size(xx) == V.shape[1] and numpy.size(yy) == V.shape[0]"
	assert(numpy.array_equal(xi.shape, yi.shape)),"(numpy.array_equal(xi.shape, yi.shape)"
#	print "xx = "
#	print xx
#	print "yy = "
#	print yy
#	print "xi = "
#	print xi
#	print "yi = "
#	print yi
#
	outV = numpy.zeros(xi.shape)

	outOfMask = numpy.logical_or(numpy.logical_or(numpy.logical_or(xi < xx[0], xi > xx[-1]), yi < yy[0]), yi > yy[-1])
	
	#outV[numpy.where(outOfMask)] = extrapval
	if numpy.any(outOfMask):
		outV[outOfMask] = extrapval

	#inMaskIDX = numpy.where(numpy.logical_not(outOfMask))
	inMaskIDX = numpy.logical_not(outOfMask)

	#if numpy.size(inMaskIDX) > 0:
	if numpy.any(inMaskIDX):

		#xSpacing = xx[1] - xx[0]
		#ySpacing = yy[1] - yy[0]
		XI = (xi[inMaskIDX] - xx[0]) / (xx[1] - xx[0])
		YI = (yi[inMaskIDX] - yy[0]) / (yy[1] - yy[0])

		if interpmethod == 'nearest':
			XI = numpy.int32(numpy.round(XI))
			YI = numpy.int32(numpy.round(YI))
			outV[inMaskIDX] = V[(YI, XI)]
		elif interpmethod == 'linear':

			XFrac = XI - numpy.floor(XI)
			#I = numpy.where(XI == xx.size - 1)
			I = (XI == xx.size - 1)
			
			#if not numpy.size(I) == 0:
			if numpy.any(I):
				XFrac[I] = 1.0
				XI[I] = xx.size - 2
				
			YFrac = YI - numpy.floor(YI)
			#I = numpy.where(YI == yy.size - 1)
			I = (YI == yy.size - 1)
			
			#if not numpy.size(I) == 0:
			if numpy.any(I):
				YFrac[I] = 1.0
				YI[I] = yy.size - 2
			
			XI = numpy.int32(XI)
			YI = numpy.int32(YI)
			outV[inMaskIDX] =	(1 - YFrac) * (1 - XFrac) * V[(YI    , XI    )] + \
								(1 - YFrac) * (    XFrac) * V[(YI    , XI + 1)] + \
								(    YFrac) * (1 - XFrac) * V[(YI + 1, XI    )] + \
								(    YFrac) * (    XFrac) * V[(YI + 1, XI + 1)]
	return outV

def maxGaussian1D(SIGMA, Derivative = 0):
	GaussianDieOff = 0.0001;

	SIGMASQ = SIGMA * SIGMA;

	W = numpy.arange(1, 501)

	FirstGaussian = numpy.exp(-(W * W) / (2 * SIGMASQ))

	#MaxWidth = find(FirstGaussian > GaussianDieOff, 1, 'last');
	MaxWidth = numpy.where(FirstGaussian > GaussianDieOff);
	
	if(numpy.size(MaxWidth) == 0):
		MaxWidth = 1;
	else:
		MaxWidth = numpy.size(MaxWidth)
	
	X = numpy.arange(-MaxWidth, MaxWidth + 1)
	Exponential = numpy.exp(-(X * X) / (2 * SIGMASQ))
	
	if Derivative == 0:
		return Exponential / (numpy.sqrt(2 * numpy.pi) * SIGMA)
	elif Derivative == 1:
		return Exponential * (-X / (numpy.sqrt(2 * numpy.pi) * SIGMA * SIGMASQ))
	
