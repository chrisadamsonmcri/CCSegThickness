#!/usr/bin/env python

import numpy
import CCSegUtils

import CCSegUtils
import matplotlib.pyplot as plt
import pylab 


def rotMatrixX(theta):
	return numpy.matrix([
		[                1,                0,                 0, 0],
		[                0, numpy.cos(theta), -numpy.sin(theta), 0],
		[                0, numpy.sin(theta),  numpy.cos(theta), 0],
		[                0,                0,                 0, 1]])

def rotMatrixY(theta):
	return numpy.matrix([
		[ numpy.cos(theta),                0, -numpy.sin(theta), 0],
		[                0,                1,                 0, 0],
		[ numpy.sin(theta),                0,  numpy.cos(theta), 0],
		[                0,                0,                 0, 1]])

def rotMatrixZ(theta):
	return numpy.matrix([
		[ numpy.cos(theta), -numpy.sin(theta),                0, 0],
		[ numpy.sin(theta),  numpy.cos(theta),                0, 0],
		[                0,                 0,                1, 0],
		[                0,                 0,                0, 1]])

def translationMatrix(T):
	return numpy.matrix([
		[1, 0, 0, T[0]], 
		[0, 1, 0, T[1]], 
		[0, 0, 1, T[2]], 
		[0, 0, 0,    1]])

#@profile
def corrCoefCost(IMG):
	if IMG.shape[0] % 2 == 0:
		#print "even"
		leftRows = numpy.arange(IMG.shape[0] / 2, dtype = numpy.uint16)
		rightRows = numpy.arange(IMG.shape[0] - 1, IMG.shape[0] / 2 - 1, -1, dtype = numpy.uint16)
	else:
		#print "odd"
		leftRows = numpy.arange(numpy.floor(IMG.shape[0] / 2), dtype = numpy.uint16)
		rightRows = numpy.arange(IMG.shape[0] - 1, numpy.ceil(IMG.shape[0] / 2), -1, dtype = numpy.uint16)
	
	leftIMG = numpy.take(IMG, leftRows, axis = 0)
	rightIMG = numpy.take(IMG, rightRows, axis = 0)

	M = numpy.logical_and(leftIMG > 0, rightIMG > 0)

	if numpy.all(M == False):
		return 0
	
	leftV = leftIMG[M]
	rightV = rightIMG[M]
	#	CCSegUtils.imshow(IMG[:, :, 100])
	#	plt.show()
	#meanLeftV = numpy.mean(leftV)
	#meanRightV = numpy.mean(rightV)

	demeanLeftV = leftV - numpy.mean(leftV)
	demeanRightV = rightV - numpy.mean(rightV)

	#varLeftV = numpy.sqrt(numpy.sum(demeanLeftV * demeanLeftV))
	#varRightV = numpy.sqrt(numpy.sum(demeanRightV * demeanRightV))
	varLeftV = numpy.sqrt(numpy.dot(demeanLeftV, demeanLeftV))
	varRightV = numpy.sqrt(numpy.dot(demeanRightV, demeanRightV))

	#C = numpy.sum(demeanLeftV * demeanRightV) / varLeftV / varRightV
	C = numpy.dot(demeanLeftV, demeanRightV) / varLeftV / varRightV
	#C = numpy.corrcoef(leftV, rightV)
	#C = C[0, 1]
	#print C
	return C

def translateIMG(IMG, T):
	#print numpy.all(IMG == 0)
	if T == 0:
		outIMG = numpy.array(IMG)
	elif T < 0:
		outIMG = numpy.zeros_like(IMG)
		outIMG[-T:] = IMG[:T]
	else:
		outIMG = numpy.zeros_like(IMG)
		outIMG[:-T] = IMG[T:]

	return outIMG

#from scipy.interpolate import RegularGridInterpolator

#mport time
#@profile
def transformIMG(IMG, IMGxx, IMGyy, IMGzz, RotY, RotZ, TransX, IMGCoordsM = None, interpmethod = 'linear'):
	RY = rotMatrixY(RotY * numpy.pi / 180.0)
	RZ = rotMatrixZ(RotZ * numpy.pi / 180.0)
	R = RY * RZ

	T = translationMatrix(numpy.array([TransX, 0, 0]))
	M = numpy.single(T * R)
	# M * IMGCoords

	# remove the need to do the meshgrid by broadcasting here

	#XI = IMGCoords[0] * M[0, 0] + IMGCoords[1] * M[0, 1] + IMGCoords[2] * M[0, 2] + M[0, 3]
	#YI = IMGCoords[0] * M[1, 0] + IMGCoords[1] * M[1, 1] + IMGCoords[2] * M[1, 2]# + M[1, 3] # these are always zero
	#ZI = IMGCoords[0] * M[2, 0] + IMGCoords[1] * M[2, 1] + IMGCoords[2] * M[2, 2]# + M[2, 3] # these are always zero
	
	#start_time = time.time()

	IMGxxT = numpy.single(numpy.reshape(numpy.atleast_3d(IMGxx), (IMGxx.size,          1,          1)))
	IMGyyT = numpy.single(numpy.reshape(numpy.atleast_3d(IMGyy), (         1, IMGyy.size,          1)))
	IMGzzT = numpy.single(numpy.reshape(numpy.atleast_3d(IMGzz), (         1,          1, IMGzz.size)))
	#print IMGxxT.shape
	#print IMGyyT.shape
	#print IMGzzT.shape

	XI = IMGxxT * M[0, 0] + IMGyyT * M[0, 1] + IMGzzT * M[0, 2] + M[0, 3]
	YI = IMGxxT * M[1, 0] + IMGyyT * M[1, 1] + IMGzzT * M[1, 2]
	ZI = IMGxxT * M[2, 0] + IMGyyT * M[2, 1] + IMGzzT * M[2, 2]
	
	T = CCSegUtils.interp3q(IMGyy, IMGxx, IMGzz, IMG, numpy.single(YI), numpy.single(XI), numpy.single(ZI), extrapval = 0, interpmethod = interpmethod)
	#return CCSegUtils.interp3q(IMGyy, IMGxx, IMGzz, IMG, YI, XI, ZI, extrapval = 0, interpmethod = 'linear')
	# do slice-by-slice to conserve memory
	# your code
	#elapsed_time = time.time() - start_time

	#outIMG = list()
	#XITArray = numpy.zeros_like(IMG)
	#YITArray = numpy.zeros_like(IMG)
	#ZITArray = numpy.zeros_like(IMG)
#	start_time = time.time()
#
#	IMGxxT = numpy.reshape(numpy.atleast_3d(IMGxx), (IMGxx.size,          1, 1))
#	IMGyyT = numpy.reshape(numpy.atleast_3d(IMGyy), (         1, IMGyy.size, 1))
#	
#	#print IMGxxT.shape
#	#print IMGyyT.shape
#	#print M
#	XIXY = IMGxxT * M[0, 0] + IMGyyT * M[0, 1] + M[0, 3]
#	YIXY = IMGxxT * M[1, 0] + IMGyyT * M[1, 1]
#	ZIXY = IMGxxT * M[2, 0] + IMGyyT * M[2, 1]
#
#	#print IMGzz.shape
#	for z in range(IMG.shape[2]):
##		
#		XI = XIXY + (IMGzz[z] * M[0, 2])
#		YI = YIXY + (IMGzz[z] * M[1, 2])
#		ZI = ZIXY + (IMGzz[z] * M[2, 2])
##		#XITArray[:, :, z] = XIT
##		#YITArray[:, :, z] = YIT
#		#ZITArray[:, :, z] = ZIT
##		
#		outIMG.append(
#		CCSegUtils.interp3q2D(IMGyy, IMGxx, IMGzz, IMG, YI, XI, ZI, extrapval = 0, interpmethod = 'linear')
#		)
#	outIMG = numpy.concatenate(outIMG, axis = 2)
	#elapsed_time2 = time.time() - start_time
	#print "3d: " + str(elapsed_time)
	#print "2d: " + str(elapsed_time2)

#	if numpy.array_equal(T, outIMG) == False:
#		I = numpy.where(T != outIMG)
#		print I
#		print T[I]
#		print outIMG[I]
#		#print XITArray[I]
#		#print YITArray[I]
#		#print ZITArray[I]
#	#fn = RegularGridInterpolator((IMGxx, IMGyy, IMGzz), IMG, bounds_error = False, fill_value = 0)

	#XYZI = numpy.concatenate((numpy.atleast_2d(numpy.ravel(XI)).T, numpy.atleast_2d(numpy.ravel(YI)).T, numpy.atleast_2d(numpy.ravel(ZI)).T), axis = 1)

	#G = fn(XYZI)

	return T
	
	#IMGCoords = numpy.concatenate((
	#	numpy.atleast_2d(numpy.ravel(IMGCoords[0])),
	#	numpy.atleast_2d(numpy.ravel(IMGCoords[1])),
	#	numpy.atleast_2d(numpy.ravel(IMGCoords[2]))
	#	), axis = 0)
	
	#IMGCoords = numpy.array(numpy.matrix(M[0:3, 0:3]) * numpy.matrix(IMGCoords)) + M[0:3, 3]

#@profile	
def transformCost(IMG, IMGxx, IMGyy, IMGzz, RotY, RotZ, TransX, return_image = False):
	TIMG = transformIMG(IMG, IMGxx, IMGyy, IMGzz, RotY, RotZ, TransX)
	if return_image:
		return (corrCoefCost(TIMG), TIMG)
	else:
		return corrCoefCost(TIMG)

def getMidSagSlice(TIMG):
	midSlice = int(numpy.floor(TIMG.shape[0] / 2))

	# extract the midsagittal slice
	if (TIMG.shape[0] % 2 == 0):
		# even number of slices
		midSagAVW = (numpy.double(TIMG[midSlice - 1, :, :]) + numpy.double(TIMG[midSlice, :, :])) / 2.0
	else:
		# odd number of slices
		midSagAVW = numpy.double(TIMG[midSlice, :, :])
	return numpy.rot90(midSagAVW, 1)

def showMidSag(TIMG, midSagAVW = None):
	SR = 2
	SC = 2
	pylab.subplot(SR, SC, 1)
	CCSegUtils.imshow(TIMG[:, :, TIMG.shape[2] // 2])
	plt.plot([0, TIMG.shape[1]], [TIMG.shape[0] / 2, TIMG.shape[0] / 2], 'w-')
	pylab.subplot(SR, SC, 2)
	T = numpy.rot90(TIMG[TIMG.shape[0] // 2], 1)
	CCSegUtils.imshow(T)
	pylab.subplot(SR, SC, 3)
	T = numpy.rot90(TIMG[:, TIMG.shape[1] // 2], 1)
	CCSegUtils.imshow(T)
	plt.plot([T.shape[1] / 2, T.shape[1] / 2], [0, TIMG.shape[0]], 'w-')
	pylab.subplot(SR, SC, 4)
	if midSagAVW is None:
		CCSegUtils.imshow(getMidSagSlice(TIMG))
	else:
		CCSegUtils.imshow(midSagAVW)
	

def paramString(P):
	if P < 0:
		prefix = 'm'
	else:
		prefix = 'p'
	return prefix + str(numpy.abs(P)).zfill(3)
