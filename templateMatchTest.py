#!/usr/bin/python

import numpy
import cv2

templateIMG = numpy.zeros((3, 3), dtype = numpy.single)

N = 7

IMG = numpy.reshape(numpy.arange(N * N), [N, N])
print IMG

normXCorrAVW = cv2.matchTemplate(numpy.single(IMG), numpy.single(templateIMG), cv2.TM_SQDIFF)

print normXCorrAVW

X = numpy.zeros([IMG.shape[0] - templateIMG.shape[0] + 1, IMG.shape[1] - templateIMG.shape[1] + 1])

for curCol in numpy.arange(IMG.shape[1] - templateIMG.shape[1] + 1):
	for curRow in numpy.arange(IMG.shape[0] - templateIMG.shape[0] + 1):
		T = IMG.take(numpy.arange(curRow, curRow + templateIMG.shape[0]), axis = 0).take(numpy.arange(curCol, curCol + templateIMG.shape[1]), axis = 1)
		X[curRow, curCol] = numpy.sum(T * T)

print X
