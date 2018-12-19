#!/usr/bin/python

import numpy
import cv2

A = numpy.zeros((20, 20))

A[0:4, 0:3] = 1
A[10:14, 0:3] = 2

contours, heirarchy = CCSegUtils.findContours(numpy.uint8(A), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
print A
