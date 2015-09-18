#!/usr/bin/python

import CCSegPipe
import numpy
import nibabel

import pylab
import CCSegUtils

import os

numpy.set_printoptions(precision = 3, formatter = {'float':lambda x: "%.3f" % x})
#A = numpy.uint8(numpy.arange(4))

inputFile = '/mnt/addo/data/cc_seg/louisee/RawT1ReorientCropped/AMY0024_T1.nii.gz'
if not os.path.isfile(inputFile):
	inputFile = '/data/addo/louisee/RawT1ReorientCropped/AMY0024_T1.nii.gz'

NII = nibabel.load(inputFile)

NIIData = numpy.array(NII.get_data())

#I = numpy.where(NIIData > 0)
X = NIIData[NIIData > 0]

X = numpy.double(X)
X = (X - numpy.min(X)) / (numpy.max(X) - numpy.min(X))
X = numpy.uint8(numpy.round(X * 255))

#T = CCSegPipe.otsu1(X)

#import cv2

#ret2, th2 = cv2.threshold(X,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
#print ret2
#print th2
THRESH1, OmegaZeros1, MuZeros1, SigmaZeros1, OmegaOnes1, MuOnes1, SigmaOnes1, OmegaTwos1, MuTwos1, SigmaTwos1, MuOneHistArray1, SigmaB1, SigmaB2 = CCSegPipe.otsu2(X, returnWorkingValues = True)

#THRESH2, OmegaZeros2, MuZeros2, SigmaZeros2, OmegaOnes2, MuOnes2, SigmaOnes2, OmegaTwos2, MuTwos2, SigmaTwos2, MuOneHistArray2, SigmaB2, OmegaMask2 = CCSegPipe.otsu2Fast(X, returnWorkingValues = True)

#print "Zeros"
#print OmegaZeros1[0:10]
#print OmegaZeros2[0:10]
#print "---"
#print MuZeros1[0:10]
#print MuZeros2[0:10]
#print "---"
#print SigmaZeros1[0:10]
#print SigmaZeros2[0:10]
#print "---"
#
#print "Twos"
#print OmegaTwos1[0, 0:10]
#print OmegaTwos2
##print "---"
##print MuTwos1[0, 0:10]
##print MuTwos2[0:10]
##print "---"
##print SigmaTwos1[0, 0:10]
##print SigmaTwos2[0:10]
#print "---"
#print "Ones"
#print OmegaOnes1[0:6, 0:6]
#print OmegaOnes2[0:6, 0:6]
#print "---"
#print MuOnes1[0:6, 0:6]
#print MuOnes2[0:6, 0:6]
#print "---"
#print SigmaOnes1[0:6, 0:6]
#print SigmaOnes2[0:6, 0:6]
##print "---"
##
#
#print OmegaTwos1[0]
#print OmegaTwos2

#rint MuOnes1
#rint MuOnes2
#print numpy.where(OmegaOnes1 == 0)

#print OmegaTwos1[0]
#print 1 - OmegaZeros1 - OmegaTwos1[0]
#print OmegaOnes1
#print OmegaOnes2

#rint SigmaB1[6, 255]
#rint SigmaB2[6, 255]


#SR = 1
#SC = 3
#pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(OmegaOnes1); pylab.colorbar();
#pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(OmegaOnes2); pylab.colorbar();
#pylab.subplot(SR, SC, 3); CCSegUtils.showIMG(OmegaOnes2 - OmegaOnes1); pylab.colorbar();
#pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(MuOnes1); pylab.colorbar();
#pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(MuOnes2); pylab.colorbar();
#pylab.subplot(SR, SC, 3); CCSegUtils.showIMG(MuOnes2 - MuOnes1); pylab.colorbar();
#pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(SigmaOnes1); pylab.colorbar();
#pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(SigmaOnes2); pylab.colorbar();
#pylab.subplot(SR, SC, 3); CCSegUtils.showIMG(MuOneHistArray2 - MuOneHistArray1); pylab.colorbar();


#SR = 1
#SC = 2
#pylab.subplot(SR, SC, 1); CCSegUtils.showIMG(SigmaB1); pylab.colorbar();
#pylab.subplot(SR, SC, 2); CCSegUtils.showIMG(SigmaB2); pylab.colorbar();

#CCSegUtils.pylabShow()

#print numpy.max(numpy.abs(MuOneHistArray2 - MuOneHistArray1))
#pylab.show()
#print SigmaZeros1
print THRESH1


