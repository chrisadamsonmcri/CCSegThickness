#!/usr/bin/python

import numpy
import scipy
import scipy.sparse

import LinAlgHelpers

n = 10

diagonals = [[-1] * (n - 1), [2] * n, [-1] * (n - 1)]

A = scipy.sparse.diags(diagonals, [-1, 0, 1])

B = numpy.concatenate((numpy.array([1]), numpy.zeros((n - 2)), numpy.array([2])))
B = numpy.double(B)
#print A.getformat()
#print A.todense()

L = LinAlgHelpers.chol(A)
#print L

#print L.todense()
Y = LinAlgHelpers.substituteSolve(L, B)
#print Y
Z = LinAlgHelpers.substituteSolve(L.T, Y)
#print Z

