import numpy

def otsu2(IMG):
	
	hist, bin_edges = numpy.histogram(IMG, bins=range(257), range=None, normed=False, weights=None, density=True)
	CumSumHist = numpy.cumsum(hist)
	
	MaxSigmaB = 0
	for I in range(256):
		for J in range(I + 1, 256):
			OmegaZero = CumSumHist[I]
			OmegaOne = CumSumHist[J] - CumSumHist[I]
			OmegaTwo = CumSumHist[-1] - CumSumHist[J]

			if OmegaZero > 0 and OmegaOne > 0 and OmegaTwo > 0:
				MuZero = (1 / OmegaZero) * numpy.sum(numpy.arange(0, I + 1, 1) * hist[0:I + 1])
				SigmaZero = (1 / OmegaZero) * numpy.sum((numpy.arange(0, I + 1, 1) - MuZero) * (numpy.arange(0, I + 1, 1) - MuZero) * hist[0:I + 1])

				MuOne = (1 / OmegaOne) * numpy.sum(numpy.arange(I, J + 1, 1) * hist[I:J + 1])
				SigmaOne = (1 / OmegaOne) * numpy.sum((numpy.arange(I, J + 1, 1) - MuOne) * (numpy.arange(I, J + 1, 1) - MuOne) * hist[I:J + 1])
				
				MuTwo = (1 / OmegaTwo) * numpy.sum(numpy.arange(J, 256, 1) * hist[J:])
				SigmaTwo = (1 / OmegaTwo) * numpy.sum((numpy.arange(J, 256, 1) - MuTwo) * (numpy.arange(J, 256, 1) - MuTwo) * hist[J:])

				SigmaB = OmegaZero * OmegaOne * OmegaTwo * (MuZero - MuOne - MuTwo) * (MuZero - MuOne - MuTwo)
				
				if SigmaB > MaxSigmaB:
					MaxSigmaB = SigmaB
					T = numpy.array([I, J])
	return T


