import numpy
import scipy
import scipy.sparse

def substituteSolve(A, B):
	n = A.shape[0]
	
	assert(numpy.array_equal((n, n), A.shape)),"A must be square"
	assert(len(B.shape) == 1),"B must be a vector"
	assert(numpy.size(B) == n),"The number of elements in B must equal the number of rows/columns in A"

	upperSize = numpy.size(scipy.sparse.triu(A, 1).nonzero())
	lowerSize = numpy.size(scipy.sparse.tril(A, -1).nonzero())

	if upperSize > 0 and lowerSize == 0:
		matrixShape = 'upper'
	elif upperSize == 0 and lowerSize > 0:
		matrixShape = 'lower'
	else:
		raise Exception("Matrix not triangular")
	
	Z = numpy.zeros_like(B)

	if matrixShape == 'upper':
		Z[-1] = B[-1] / A[-1, -1]
		for z in range(n - 2, -1, -1):
			Z[z] = (B[z] - A[z, (z + 1):].dot(Z[(z + 1):])) / A[z, z]
			#pass#for z in 
		#pass
	elif matrixShape == 'lower':
		Z[0] = B[0] / A[0, 0]
		for z in range(1, n):
			Z[z] = (B[z] - A[z, 0:z].dot(Z[0:z])) / A[z, z]
	return Z


def chol(A):
	
	#L = scipy.sparse.lil_matrix(scipy.sparse.tril(A))
	L = scipy.sparse.tril(A).todense()

	n = L.shape[0]
	
	for k in range(n):
		L[k, k] = numpy.sqrt(L[k, k])
		if k < (n - 1):
			if L[k, k] != 0.0:
				L[k + 1:, k] = L[k + 1:, k] / L[k, k]

			for j in range(k + 1, n):
				if L[j, k] != 0:
					#print str(j) + ' of ' + str(n)
					L[j:, j] = L[j:, j] - L[j, k] * L[j:, k]
		
	return scipy.sparse.csr_matrix(L)
