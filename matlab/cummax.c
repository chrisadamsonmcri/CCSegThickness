#include "mex.h"

/* Input Arguments */

#define	IN_A	prhs[0]

/* Output Arguments */

#define	OUT_A	plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[])
{
	/* Input Arguments */
	const mxArray* InA;

	/* Output arguments */
	mxArray* OutA; double* OutAptr;
	
	int i;
	int NumElements;

	if (nrhs != 1) { mexErrMsgTxt("requires 1 input arguments: [B] = cummax(A)"); }
	else if (nlhs != 1) { mexErrMsgTxt("requires 1 output arguments: [B] = cummax(A)"); }
	
	InA = IN_A;
	/* Then we check the types of the arguments */

	if(!mxIsDouble(InA) || mxIsComplex(InA)) { mexErrMsgTxt("requires InA to be of type double and non-complex"); }
	
	OutA = mxDuplicateArray(InA);
	OutAptr = mxGetPr(OutA);
	NumElements = mxGetNumberOfElements(OutA);
	for(i = 1; i < NumElements; i++)
	{
		if(OutAptr[i] < OutAptr[i - 1])
		{
			OutAptr[i] = OutAptr[i - 1];
		}
	}

	OUT_A = OutA;
	return;
}
