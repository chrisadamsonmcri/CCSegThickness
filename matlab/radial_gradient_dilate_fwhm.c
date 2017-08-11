#include "mex.h"
#include "mex_utils.h"

#include <math.h>
/* Input Arguments */

#define	IN_BW					prhs[0]
#define	IN_ANGLES				prhs[1]
#define	IN_TRADITIONALDILATION	prhs[2]

/*#define	IN_		prhs[1]
#define	IN_VARIABLEI		prhs[2]
#define IN_CLASSIFICATION	prhs[3]
#define	IN_X				prhs[4]*/

/* Output Arguments */

#define	Y_OUT				plhs[0]

#define MAX(A, B)	(A > B ? A : B)
#define MIN(A, B)	(A < B ? A : B)

#define BWPIXEL(I, J) (InBWptr[(I) + NumRowsBW * (J)])
#define AnglesPIXEL(I, J) (InAngleptr[(I) + NumRowsAngles * (J)])
#define CurSTRELPIXEL(I, J) (CurSTRELptr[(I) + NumRowsCurSTREL * (J)])

mxArray* gaussian_fwhm2d_strel_mask(const mxArray* SIGMA)
{
	/* mex Version of gaussian_fwhm2d, does not do input checking
	 returns a mask which is 1 if the element is larger than half the maximum */
	double* SIGMAptr;

	double DetSIGMA;
	mxArray* PrecisionMatrix; double* PrecisionMatrixptr;
	
	mxArray* xx; double* xxptr;
	mxArray* yy; double* yyptr;
	
	int XWidth;
	int YWidth;
	
	/* Iteration variable */
	int z;
	
	int CurI, CurJ;

	double Maximum;
	double CurF;

	/* For cropping */
	int MinI, MaxI, MinJ, MaxJ;
	
	int NumRowsCroppedF, NumColsCroppedF;
	int NumRowsF, NumColsF;

	/* Output */
	mxArray* F; mxLogical* Fptr;
	mxArray* CroppedF; mxLogical* CroppedFptr;

	SIGMAptr = mxGetPr(SIGMA);
	
	DetSIGMA = SIGMAptr[0] * SIGMAptr[3] - SIGMAptr[1] * SIGMAptr[2];
	
	PrecisionMatrix = mxCreateDoubleMatrix(2, 2, mxREAL);
	PrecisionMatrixptr = mxGetPr(PrecisionMatrix);

	PrecisionMatrixptr[0] = SIGMAptr[3] / DetSIGMA;
	PrecisionMatrixptr[1] = -SIGMAptr[1] / DetSIGMA;
	PrecisionMatrixptr[2] = -SIGMAptr[2] / DetSIGMA;
	PrecisionMatrixptr[3] = SIGMAptr[0] / DetSIGMA;
	
	XWidth = (int)ceil(fabs(SIGMAptr[0]) / 3);
	YWidth = (int)ceil(fabs(SIGMAptr[3]) / 3);
	
	NumRowsF = (YWidth << 1) + 1;
	NumColsF = (XWidth << 1) + 1;
	xx = mxCreateDoubleMatrix(1, NumRowsF, mxREAL);
	yy = mxCreateDoubleMatrix(1, NumColsF, mxREAL);
	xxptr = mxGetPr(xx);
	yyptr = mxGetPr(yy);
	
	for(z = 0; z < NumColsF; z++)
	{
		xxptr[z] = -XWidth + z;
	}
	
	for(z = 0; z < NumRowsF; z++)
	{
		yyptr[z] = -YWidth + z;
	}

	Maximum = 1.0 / (2.0 * M_PI * sqrt(DetSIGMA));
	/*mexPrintf("Maximum: %f\n", Maximum);
	mexPrintf("XWidth: %f\n", Maximum);
	mexPrintf("Maximum: %f\n", Maximum);*/

	F = mxCreateLogicalMatrix(NumRowsF, NumColsF);
	Fptr = mxGetLogicals(F);

	MaxI = -1;
	MaxJ = -1;
	MinI = NumRowsF;
	MinJ = NumColsF;

	z = 0;
	for(CurJ = 0; CurJ < NumColsF; CurJ++)
	{
		for(CurI = 0; CurI < NumRowsF; CurI++)
		{
			CurF = Maximum * exp(-0.5 * (
						(xxptr[CurJ] * PrecisionMatrixptr[0] + yyptr[CurI] * PrecisionMatrixptr[1]) * xxptr[CurJ] +
						(xxptr[CurJ] * PrecisionMatrixptr[2] + yyptr[CurI] * PrecisionMatrixptr[3]) * yyptr[CurI]));
			/*
			   [a b]	[c d  [a]
						 d e] [b]
			   = (a*c + b*d) * a + (a*d + b*e) * b
		   	*/
			Fptr[z] = (CurF > Maximum / 2.0);
			if(Fptr[z])
			{
				if(CurI < MinI) MinI = CurI;
				if(CurI > MaxI) MaxI = CurI;
				if(CurJ < MinJ) MinJ = CurJ;
				if(CurJ > MaxJ) MaxJ = CurJ;
			}
			z++;
		}
	}
	mxDestroyArray(xx);
	mxDestroyArray(yy);
	mxDestroyArray(PrecisionMatrix);
	NumRowsCroppedF = MaxI - MinI + 1;
	NumColsCroppedF = MaxJ - MinJ + 1;
	mexPrintf("(%d, %d)\n", NumRowsCroppedF, NumColsCroppedF);
	/*if(MaxI != -1 && MaxJ != -1)
	{
		CroppedF = mxCreateLogicalMatrix(NumRowsCroppedF, NumColsCroppedF);
		CroppedFptr = mxGetLogicals(CroppedF);
		
		z = 0;
		for(CurJ = 0; CurJ < NumColsCroppedF; CurJ++)
		{
			for(CurI = 0; CurI < NumRowsCroppedF; CurI++)
			{
				CroppedFptr[z] = Fptr[CurI + MinI + NumRowsF * (CurJ + MinJ)];
				z++;
			}
		}
	}
	else
	{
		CroppedF = mxCreateDoubleMatrix(1, 1, mxREAL);
	}*/
	/*mxDestroyArray(F);	
	return CroppedF;*/
	return F;
}

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[])
{
	/* Input Arguments */
	const mxArray* InBW; mxLogical* InBWptr;
	const mxArray* InAngles; double* InAnglesptr;
	mxLogical TraditionalDilation;

	/* Output arguments */
	mxArray* D; mxLogical* Dptr;
	mxArray* STRELCell;

	mxArray* CurSIGMA; double* CurSIGMAptr;
	mxArray* CurSTREL; mxLogical* CurSTRELptr;
	
	/* Dimension variables */
	int NumRowsBW, NumColsBW;
	int NumRowsAngles, NumColsAngles;
	int NumRowsCurSTREL, NumColsCurSTREL;
	
	/* Creating SIGMA */
	double SQ;
	double AngleWeighting;
	
	/* For the correlation */
	int CurSTRELBWRow, CurSTRELBWCol; 
	int CurSTRELIDX;
	int CurSTRELRow, CurSTRELCol;
	int CurSTRELCentreRow, CurSTRELCentreCol;

	/* Iterates through elements of the main image */
	int z;
	int CurRow, CurCol;
	
	mxLogical Debug;

	Debug = 0;
	if (nrhs != 3) { error_and_die(__FILE__, __LINE__, "requires 1 input arguments: [F] = radial_gradient_dilate_fwhm(SIGMA)"); }
	else if (nlhs != 2) { error_and_die(__FILE__, __LINE__, "requires 1 output arguments: [F] = radial_gradient_dilate_fwhm(SIGMA)"); }
	
	InBW = IN_BW;
	InAngles = IN_ANGLES;
	TraditionalDilation = *mxGetLogicals(IN_TRADITIONALDILATION);

	/* Then we check the types of the arguments */

	if(!mxIsLogical(InBW)) { error_and_die(__FILE__, __LINE__, "requires InBW to be of type logical"); }
	if(!mxIsDouble(InAngles) || mxIsComplex(InAngles)) { error_and_die(__FILE__, __LINE__, "requires InAngles to be of type double and non-complex"); }
	
	InBWptr = mxGetLogicals(InBW);
	InAnglesptr = mxGetPr(InAngles);

	NumRowsBW = mxGetM(InBW);
	NumColsBW = mxGetN(InBW);
	
	NumRowsAngles = mxGetM(InAngles);
	NumColsAngles = mxGetN(InAngles);

	if(NumRowsAngles != NumRowsBW || NumColsAngles != NumColsBW) { error_and_die(__FILE__, __LINE__, "requires InBW to have the same dinensions as InAngles"); }
	
	CurSIGMA = mxCreateDoubleMatrix(2, 2, mxREAL);
	CurSIGMAptr = mxGetPr(CurSIGMA);

	D = mxCreateLogicalMatrix(NumRowsBW, NumColsBW);
	Dptr = mxGetLogicals(D);
	
	STRELCell = mxCreateCellMatrix(NumRowsBW, NumColsBW);

	z = 0;
	for(CurCol = 0; CurCol < NumColsBW; CurCol++)
	{
		for(CurRow = 0; CurRow < NumRowsBW; CurRow++)
		{
			/*compute strel mask*/
			/* MATLAB pseudocode for generating the SIGMA */
			/*
			R = [25, 25] + 25 * abs([cos(Angle), sin(Angle)]);

			%[GX(I, J), GY(I, J)]

			AngleWeighting = 0.9 * cos(2 * (Angle + 45 * pi / 180));

			SQ = sqrt(R(1) * R(2)) * AngleWeighting;

			SIGMA = [R(1), SQ; SQ, R(2)];*/
			if(InBWptr[z] || TraditionalDilation)
			{
				CurSIGMAptr[0] = 200 + 200 * fabs(cos(InAnglesptr[z]));
				CurSIGMAptr[3] = 200 + 200 * fabs(sin(InAnglesptr[z]));

				AngleWeighting = 0.9 * cos(2.0 * (InAnglesptr[z] + 45.0 * M_PI / 180.0));
				SQ = -sqrt(CurSIGMAptr[0] * CurSIGMAptr[3]) * AngleWeighting;
				CurSIGMAptr[1] = SQ;
				CurSIGMAptr[2] = SQ;

				CurSTREL = gaussian_fwhm2d_strel_mask(CurSIGMA);
				if(!mxIsDouble(CurSTREL))
				{
					CurSTRELptr = mxGetLogicals(CurSTREL);
					/*do correlation*/
					
					NumRowsCurSTREL = mxGetM(CurSTREL);
					NumColsCurSTREL = mxGetN(CurSTREL);
					
					CurSTRELCentreRow = (NumRowsCurSTREL + 1) >> 1;
					CurSTRELCentreCol = (NumColsCurSTREL + 1) >> 1;

					CurSTRELIDX = 0;
					for(CurSTRELCol = 0; CurSTRELCol < NumColsCurSTREL; CurSTRELCol++)
					{
						for(CurSTRELRow = 0; CurSTRELRow < NumRowsCurSTREL; CurSTRELRow++)
						{
							if(CurSTRELptr[CurSTRELIDX])
							{
								CurSTRELBWRow = CurRow + CurSTRELRow - CurSTRELCentreRow;
								CurSTRELBWRow = MIN(CurSTRELBWRow, NumRowsBW - 1);
								CurSTRELBWRow = MAX(CurSTRELBWRow, 0);
								CurSTRELBWCol = CurCol + CurSTRELCol - CurSTRELCentreCol;
								CurSTRELBWCol = MIN(CurSTRELBWCol, NumColsBW - 1);
								CurSTRELBWCol = MAX(CurSTRELBWCol, 0);
								
								if(Debug)
								{
									if(CurRow == 0 && CurCol == 0)
									{
										mexPrintf("(%d, %d), (%d, %d), (%d, %d)\n", CurRow, CurCol, CurRow + CurSTRELRow - CurSTRELCentreRow, CurCol + CurSTRELCol - CurSTRELCentreCol, CurSTRELBWRow, CurSTRELBWCol);
									}
								}
								/* Traditional dilation */
								if(TraditionalDilation)
								{
									if(BWPIXEL(CurSTRELBWRow, CurSTRELBWCol))
									{
										Dptr[z] = 1;
									}
								}
								else
								{
									Dptr[CurSTRELBWRow + NumRowsBW * CurSTRELBWCol] = 1;
								}
							}
							CurSTRELIDX++;
						}
					}
					/*mxSetCell(STRELCell, z, CurSTREL);*/
				}
			/*mxDestroyArray(CurSTREL);*/
			}
			z++;
		}
	}	
	plhs[0] = D;/*gaussian_fwhm2d_strel_mask(SIGMA);*/
	plhs[1] = STRELCell;/*_fwhm2d_strel_mask(SIGMA);*/

	mxDestroyArray(CurSIGMA);
	return;
}
