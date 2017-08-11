/* Copyright 1984-2002 The MathWorks, Inc.  */

/*

  stream3c.c   2D and 3D streamline MEX file
  The calling syntax is:

      verts = stream3c( x,y,z,    u,v,w,  sx,sy,sz, step, maxvert)
      verts = stream3c( [],[],[], u,v,w,  sx,sy,sz, step, maxvert)
      verts = stream3c( x,y,[],   u,v,[], sx,sy,[], step, maxvert)
      verts = stream3c( [],[],[], u,v,[], sx,sy,[], step, maxvert)

*/

static char rcsid[] = "$Revision: 1.8 $";

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "mex.h"
 

/* Input Arguments */


#define	X_IN	    prhs[0]
#define	Y_IN	    prhs[1]
#define	THETA_IN	prhs[2]
#define	SX_IN	    prhs[3]
#define	SY_IN	    prhs[4]
#define	STEP_IN	    prhs[5]
#define	MAXVERT_IN  prhs[6]
#define	WITHORAGAINST_IN  prhs[7]

#ifndef NULL
#define NULL 0
#endif

/* Output Arguments */

#define	VERTS_OUT plhs[0]


int vertLimit;
double *verts;
int allocIncrement=2000;

/************************* allocVertArray **************************/
/* 
 * This subroutine allocates memory for the vertex array
 */
void 
allocVertArray()
{
    vertLimit += allocIncrement;
    if (verts==NULL) 
      {
	  /* use malloc for the first time */
	  if ((verts = (double *) malloc(vertLimit * 3 * sizeof(double))) == NULL) 
	    {
		mexErrMsgTxt("not enough memory to store vertices\n");
	    }
      }
    else 
      {
	  double *oldVerts = verts;
	  /* use realloc from now on */
	  if ((verts = (double *) realloc(verts, vertLimit * 3 * sizeof(double))) == NULL) 
	    {
		if (oldVerts)
		  free(oldVerts);
		mexErrMsgTxt("not enough memory to store vertices\n");
	    }
      }
}
	
#define GETX2(X) xvec[(X)]
#define GETY2(Y) yvec[(Y)]

#define GETU2(X,Y) ugrid[(Y) + ydim*(X)]
#define GETV2(X,Y) vgrid[(Y) + ydim*(X)]
#define GETW2(X,Y) wgrid[(Y) + ydim*(X)]
#define GETA2(X,Y) anglegrid[(Y) + ydim*(X)]

/*
 * 2D streamline(u,v,sx,sy)
 */
int
traceStreamUV(double *anglegrid,
	      int xdim, int ydim, 
	      double sx, double sy,
	      double step, int maxVert)
{
    int numVerts = 0;
    double x = sx-1, y = sy-1;
    int ix, iy;
    double xfrac, yfrac, ui, vi;
    double a,b,c,d, imax;

    while(1)
      {
	  /* if outside of the volume, done */
          if (x<0 || x>xdim-1 ||
              y<0 || y>ydim-1 ||
	      numVerts>=maxVert)
	    break;
	  
	  ix = (int)x;
	  iy = (int)y;
	  if (ix==xdim-1) ix--;
	  if (iy==ydim-1) iy--;
	  xfrac = x-ix;
	  yfrac = y-iy;

	  /* weights for linear interpolation */
	  a=(1-xfrac)*(1-yfrac);
	  b=(  xfrac)*(1-yfrac);
	  c=(1-xfrac)*(  yfrac);
	  d=(  xfrac)*(  yfrac);
	  
	  /* allocate if necessary and build up output array */
	  if ( numVerts+1 >= vertLimit ) 
	    allocVertArray();
	  verts[2*numVerts + 0] = x+1;
	  verts[2*numVerts + 1] = y+1;

	  /* if already been here, done */
	  if (numVerts>=2)
	    if (verts[2*numVerts + 0] == verts[2*(numVerts-2) + 0] &&
		verts[2*numVerts + 1] == verts[2*(numVerts-2) + 1])
	      break;
	  
	  numVerts++;
	  
		       
	  /* interpolate the vector data */
	  ui = 
	    GETU2(ix,  iy  )*a + GETU2(ix+1,iy  )*b +
	    GETU2(ix,  iy+1)*c + GETU2(ix+1,iy+1)*d;
	  vi = 
	    GETV2(ix,  iy  )*a + GETV2(ix+1,iy  )*b +
	    GETV2(ix,  iy+1)*c + GETV2(ix+1,iy+1)*d;
	  
	  /* calculate step size, if 0, done */
	  if (fabs(ui)>fabs(vi)) imax=fabs(ui); else imax=fabs(vi);
	  if (imax==0) break;
	  
	  imax = step/imax;
	  
	  ui *= imax;
	  vi *= imax;

	  /* update the current position */
	  x += ui;
	  y += vi;
	  
      }
    
    return(numVerts);
}

/*
 * 2D streamline(x,y,u,v,sx,sy)
 */
int
traceStreamXYUV(double *xvec, double *yvec,
		double *anglegrid,
	       int xdim, int ydim, 
	       double sx, double sy,
	      double step, int maxVert)
{
    int numVerts = 0;
    double x = sx-1, y = sy-1;
    int ix, iy;
    double x0,x1,y0,y1,xi,yi,dx,dy;
    double xfrac, yfrac, ui, vi;
    double a,b,c,d, imax;

    while(1)
      {
          if (x<0 || x>xdim-1 ||
              y<0 || y>ydim-1 ||
	      numVerts>=maxVert)
	    break;
	  
	  ix = (int)x;
	  iy = (int)y;
	  if (ix==xdim-1) ix--;
	  if (iy==ydim-1) iy--;
	  xfrac = x-ix;
	  yfrac = y-iy;

	  a=(1-xfrac)*(1-yfrac);
	  b=(  xfrac)*(1-yfrac);
	  c=(1-xfrac)*(  yfrac);
	  d=(  xfrac)*(  yfrac);
	  
	  if ( numVerts+1 >= vertLimit ) 
	    allocVertArray();
	  
	  x0 = GETX2(ix); x1 = GETX2(ix+1);
	  y0 = GETY2(iy); y1 = GETY2(iy+1);
	  xi = x0*(1-xfrac) + x1*xfrac;
	  yi = y0*(1-yfrac) + y1*yfrac;
	  
	  verts[2*numVerts + 0] = xi;
	  verts[2*numVerts + 1] = yi;
	  if (numVerts>=2)
	    if (verts[2*numVerts + 0] == verts[2*(numVerts-2) + 0] &&
		verts[2*numVerts + 1] == verts[2*(numVerts-2) + 1])
	      break;
	  numVerts++;
	  
	  ui = 
	    GETU2(ix,  iy  )*a + GETU2(ix+1,iy  )*b +
	    GETU2(ix,  iy+1)*c + GETU2(ix+1,iy+1)*d;
	  vi = 
	    GETV2(ix,  iy  )*a + GETV2(ix+1,iy  )*b +
	    GETV2(ix,  iy+1)*c + GETV2(ix+1,iy+1)*d;
	  
	  dx = x1-x0;
	  dy = y1-y0;
	  if (dx) ui /= dx;
	  if (dy) vi /= dy;
	  
	  if (fabs(ui)>fabs(vi)) imax=fabs(ui); else imax=fabs(vi);
	  if (imax==0) break;
	  
	  imax = step/imax;
	  
	  ui *= imax;
	  vi *= imax;
	  
	  x += ui;
	  y += vi;
	  
      }
    
    return(numVerts);
}


void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[])
{
	double *x, *y, *angle, *sx, *sy, *step, *maxVert;
	int xSize, ySize;
	const int *dims;
	int numVerts;
	int i, NumPoints;
	
	mxLogical withOrAgainst;
	mxArray* VOutTemp; double* VOutTempptr;
	
	/* Check for proper number of arguments */

	if (nrhs != 8)
	{
		mexErrMsgTxt("mystream2c_vector requires 8 input arguments. [verts] = stream3c(x,y,a,sx,sy,step,maxVert,withOrAgainst);");
	}
	else if (nlhs != 1) 
	{
		mexErrMsgTxt("mystream2c_vector requires 1 output arguments. [verts] = stream3c(x,y,a,sx,sy,step,maxVert,withOrAgainst)");
	}

	x = mxGetPr( X_IN );
	y = mxGetPr( Y_IN );
	angle = mxGetPr( A_IN );
	/*fieldgrid = mxGetPr( FIELDGRID_IN );*/
	sx = mxGetPr( SX_IN );
	sy = mxGetPr( SY_IN );
	step = mxGetPr( STEP_IN    );
	maxVert = mxGetPr( MAXVERT_IN );
	withOrAgainst = *mxGetLogicals(WITHORAGAINST_IN);

	NumPoints = mxGetNumberOfElements(SX_IN);

	dims = mxGetDimensions(A_IN);

	VERTS_OUT = mxCreateCellMatrix(NumPoints, 1);
	
	xSize = dims[1];
	ySize = dims[0];

	if (xSize <= 1 || ySize <= 1)
	{
		mexErrMsgTxt("stream2c_vector requires that all three dimensions be greater than 1");
	}

	for (i = 0; i < NumPoints; i++)
	{
		vertLimit = numVerts = 0;
		verts = NULL;

		if(x != NULL && y != NULL)
		{
			numVerts = traceStreamXYUV(x, y, angle, xSize, ySize, sx[i], sy[i], *step, (int)(*maxVert), withOrAgainst);
		}
		else
		{
			numVerts = traceStreamUV(angle, xSize, ySize, sx[i], sy[i], *step, (int)(*maxVert), withOrAgainst);
		}
		if (verts)
		{
		
			VOutTemp = mxCreateDoubleMatrix(2, numVerts, mxREAL);
			VOutTempptr = mxGetPr(VOutTemp);
			
			/* compy the results into the output parameters*/
			memcpy((char *)VOutTempptr, (char *)verts, numVerts*2*sizeof(double));

			mxSetCell(VERTS_OUT, i, VOutTemp);
			
			/*mxDestroyArray(VOutTemp);*/

			free(verts);
		}
	}
	return;	
}


