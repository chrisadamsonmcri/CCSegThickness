import numpy
cimport numpy

def streamline2D(numpy.ndarray xx, numpy.ndarray yy, numpy.ndarray FX, numpy.ndarray FY, double sx, double sy, double step, int maxVert, numpy.ndarray contourX, numpy.ndarray contourY):
	#print xx.shape[0]
	#print xx.shape[1]
	
	cdef int numVerts = 0

	# these are the current x and y positions
	cdef double x = sx
	cdef double y = sy
	# these are the indices, the left cell of the current x and y positions
	cdef int xi
	cdef int yi
	
	cdef int xdim = numpy.size(xx)
	cdef int ydim = numpy.size(yy)
	# x and y coordinates of the cell that x, y lies in, x0, y0 are the left upper, x1, y1 are the right lower
	cdef double x0, y0, x1, y1

	# current coordinate point index
	#def int i

	cdef double ui, vi, dx, dy, imax

	# linear interpolant weights
	cdef double a, b, c, d
	
	# contourXmin[I] = min(contourX[I], contourX[I + 1])
	# contourXmin[I] = min(contourX[I], contourX[I + 1])
	# contourYmax[I] = max(contourY[I], contourY[I + 1])
	# contourYmax[I] = max(contourY[I], contourY[I + 1])

	cdef numpy.ndarray contourXmin
	cdef numpy.ndarray contourXmax
	cdef numpy.ndarray contourYmin
	cdef numpy.ndarray contourYmax

	cdef numpy.ndarray outContour = numpy.zeros([2, maxVert], dtype = numpy.double)
	# intersection test variables

	cdef double LineX1, LineX2, LineX3, LineX4
	cdef double LineY1, LineY2, LineY3, LineY4
	cdef double uAnum, uBnum, uden
	cdef int i
	cdef int intersected

	contourXmin = numpy.minimum(contourX[:-1], contourX[1:])
	contourXmax = numpy.maximum(contourX[:-1], contourX[1:])
	contourYmin = numpy.minimum(contourY[:-1], contourY[1:])
	contourYmax = numpy.maximum(contourY[:-1], contourY[1:])

	while True:
		if x < 0 or x > xdim - 1 or y < 0 or y > ydim - 1 or numVerts >= maxVert:
			break;
		
		#print "(x, y): " + str(x) + " " + str(y)
		ix = int(x)
		iy = int(y)
		
		if ix == xdim - 1:
			ix = ix - 1
		if iy == ydim - 1:
			iy = iy - 1
		#print "(ix, iy): " + str(ix) + " " + str(iy)

		xfrac = x - float(ix)
		yfrac = y - float(iy)

		#print "(xfrac, yfrac): " + str(xfrac) + " " + str(yfrac)

		a = (1 - xfrac) * (1 - yfrac)
		b = (    xfrac) * (1 - yfrac)
		c = (1 - xfrac) * (    yfrac)
		d = (    xfrac) * (    yfrac)
		
		x0 = xx[ix]; x1 = xx[ix + 1]
		y0 = yy[iy]; y1 = yy[iy + 1]
		#print "(x0, x1): " + str(x0) + " " + str(x1)
		#print "(y0, y1): " + str(y0) + " " + str(y1)
		
		# this is the line-line intersect section where we look for intersections with 
		outContour[0, numVerts] = x0 * (1.0 - xfrac) + x1 * xfrac
		outContour[1, numVerts] = y0 * (1.0 - yfrac) + y1 * yfrac
		
		intersected = 0
		if numVerts >= 1:
			oldx = outContour[0, numVerts - 1]
			newx = outContour[0, numVerts    ]
			oldy = outContour[1, numVerts - 1]
			newy = outContour[1, numVerts    ]
		
			for i in range(numpy.size(contourXmin)):
				if not (\
					(oldx <= contourXmin[i] and newx <= contourXmin[i]) or \
					(oldx >= contourXmax[i] and newx >= contourXmax[i]) or \
					(oldy <= contourYmin[i] and newy <= contourYmin[i]) or \
					(oldy >= contourYmax[i] and newy >= contourYmax[i])):
						#IntersectionTestsptr[CurStartPoint]++;
						
						#print "contour: " + str(i)
						LineX1 = oldx;
						LineX2 = newx;
						LineX3 = contourX[i];
						LineX4 = contourX[i + 1];
						
						LineY1 = oldy;
						LineY2 = newy;
						LineY3 = contourY[i];
						LineY4 = contourY[i + 1];
						
						uAnum = (LineX4 - LineX3) * (LineY1 - LineY3) - (LineY4 - LineY3) * (LineX1 - LineX3);
						uBnum = (LineX2 - LineX1) * (LineY1 - LineY3) - (LineY2 - LineY1) * (LineX1 - LineX3);

						uden = (LineY4 - LineY3) * (LineX2 - LineX1) -  (LineX4 - LineX3) * (LineY2 - LineY1);
						#print uden
						if (uden != 0):
							uA = uAnum / uden
							uB = uBnum / uden
							
							xIntersection = oldx + uA * (newx - oldx)
							yIntersection = oldy + uA * (newy - oldy)
							
							#print str(uA) + " " + str(uB)
							if (uA >= 0.0 and uA <= 1.0 and uB >= 0.0 and uB <= 1.0):
								outContour[0, numVerts] = xIntersection
								outContour[1, numVerts] = yIntersection
								#Intersectedptr[CurStartPoint] = true;
								#LineNumIntersectedptr[CurStartPoint] = i;
								#uBValueptr[CurStartPoint] = uB;
								intersected = 1
								numVerts = numVerts + 1
								break
							#if uA >= 0.0 and uA <= 1.0 and uB >= 0.0 and uB <= 1.0):
						#if (uden != 0):
				#if not (\
			#for i in range(numpy.size(contourXmin)):
		#if numVerts >= 1:

		if intersected == 1:
			break;

		if numVerts >= 2:
			if outContour[0, numVerts] == outContour[0, numVerts - 1] and outContour[1, numVerts] == outContour[1, numVerts - 1]:
				break;

		numVerts = numVerts + 1

		ui = FX[iy, ix] * a + FX[iy, ix + 1] * b + FX[iy + 1, ix] * c + FX[iy + 1, ix + 1] * d
		vi = FY[iy, ix] * a + FY[iy, ix + 1] * b + FY[iy + 1, ix] * c + FY[iy + 1, ix + 1] * d
		#print "1 (ui, vi): " + str(ui) + " " + str(vi)
		
		dx = x1 - x0
		dy = y1 - y0
		
		#print "(dx, dy): " + str(dx) + " " + str(dy)

		if dx > 0:
			ui = ui / dx
		if dy > 0:
			vi = vi / dy
		#print "2 (ui, vi): " + str(ui) + " " + str(vi)
		
		if numpy.abs(ui) > numpy.abs(vi):
			imax = numpy.abs(ui)
		else:
			imax = numpy.abs(vi)

		if imax == 0:
			break;

		imax = step / imax

		ui = ui * imax
		vi = vi * imax
		#print "3 (ui, vi): " + str(ui) + " " + str(vi)

		x = x + ui
		y = y + vi

	return (outContour[:, 0:numVerts], intersected)
