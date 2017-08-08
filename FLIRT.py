import numpy

#function [rotmat, skew, scales, transl, angles] = fsl_decomp_aff(A)
def fsl_decomp_aff(A):
#	aff3 = A(1:3, 1:3);
	aff3 = A.take([0, 1, 2], axis=0).take([0, 1, 2], axis=1)
#	x = A(1:3, 1);
	x = A[0:3, 0]
#	y = A(1:3, 2);
	y = A[0:3, 1]
#	z = A(1:3, 3);
	z = A[0:3, 2]
#	sx = sqrt(sum(x .* x));
	sx = numpy.sqrt(numpy.sum(numpy.multiply(x, x)))
#	% sy = std::sqrt( dot(y,y) - (Sqr(dot(x,y)) / Sqr(sx)) );
#	DXY = dot(x, y);
	DXY = numpy.sum(numpy.multiply(x, y))
#	sy = sqrt(dot(y, y) - (DXY * DXY / (sx * sx)));
	sy = numpy.sqrt(numpy.sum(numpy.multiply(y, y)) - (DXY * DXY / (sx * sx)))

#	a = dot(x, y) / (sx * sy);
	a = numpy.sum(numpy.multiply(x, y)) / (sx * sy);

	x0 = x / sx;
	
	y0 = y / sy - a * x0;
#	% sz = std::sqrt(dot(z,z) - Sqr(dot(x0,z)) - Sqr(dot(y0,z)));
	Dx0z = numpy.sum(numpy.multiply(x0, z));
	Dy0z = numpy.sum(numpy.multiply(y0, z));
#	sz = sqrt(dot(z, z) - Dx0z * Dx0z - Dy0z * Dy0z);
	sz = numpy.sqrt(numpy.sum(numpy.multiply(z, z)) - Dx0z * Dx0z - Dy0z * Dy0z);
#
#	b = dot(x0, z) / sz;
	b = numpy.sum(numpy.multiply(x0, z)) / sz;
#	c = dot(y0, z) / sz;
	c = numpy.sum(numpy.multiply(y0, z)) / sz;

#	% params(7) = sx;  params(8) = sy;  params(9) = sz;
#	scales = diag([sx, sy, sz]);
	scales = numpy.matrix(numpy.diag([sx, sy, sz]))
#	% diag(scales,diagvals);
#	% Real skewvals[] = {1,a,b,0 , 0,1,c,0 , 0,0,1,0 , 0,0,0,1}; 
#	skew = [1, a, b, 0; 0, 1, c, 0; 0, 0, 1, 0; 0, 0, 0, 1];
	skew = numpy.matrix(numpy.vstack([[1, a, b, 0], [0, 1, c, 0], [0, 0, 1, 0], [0, 0, 0, 1]]))
#	% params(10) = a;  params(11) = b;  params(12) = c; 
#	% Matrix rotmat(3,3);
	rotmat = numpy.matrix(numpy.zeros((3, 3)))
#	rotmat = aff3 * inv(scales) * inv(skew(1:3, 1:3));
	rotmat = aff3 * scales.I * skew[0:3, 0:3].I

#	transl = A(1:3, 4);
	transl = A[0:3, 3]
#	print "rotmat"
#	print rotmat
#	print "transl"
#	print transl
#	% ColumnVector transl(3);
#	% transl = affmat.SubMatrix(1,3,1,3)*centre + affmat.SubMatrix(1,3,4,4)
#	% 	 - centre;
#
#	angles = zeros(1, 3);
	angles = numpy.zeros((3, 1))
#
#	% from miscmaths/miscmaths.cc
#	% construct euler angles
#	cy = sqrt(rotmat(1, 1) * rotmat(1, 1) + rotmat(1, 2) * rotmat(1, 2));
	cy = numpy.sqrt(rotmat[0, 0] * rotmat[0, 0] + rotmat[0, 1] * rotmat[0, 1])

#	if cy < 1e-4
	if cy < 1e-4:
#		cx = rotmat(2,2);
		cx = rotmat[1, 1]
#		sx = -rotmat(3,2);
		sx = -rotmat[2, 1]
#		sy = -rotmat(1,3);
		sy = -rotmat[0, 2]
#		angles(1) = atan2(sx, cx);
		angles[0] = numpy.arctan2(sx, cx)
#		angles(2) = atan2(sy, 0);
		angles[1] = numpy.arctan2(sy, 0)
#		angles(3) = 0.0;
	else:
#		cz = rotmat(1,1)/cy;
		cz = rotmat[0, 0] / cy
#		sz = rotmat(1,2)/cy;
		sz = rotmat[0, 1] / cy
#		cx = rotmat(3,3)/cy;
		cx = rotmat[2, 2] / cy
#		sx = rotmat(2,3)/cy;
		sx = rotmat[1, 2] / cy
#		sy = -rotmat(1,3);
		sy = -rotmat[0, 2]
#		angles(1) = atan2(sx,cx);
		angles[0] = numpy.arctan2(sx,cx)
#		angles(2) = atan2(sy,cy);
		angles[1] = numpy.arctan2(sy,cy)
#		angles(3) = atan2(sz,cz);
		angles[2] = numpy.arctan2(sz,cz)
#	end
#	print "angles"
#	print angles
	return (rotmat, skew, scales, transl, angles)
#end
