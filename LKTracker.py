import numpy
import pylab
import scipy.sparse

import CCSegUtils

def jacobianAffine(IMGSize):

#	% [Jacobian] = lk_jacobian_affine(IMGSize)
#	%
#	% DESCRIPTION
#	%   Computes the Jacobian dW/dp for the affine warp
#	%   J = [ x 0 y 0 1 0;
#	%         0 x 0 y 0 1];  
#	%
#	% PARAMETERS
#	%   IMGSize [2]: the dimensions of the image [M, N]
#	%
#	% RETURNS
#	%   Jacobian [2, M * N, 6]: row 1 is dWx/dp, row 2 is dWy/dp, slice 1 is dW/dp1
#	%
#	% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
#	% See LICENSE for full license information.
#	%

#	if(numel(IMGSize) ~= 2)
#		error('IMGSize needs to have 2 elements');
#	end
	assert(numpy.size(IMGSize) == 2),"Affine Jacobian calculation, image shape is not 2D"

#	NumPixels = prod(IMGSize);
	numPixels = numpy.prod(IMGSize)

#	[X, Y] = meshgrid(1:IMGSize(2), 1:IMGSize(1));
	X, Y = numpy.meshgrid(numpy.arange(1, IMGSize[1] + 1), numpy.arange(1, IMGSize[0] + 1))

	# create a list of None values
	#	Jacobian = cell(6);
	templateJacobian = [None,] * 6

#	Jacobian{1} = cat(1, X(:)', zeros(1, NumPixels));
	#print numpy.atleast_2d(numpy.ravel(X)).shape
	#print numpy.zeros([1, numPixels]).shape
	templateJacobian[0] = numpy.concatenate((numpy.atleast_2d(numpy.ravel(X)), numpy.zeros([1, numPixels])), axis = 0)
#	Jacobian{2} = cat(1, zeros(1, NumPixels), X(:)');
	templateJacobian[1] = numpy.concatenate((numpy.zeros([1, numPixels]), numpy.atleast_2d(numpy.ravel(X))), axis = 0)
#	Jacobian{3} = cat(1, Y(:)', zeros(1, NumPixels));
	templateJacobian[2] = numpy.concatenate((numpy.atleast_2d(numpy.ravel(Y)), numpy.zeros([1, numPixels])), axis = 0)
#	Jacobian{4} = cat(1, zeros(1, NumPixels), Y(:)');
	templateJacobian[3] = numpy.concatenate((numpy.zeros([1, numPixels]), numpy.atleast_2d(numpy.ravel(Y))), axis = 0)
#	Jacobian{5} = cat(1, ones(1, NumPixels), zeros(1, NumPixels));
	templateJacobian[4] = numpy.concatenate((numpy.ones([1, numPixels]), numpy.zeros([1, numPixels])), axis = 0)
#	Jacobian{6} = cat(1, zeros(1, NumPixels), ones(1, NumPixels));
	templateJacobian[5] = numpy.concatenate((numpy.zeros([1, numPixels]), numpy.ones([1, numPixels])), axis = 0)
	#	Jacobian = cat(3, Jacobian{:});
	return templateJacobian

def jacobianAffineDisplay(jacobian, IMGSize):

	#% Displays the jacobian as an image in the current plot using imshow

	#IMG = cell(2, 6);
	# list of two empty lists
	IMG = list()
	IMG.append([])
	IMG.append([])

	for z in range(len(jacobian)):
	#for CurParameter = 1:6
		IMG[0].append(numpy.reshape(jacobian[z][0, :], IMGSize))
		IMG[1].append(numpy.reshape(jacobian[z][1, :], IMGSize))

	IMG = numpy.concatenate((numpy.concatenate(IMG[0], axis = 1), numpy.concatenate(IMG[1], axis = 1)), axis = 0)
	#IMG = cell2mat(IMG);
	CCSegUtils.showIMG(IMG)

def steepestDescentImagesAffine(FX, FY, jacobian):

#	% [SDImages] = lk_sd_images_affine(FX, FY)
#	%
#	% DESCRIPTION
#	%   Lucas-Kanade algorithm
#	%   Computes the steepest descent images using the affine transformation
#	%
#	% M, N refer to the size of the template
#	% PARAMETERS
#	%   FX [M, N]: the warped X gradient of the image dI/dx
#	%   FY [M, N]: the warped X gradient of the image dI/dy
#	%   Jacobian [2, M * N, 6]: the Jacobian
#	%
#	% RETURNS
#	%   SDImages [M * N, 6]: the steepest descent images grad I * dW/dp
#	%   column i is (dIx + dIy)dW/dp{i}
#	%
#	% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
#	% See LICENSE for full license information.
#	%
	
	IMGSize = numpy.size(FX)

	SDImages = numpy.zeros([IMGSize, len(jacobian)])
	
	gradArray = numpy.concatenate((numpy.atleast_2d(numpy.ravel(FX)).T, numpy.atleast_2d(numpy.ravel(FY)).T), axis = 1).T 
	
	for z in range(len(jacobian)):
		SDImages[:, z] = numpy.sum(gradArray * jacobian[z], axis = 0).T
	return SDImages
#
#	[M, N] = size(FX);
#	IMGSize = M * N;
#
#	SDImages = zeros(IMGSize, 6);
#
#	for z = 1:6
#		SDImages(:, z) = sum([FX(:), FY(:)]' .* Jacobian(:, :, z))';
#		
#		%SDImages(:, z) = ([FX(:), FY(:)] * Jacobian(:, :, z))';
#	end
#
	#mshow(IMG, []);

def steepestDescentImagesAffineDisplay(SDImages, IMGSize):
	
	IMG = list()

	for z in range(SDImages.shape[1]):
		IMG.append(numpy.reshape(numpy.take(SDImages, [z], axis = 1), IMGSize))
	
	IMG = numpy.concatenate(IMG, axis = 1)

	CCSegUtils.showIMG(IMG);


def warpImageAffine(targetIMG, templateIMGSize, affineParameters):
#
#function [WarpedIMG, TX, TY] = lk_warp_image_affine(targetIMG, TemplateIMGSize, Parameters)
#
#% [WarpedImage] = lk_warp_image_affine(targetIMG, TemplateIMGSize, Parameters)
#% 
#% DESCRIPTION
#%   Warps an image onto a template using the affine parameters
#%   Parameters = [p1, p2, p3, p4, p5, p6];
#%   Assumes that the template coordinates are [1:size(TemplateIMGSize(2)), 1:size(TemplateIMGSize(1))]
#%
#% PARAMETERS
#%   targetIMG [InputM, InputN]: the input image
#%   TemplateIMGSize [2]: [M, N] the size of the template
#%
#% RETURNS
#%   WarpedImage [M, N]: the portion of the warped image that was in the boundary of the template
#%
#% In the algorithm, the template is the input and the TemplateIMGSize comes
#% from the I dont know
#%
#% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
#% See LICENSE for full license information.
#%

#
#if(numel(Parameters) ~= 6)
#    error('Parameters should have 6 elements');
#end
	assert(numpy.size(affineParameters) == 6),"Affine parameters should have 6 elements"

#[TemplateX, TemplateY] = meshgrid(1:TemplateIMGSize(2), 1:TemplateIMGSize(1));
	templateX, templateY = numpy.meshgrid(numpy.arange(1, templateIMGSize[1] + 1), numpy.arange(1, templateIMGSize[0] + 1))

#XY = cat(1, TemplateX(:)', TemplateY(:)', ones(1, prod(TemplateIMGSize)));
	XY = numpy.concatenate((numpy.atleast_2d(numpy.ravel(templateX)), numpy.atleast_2d(numpy.ravel(templateY)), numpy.ones([1, numpy.prod(templateIMGSize)])), axis = 0)
	
	transformMatrix = numpy.matrix([[1 + affineParameters[0], affineParameters[2], affineParameters[4]], [affineParameters[1], 1 + affineParameters[3], affineParameters[5]], [0, 0, 1]])

	TXY = transformMatrix * XY
	#TXY = numpy.take(TXY, [0, 1], axis = 0)

	TX = numpy.array(numpy.reshape(numpy.take(TXY, [0], axis = 0), templateIMGSize))
	TY = numpy.array(numpy.reshape(numpy.take(TXY, [1], axis = 0), templateIMGSize))

	warpedIMG = CCSegUtils.interp2q(numpy.arange(1, targetIMG.shape[1] + 1), numpy.arange(1, targetIMG.shape[0] + 1), targetIMG, TX, TY, extrapval = 0)
	
	return (warpedIMG, TX, TY)

#
#
#TransformationMatrix = [1 + Parameters(1), Parameters(3), Parameters(5); ...
#    Parameters(2), 1 + Parameters(4), Parameters(6); ...
#    0, 0, 1];
#
#TXY = TransformationMatrix * XY;
#
#TXY = TXY(1:2, :);
#
#TX = reshape(TXY(1, :), TemplateIMGSize);
#TY = reshape(TXY(2, :), TemplateIMGSize);
#
#%WarpedIMG = interp2(targetIMG, TX, TY, 'linear', 0);
#%keyboard;
#WarpedIMG = interp2q_linear_fast_c(1:size(targetIMG, 2), 1:size(targetIMG, 1), targetIMG, double(TX), double(TY), 0);
#%WarpedIMG(isnan(WarpedIMG)) = 0;
#%WarpedIMG = interp2_linear_fast(targetIMG, TX, TY);
#%WarpedIMG(isnan(WarpedIMG)) = 0;
#%keyboard;

def weightedAffineInvCompWarpCost(targetIMG, templateIMG, templateWeight, curParameters, displayStuff, targetIMGMask = None):
	#function [ImageWarpedToTemplate, TX, TY, ErrorIMG, CostValue] = lk_weighted_run_affine_inv_comp_warpcost(targetIMG, templateIMG, TemplateW, CurParameters, displayStuff)
#	[ImageWarpedToTemplate, TX, TY] = lk_warp_image_affine(targetIMG, size(templateIMG), CurParameters);
	targetIMGToTemplate, TX, TY = warpImageAffine(targetIMG, templateIMG.shape, curParameters)
	if displayStuff == True:
		pylab.subplot(4, 3, 5);	CCSegUtils.showIMG(targetIMGToTemplate); pylab.title('template coordinates warped to image');
		pylab.subplot2grid((4, 3), (1, 2), rowspan = 2, colspan = 1); pylab.cla(); CCSegUtils.showIMG(targetIMG);
			
		pylab.plot(TX[:, 0], TY[:, 0], 'b-')
		pylab.plot(TX[0, :], TY[0, :], 'b-')
		pylab.plot(TX[:, -1], TY[:, -1], 'b-')
		pylab.plot(TX[-1, :], TY[-1, :], 'b-')
		pylab.title('Coordinates on target')
	#print "oiajfdoijadsf"

	errorIMG = targetIMGToTemplate - templateIMG

	
	LKCost = numpy.sum(errorIMG * templateWeight * errorIMG)
	
	# find out if any coordinates are not in the mask
	if targetIMGMask != None:
		T = CCSegUtils.interp2q(numpy.arange(1, targetIMG.shape[1] + 1), numpy.arange(1, targetIMG.shape[0] + 1), targetIMGMask, TX, TY, extrapval = 0)
		if numpy.any(T == 0):
			LKCost = numpy.inf
	
	return (targetIMGToTemplate, TX, TY, errorIMG, LKCost)
	#
#	ErrorIMG = ImageWarpedToTemplate - templateIMG;
#
#	% simple SSD cost
#	CostValue = (ErrorIMG(:)' .* TemplateW(:)') * ErrorIMG(:);	

def coordsOfAffineWarpedTemplate(parameters, IMG, template):
	
#function [TX, TY, InterpX, InterpY] = coords_template_lk_img(Parameters, IMG, Template)

#% returns the coordinates of thetemplate transformed by the LK affine transformation matrix Parameters in the space of the target image IMG
#% InterpX and InterpY are the coordinates of the IMG in the space of Template so that
#% interp2(IMG, InterpX, InterpY) will warp the template to the space of IMG
	templateX, templateY = numpy.meshgrid(numpy.arange(1, template.shape[1] + 1), numpy.arange(1, template.shape[0] + 1))

	XY = numpy.concatenate((numpy.atleast_2d(numpy.ravel(templateX)), numpy.atleast_2d(numpy.ravel(templateY)), numpy.ones([1, numpy.size(template)])), axis = 0)

	#transformMatrix = numpy.concatenate(([1 + parameters[0], parameters[2], parameters[4]], [parameters[1], 1 + parameters[3], parameters[5]], [0, 0, 1]), axis = 0)
	transformMatrix = numpy.matrix([[1 + parameters[0], parameters[2], parameters[4]], [parameters[1], 1 + parameters[3], parameters[5]], [0, 0, 1]])
	#print transformMatrix.shape
	#print XY.shape
	TXY = numpy.matrix(transformMatrix) * numpy.matrix(XY)
	
	TX = numpy.reshape(numpy.take(TXY, [0], axis = 0), template.shape)
	TY = numpy.reshape(numpy.take(TXY, [1], axis = 0), template.shape)
	del TXY

	IMGX, IMGY = numpy.meshgrid(numpy.arange(1, IMG.shape[1] + 1), numpy.arange(1, IMG.shape[0] + 1))
	#rint numpy.atleast_2d(numpy.ravel(IMGX)).shape
	#rint numpy.atleast_2d(numpy.ravel(IMGY)).shape
	#rint numpy.ones([1, numpy.size(IMG)]).shape

	IMGXY = numpy.concatenate((numpy.atleast_2d(numpy.ravel(IMGX)), numpy.atleast_2d(numpy.ravel(IMGY)), numpy.ones([1, numpy.size(IMG)])), axis = 0)

	ITXY = scipy.linalg.solve(transformMatrix, IMGXY)

	InterpX = numpy.reshape(numpy.take(ITXY, [0], axis = 0), IMG.shape)
	InterpY = numpy.reshape(numpy.take(ITXY, [1], axis = 0), IMG.shape)
	
	return (TX, TY, InterpX, InterpY)
#% dot product method doesnt seem to work
#% % find the change in vector along the X direction
#% XChangeVector = [TX(1, 2) - TX(1, 1), TY(1, 2) - TY(1, 1)];
#% % find the change in vector along the Y direction
#% YChangeVector = [TX(2, 1) - TX(1, 1), TY(2, 1) - TY(1, 1)];
#% 
#% UnitXChangeVector = XChangeVector ./ norm(XChangeVector);
#% UnitYChangeVector = YChangeVector ./ norm(YChangeVector);	
#% 
#% % find the change in vector along the X direction
#% % scalar resolute in the direction of XChangeVector
#% 
#% [IMGX, IMGY] = meshgrid(1:size(IMG, 2), 1:size(IMG, 1));
#% 
#% XDot = (IMGX - TX(1, 1)) .* UnitXChangeVector(1) + (IMGY - TY(1, 1)) .* UnitXChangeVector(2);
#% YDot = (IMGX - TX(1, 1)) .* UnitYChangeVector(1) + (IMGY - TY(1, 1)) .* UnitYChangeVector(2);
#% 
#% InterpX = XDot ./ norm(XChangeVector);
#% InterpY = YDot ./ norm(YChangeVector);
#
#
#% %%
#% clf;
#% subplot 221;
#% imshow(IMG, []);
#% hold on;
#% plot(TX, TY, '*');
#% IMGI = interp2(Template, InterpX, InterpY, 'linear', 0);
#% subplot 222;
#% imshow(IMGI, []);
#% T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
#% T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
#% T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
#% T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');
#% IMGI = interp2(Template, ITX, ITY, 'linear', 0);
#% subplot 223;
#% imshow(IMGI, []);
#% T = line([TX(1, 1) TX(end, 1)], [TY(1, 1), TY(end, 1)]); set(T, 'Color', 'r');
#% T = line([TX(1, end) TX(end, end)], [TY(1, end), TY(end, end)]); set(T, 'Color', 'r');
#% T = line([TX(1, 1) TX(1, end)], [TY(1, 1), TY(1, end)]); set(T, 'Color', 'r');
#% T = line([TX(end, 1) TX(end, end)], [TY(end, 1), TY(end, end)]); set(T, 'Color', 'r');


def weightedAffineInvComp(targetIMG, templateIMG, templateWeight, initialParameters, numIterations, targetIMGMask = None):

# Python port of lk_run_affine_for_comp
#% [varargout] = lk_run_affine_for_comp(targetIMG, templateIMG, InitialParameters, NumIterations)
#%
#% LUCAS-KANADE algorithm Inverse Compositional method, affine
#% transformation
#% Registers the image to the template, using composition to update the
#% parameters
#%
#% PARAMETERS
#%   targetIMG [M, N]: the input image
#%   templateIMG [TM, TN]: the target image
#%   InitialParameters [6, 1]: the initial transformation
#%   NumIterations [1]: the number of iterations to perform
#%
#% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
#% See LICENSE for full license information.
#%

	displayStuff = False

	if displayStuff == True:
		pylab.clf()
		pylab.subplot(4, 3, 1); CCSegUtils.showIMG(targetIMG); pylab.title("Image");
		pylab.subplot(4, 3, 2); CCSegUtils.showIMG(templateIMG); pylab.title("Template");

	FY, FX = numpy.gradient(templateIMG)

	templateJacobian = jacobianAffine(templateIMG.shape)

	if displayStuff == True:
		pylab.subplot(4, 3, 3); CCSegUtils.showIMG(numpy.concatenate((FX, FY), axis = 1)); pylab.title('Template gradients');
		pylab.subplot(4, 3, 4); jacobianAffineDisplay(templateJacobian, templateIMG.shape); pylab.title('Template Jacobian');
	
#	% compute steepest descent images
	#[SDImages] = lk_sd_images_affine(FX, FY, Jacobian);
	SDImages = steepestDescentImagesAffine(FX, FY, templateJacobian)

#	TemplateWMatrix = spdiags(TemplateW(:), 0, numel(TemplateW), numel(TemplateW));
	templateWeightMatrix = scipy.sparse.spdiags(numpy.ravel(templateWeight), numpy.array([0]), numpy.size(templateWeight), numpy.size(templateWeight))

#	H = full(SDImages' * TemplateWMatrix * SDImages);
	Hessian = numpy.matrix(SDImages).T * templateWeightMatrix * numpy.matrix(SDImages)
	#print Hessian	
#	SDImages = full(TemplateWMatrix * SDImages);
	SDImages = templateWeightMatrix * numpy.matrix(SDImages)

	#rint templateWeightMatrix

	if displayStuff == True:
		pylab.subplot(4, 3, 8); steepestDescentImagesAffineDisplay(SDImages, templateIMG.shape); pylab.title('Steepest descent images');
		pylab.subplot(4, 3, 10); CCSegUtils.showIMG(Hessian); pylab.title('Hessian');

	curParameters = initialParameters
	(targetIMGToTemplate, TX, TY, errorIMG, LKCost) = weightedAffineInvCompWarpCost(targetIMG, templateIMG, templateWeight, curParameters, displayStuff, targetIMGMask = targetIMGMask)
	
	if numpy.isinf(LKCost):
		return (curParameters, LKCost)

	del targetIMGToTemplate; del TX; del TY;
	#[~, ~, ~, ErrorIMG, CostValue] = lk_weighted_run_affine_inv_comp_warpcost(targetIMG, templateIMG, TemplateW, CurParameters, displayStuf    f);
	
#	for CurIter = 1:NumIterations
	for curIter in range(numIterations):

		if(displayStuff == True):
			pylab.subplot(4, 3, 7);
			CCSegUtils.showIMG(errorIMG); pylab.title('Error Image');
		
		SDUpdate = scipy.linalg.solve(Hessian, numpy.matrix(SDImages).T * numpy.matrix(numpy.atleast_2d(numpy.ravel(errorIMG))).T)
		
		if(displayStuff == True):
			pylab.subplot(4, 3, 11);
			pylab.cla()
			pylab.bar(numpy.arange(0, numpy.size(SDUpdate)), SDUpdate); pylab.title('SD Updates');

#		compose the new warp
#		make the transformation matrix used to transform the coordinates
		
		curParametersMatrix = numpy.concatenate((numpy.reshape(curParameters, (2, 3), order='F'), numpy.matrix([0, 0, 1])), axis = 0)
		curParametersMatrix[0, 0] = curParametersMatrix[0, 0] + 1.0
		curParametersMatrix[1, 1] = curParametersMatrix[1, 1] + 1.0

#		CurParametersMatrix = [reshape(CurParameters, 2, 3); 0 0 1];
#		CurParametersMatrix(1, 1) = CurParametersMatrix(1, 1) + 1;
#		CurParametersMatrix(2, 2) = CurParametersMatrix(2, 2) + 1;

		SDUpdateMatrix = numpy.concatenate((numpy.reshape(SDUpdate, (2, 3), order = 'F'), numpy.matrix([0, 0, 1])), axis = 0)
		SDUpdateMatrix[0, 0] = SDUpdateMatrix[0, 0] + 1.0
		SDUpdateMatrix[1, 1] = SDUpdateMatrix[1, 1] + 1.0
		#print curParameters
		#print curParametersMatrix
		#print "SDUpdateMatrix"
		#print SDUpdateMatrix
		
		# matrix division, solve composedMatrix * SDUpdateMatrix = curParametersMatrix for composedMatrix
		# numpy has no equivalent, so solve for
		# SDUpdateMatrix^T * composedMatrix^T = curParametersMatrix
		# then traspose the result

		composedMatrix = scipy.linalg.solve(SDUpdateMatrix.T, curParametersMatrix.T).T
		
		composedMatrix[0, 0] = composedMatrix[0, 0] - 1.0
		composedMatrix[1, 1] = composedMatrix[1, 1] - 1.0
		#print "Composed Matrix"
		#print composedMatrix	
		
		composedMatrix = numpy.take(composedMatrix, [0, 1], axis = 0)
		curParameters = numpy.ravel(composedMatrix, order = 'F')
		#print curParameters

		if(displayStuff == True):
			pylab.subplot(4, 3, 12)
			pylab.cla()
			pylab.bar(numpy.arange(0, numpy.size(curParameters)), curParameters); pylab.title('Parameters');
			F = pylab.gcf()
			F.set_size_inches((20, 10), forward = True)
		
		(targetIMGToTemplate, TX, TY, errorIMG, LKCost) = weightedAffineInvCompWarpCost(targetIMG, templateIMG, templateWeight, curParameters, displayStuff, targetIMGMask = targetIMGMask)
		if numpy.isinf(LKCost):
			#print "InfCost"
			return (curParameters, LKCost)
		if displayStuff == True:
			pylab.draw()
			pylab.show(block = False)
		#quit()

		#composedMatrix = curParametersMatrix
#		SDUpdateMatrix = [reshape(SDUpdate, 2, 3); 0 0 1];
#		SDUpdateMatrix(1, 1) = SDUpdateMatrix(1, 1) + 1;
#		SDUpdateMatrix(2, 2) = SDUpdateMatrix(2, 2) + 1;
#		%SDUpdateMatrix = inv(SDUpdateMatrix);
#		
#		%ComposedMatrix = CurParametersMatrix * SDUpdateMatrix;
#		ComposedMatrix = CurParametersMatrix / SDUpdateMatrix;
#		ComposedMatrix = ComposedMatrix(1:2, :);
#		CurParameters = ComposedMatrix(:);
#		CurParameters(1) = CurParameters(1) - 1;
#		CurParameters(4) = CurParameters(4) - 1;
#		%CurParameters = CurParameters + SDUpdate;
#		
#		ParameterHistory(:, CurIter) = CurParameters;
#
#		if(displayStuff == true)
#			subplot(4, 3, 12);
#			bar(CurParameters); title('Parameters');
#			drawnow;
#		end
#		%disp([num2str(CurIter) ': ' num2str(ObjectiveFunction(CurIter))]);
#		%urParameters
#		%keyboard;
#		[~, ~, ~, ErrorIMG, CostValue] = lk_weighted_run_affine_inv_comp_warpcost(targetIMG, templateIMG, TemplateW, CurParameters, displayStuff);
#	end
#
	#pylab.show()
	return (curParameters, LKCost)
#
#	
