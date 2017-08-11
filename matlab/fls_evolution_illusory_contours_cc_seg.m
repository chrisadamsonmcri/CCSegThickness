function [LevelSetPSI, NumIter] = fls_evolution_illusory_contours_cc_seg(D, OriginalPSI, CurvatureWeight, ExclusionWeight, H, NumIter, ParameterMU, ParameterBETA, ParameterLAMBDA, BackgroundImage)

% [LevelSetPSI] = fls_illusory_contours_evolution(D, OriginalPSI, H, ParameterMU, ParameterBETA, ParameterLAMBDA, NumIter)
%
% DESCRIPTION
%	Fast Level Sets, Illusory Contour Evolution script
%	Implements the method proposed by Dong
%
% OriginalPSI is the starting shape that we are using
% U is 1 and -1 for outside and inside the objects

%[D] = fls_signed_distance_crossing_time(U, H);

if nargin < 8
	NumIter = 1000;
end

% if nargin > 7
% 	ImageSpeedTermPresent = true;
% else
% 	ImageSpeedTermPresent = false;
% end
% ImageSpeedTermPresent = false;
UpwindScheme = 'ENO3';

[BETA, GAMMA, DilationFactor] = fls_fastlocallevelset_beta_gamma(H, UpwindScheme);
[LevelSetPSI] = fls_initialise_psi(OriginalPSI, H, UpwindScheme);

%FirstPSI = LevelSetPSI;
LevelSetPSI = min(max(LevelSetPSI, -2 * GAMMA), 2 * GAMMA);

%PVD = isosurface(D, 0);
%LevelSetPSI = LevelSetPSI * H;

%ALPHA = 0.5;

%CurvatureD = ls_curvature(D, H, true(size(D)), 'gaussian');

%AbsoluteDistanceFunction = abs(D);

if length(H) == 2
	RealH = [H, 1];
else
	RealH = H;
end

RealHCell = num2cell(RealH(:));

%[FID, FJD
% FID = 123;
% FJD = 177;
% FKD = 86;

%D(FID - 1:FID + 1, FJD - 1:FJD + 1, FKD - 1:FKD + 1)
%%
%imshow(CurvatureD(:, :, 202), []); colorbar;
%%
%$[GradAbsDX, GradAbsDY, GradAbsDZ] = numeric_gradient_c(AbsoluteDistanceFunction, true(size(D)), RealH(1), RealH(2), RealH(3));
%GradCentreAbsD = cell(1, 3);
%[GradCentreAbsD{1:3}] = deal(GradAbsDX, GradAbsDY, GradAbsDZ);
%clear GradAbsDX GradAbsDY GradAbsDZ;
%[GradMinusAbsD, GradPlusAbsD] = fls_gradient_tube(AbsoluteDistanceFunction, H, true(size(D)), UpwindScheme, ones(size(D)));
%clear GChosen;

%[X, Y] = meshgrid(1:size(D, 2), 1:size(D, 1));
%ND = 2;
%GradAbsDChosen = cell(1, ND);
%for CurDim = 1:ND
%	GradAbsDChosen{CurDim} = min(GradMinusAbsD{CurDim}, 0) + max(GradPlusAbsD{CurDim}, 0);
%end
%keyboard;
%%
%imshow(D(:, :, 21) > 0, []);

%%

%SZ = size(OriginalPSI);
%NumRows = int32(SZ(1));
%NumInSlice = int32(SZ(1) * SZ(2));

%SZ = SZ([2 1 3:length(SZ)]);
%LevelSetPSI = OriginalPSI;%-ls_shape_rectangle(XYZ, [xx(10), yy(10)], [xx(end - 10), yy(end - 10)]);

%ACutoff = zeros(size(AbsoluteDistanceFunction));
%Epsilon = 4 * max(H);
%ACutoff(AbsoluteDistanceFunction <= Epsilon) = 1;
%M = AbsoluteDistanceFunction > Epsilon & AbsoluteDistanceFunction <= 2 * Epsilon;
%ACutoff(M) = (AbsoluteDistanceFunction(M) - 2 * Epsilon).^2 .* (2 * AbsoluteDistanceFunction(M) - Epsilon) ./ (Epsilon * Epsilon * Epsilon);


% 
% Filter = gaussian_filter_max_1d(2);
% for CurDim = 1:ndims(A)
% 	Y = ones(1, ndims(A));
% 	Y(CurDim) = length(Filter);
% 	A = imfilter(A, reshape(Filter, Y), 'same', 'conv', 'replicate');
% 	%A = A .* (tanh(-5 * ImageSpeedTerm + 2.8) / 2 + 0.5);
% end
% %A = smooth3(A, 'gaussian');
%clear GChosen;

%HeavisideD = 0.5 + atan(50 * D) / pi;%Hpsi=1/2*(1+2/pi*atan(D));

%P = [];
ND = ndims(LevelSetPSI);
%xx = (1:size(OriginalPSI, 2)) * RealH(1);
%yy = (1:size(OriginalPSI, 1)) * RealH(2);
%zz = (1:size(OriginalPSI, 3)) * RealH(3);
%[X, Y, Z] = meshgrid(xx, yy, zz);

%PosD = double(D > 0);

%SliceH = slice(X, Y, Z, PosD, xx(floor(length(xx) / 2)), yy(floor(length(yy) / 2)), zz(floor(length(zz) / 2)));
%[GradPSIX, GradPSIY, GradPSIZ] = numeric_gradient_c(LevelSetPSI, true(size(D)), RealH(1), RealH(2), 1);
%keyboard;
%OldMask = [];
%CH = [];
%LastPSIBeforeReinit = LevelSetPSI;
% SE = strel('arbitrary', ones(repmat(3, 1, ndims(LevelSetPSI))));
% XC = (abs(LevelSetPSI) < BETA);
% MaskJustBelowBETA = (imerode(XC, SE) == 0 & XC == 1);
% clear XC SE;

if(length(H) == 2)
	CC = contourc(abs(LevelSetPSI), [BETA BETA]);
	CC = contour2cell(CC);
	JustBelowBETAV = cat(1, CC{:});
	clear CC;
elseif(length(H) == 3)
% 	tic;
% 	[F, JustBelowBETAV] = isosurface(abs(LevelSetPSI), BETA);
% 	toc;
% 	tic;
	[~, JustBelowBETAV] = isosurface(abs(LevelSetPSI), BETA, 'noshared');
	JustBelowBETAV = unique(JustBelowBETAV, 'rows');
% 	toc;
	%clear F;
end
%keyboard;
%LastPSI = LevelSetPSI;

GridPoints = cell(1, ND);

DimensionLookup = [2, 1, 3:ND];

for CurDim = 1:ND
	GridPoints{CurDim} = 0:H(CurDim):((size(LevelSetPSI, DimensionLookup(CurDim)) - 1) * H(CurDim));
end
%keyboard;
%ReInitTrigger = false;

LastSeg = (LevelSetPSI > 0);
%[Mask, MaskIndices, AMask, DMask, NumInMask, VariableImage, IdxSubs] = init_after_reinit(LevelSetPSI, D, H, ParameterMU, GAMMA);
%[Mask, MaskIndices, AMask, DMask, HeavisideDMask, GradCentreAbsD, GradCentreAbsDMaskA, NumInMask, VariableImage, IdxSubs] = init_after_reinit(LevelSetPSI, D, H, ParameterMU, GAMMA);
[Mask, MaskIndices, DMask, HeavisideDMask, GradCentreAbsDMaskA, NumInMask, VariableImage, IdxSubs, AccelerationFunction] = init_after_reinit(LevelSetPSI, D, ExclusionWeight, CurvatureWeight, H, ParameterMU, GAMMA);

%AllMasks = Mask;
%keyboard;
for z = 1:NumIter
	
	PSIMask = LevelSetPSI(Mask);
	%AbsPSI = abs(LevelSetPSI);
	%disp(['Iteration ' num2str(z) ', ' num2str(NumInMask) ' voxels in ]);
	%disp(['Iteration ' num2str(z) ': ' num2str(NumInMask) ' voxels in mask, ' num2str(NumInMask / numel(Mask) * 100) '%']);
	
	% load D
% 	GradCentreAbsD = cell(1, 3);
% 	
% 	if(isa(D, 'char'))
% 		load(D, 'FileD');
% 		FileD = abs(FileD);
% 		[GradCentreAbsD{1}, GradCentreAbsD{2}, GradCentreAbsD{3}] = numeric_gradient_vector_return_c(FileD, Mask, RealH(1), RealH(2), RealH(3));
% 		clear FileD;
% 		load(D, 'FileD');
% 		DMask = FileD(Mask);
% 		clear FileD;
% 	else
% 		[GradCentreAbsD{1}, GradCentreAbsD{2}, GradCentreAbsD{3}] = fls_illusory_contours_gradabs_c(D, Mask, RealH(1), RealH(2), RealH(3));
% 		DMask = D(Mask);
% 	end
	%[GradCentreAbsD{1}, GradCentreAbsD{2}, GradCentreAbsD{3}] = numeric_gradient_vector_return_c(abs(D), Mask, RealH(1), RealH(2), RealH(3));
	%[GradCentreAbsD{1}, GradCentreAbsD{2}, GradCentreAbsD{3}] = fls_illusory_contours_gradabs_c(D, Mask, RealH(1), RealH(2), RealH(3));
	
% 	if(ReInitTrigger == true)
% 		keyboard;
% 	end
	
% 	if(isequal(size(MaskJustBMaskIndices = find(Mask);
	%MaskIndices = int32(MaskIndices);
	%NumInMask = length(MaskIndices);
	%disp(['Iteration ' num2str(z) ', ' num2str(NumInMask) ' voxels in ]);
	%disp(['Iteration ' num2str(z) ': ' num2str(NumInMask) ' voxels in mask, ' num2str(NumInMask / numel(Mask) * 100) '%']);
% 		MaskJustBelowBETA = MaskJustBelowBETA(Mask);
% 	end
	%Mask = true(size(LevelSetPSI));
	%SE = strel('arbitrary', ones(5 * ones(1, ndims(LevelSetPSI))));
	%Mask = imdilate(Mask, SE);
	
	%[ID, JD, KD] = find3_c(Mask);
	%EikonalMask = fast_cube_dilate(Mask, ID, JD, KD);
	
	%if(isempty(OldMask))1.0
	%	OldMask = Mask;
	%else
		%if(~isequal(OldMask, Mask))
	%		disp('doing it because the mask changed');
	%		LevelSetPSI = fls_solve_eikonal_equation(LevelSetPSI, H, 'ENO1', Mask, 1.5);
	%	end
	%	OldMask = Mask;
	%end
	
	%EikonalMask = Mask;
	%clear ID JD KD;
	
% 	hold off;
% 	imshow(LevelSetPSI(:, :, 5), []);
% 	hold on;
% 	[CCC, CH] = contour(LevelSetPSI(:, :, 3), [0 0]);
% 	set(CH, 'Color', 'b');
% 	drawnow;
% 	if(z == 10)
% 		keyboard;
% 	end
	%if((NumInMask / numel(Mask) * 100) > 15)
	%	keyboard;
	%end
	% pre cut static variables to the mask, to save time later on
	
	%HeavisideDMask = HeavisideD(MaskIndices);
	
	%ImageSpeedTermMask = 10 * ImageSpeedTerm(MaskIndices);
	%GradMinusAbsDMaskA = cell(1, ND);
	
	%$[GradAbsDX, GradAbsDY, GradAbsDZ] = numeric_gradient_c(AbsoluteDistanceFunction, true(size(D)), RealH(1), RealH(2), RealH(3));

	%clear GradCentreAbsD;
	
	%for CurDim = 1:ND
	%	IdxSubs{CurDim} = int32(IdxSubs{CurDim});
	%end
	%$disp(num2str(z));
	%keyboard;
	[GradMinusPSI, GradPlusPSI] = fls_gradient_tube_vector_return(LevelSetPSI, H, VariableImage, UpwindScheme, 1);
	

% 	keyboard;
    GradCentrePSI = cell(1, 3);
	%GradCentrePSIF = cell(1, 3);
	%GradCentrePSI2 = cell(1, 3);
	%GradCentrePSIX = cell(1, 3);
	%GradCentrePSIY = cell(1, 3);
	%GradCentrePSIZ = cell(1, 3);
	%[GradCentrePSI{1}, GradCentrePSI{2}, GradCentrePSI{3}] = numeric_gradient_c(LevelSetPSI, Mask, RealH(1), RealH(2), RealH(3));
	%tic;
	%[GradCentrePSI{1}, GradCentrePSI{2}, GradCentrePSI{3}] = numeric_gradient_vector_return_c(LevelSetPSI, Mask, RealH(1), RealH(2), RealH(3));
	%toc;
	
	%tic;
	
	[GradCentrePSI{:}] = numeric_gradient_vector_varimage_fourth_c(PSIMask, VariableImage, RealHCell{:});
	GradCentrePSI = GradCentrePSI(1:ND);
	%[GradCentrePSIF{:}] = numeric_gradient_vector_varimage_fourth_c(PSIMask, VariableImage, RealHCell{:});
	
	%subplot 221; T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSI{1}; imshow(T, []);
	%subplot 222; T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSI{2}; imshow(T, []);
	%subplot 223; T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSIF{1}; imshow(T, []);
	%subplot 224; T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSIF{2}; imshow(T, []);
	%toc;
	%keyboard;
% 	[GradCentrePSIX{1}, GradCentrePSIX{2}, GradCentrePSIX{3}] = numeric_gradient_c(GradCentrePSI{1}, Mask, RealH(1), RealH(2), RealH(3));
% 	[GradCentrePSIY{1}, GradCentrePSIY{2}, GradCentrePSIY{3}] = numeric_gradient_c(GradCentrePSI{1}, Mask, RealH(1), RealH(2), RealH(3));
% 	[GradCentrePSIZ{1}, GradCentrePSIZ{2}, GradCentrePSIZ{3}] = numeric_gradient_c(GradCentrePSI{1}, Mask, RealH(1), RealH(2), RealH(3));

	% for loop was faster than the method below, presumably less memory
	% allocation
%     tic;
	GradMAGMean = zeros(NumInMask, 1);
	
	for CurDim = 1:ND
		GradMAGMean = GradMAGMean + ...
			    (GradCentrePSI{CurDim} .* GradCentrePSI{CurDim});
	end
	GradMAGMean = realsqrt(GradMAGMean);
% 	toc;
% 	tic;
% 	T = cat(2, GradCentrePSI{:});
% 	GradMAGMean = realsqrt(sum(T .* T, 2));
% 	toc;

	%keyboard;
	
% 	tic;
% 	CurvatureMinModPSI = zeros(NumInMask, 1);
% 	for CurDim = 1:ND
%         % do minmod approximation for other dimensions with respect to the
%         % current dimension
%         % this is to compute the normalised gradients
%         
% 		% so here we are computing the forward half-point central
% 		% difference for each dimension using the other dimensions minmod
% 		% approximation
% 		
% 		CurGradMAG = zeros(NumInMask, 1);
%         for OtherDim = 1:ND
% 			if OtherDim ~= CurDim
%                 % do the minmod approximation
% 				IdxOtherS1ubs = IdxSubs;
% 				% we are moving to the right in OtherDim, so add one to the
% 				% indices in that dimension
% 				IdxOtherSubs{DimensionLookup(OtherDim)} = min(IdxOtherSubs{DimensionLookup(OtherDim)} + 1, size(Mask, DimensionLookup(OtherDim)));
% 				%OtherIndices = sub2ind(size(Mask), IdxOtherSubs{:});
% 				%inline sub2ind
% 				%OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 				switch(ND)
% 					case 1
% 						OtherIndices = IdxOtherSubs{1};
% 					case 2
% 						OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows;
% 					case 3
% 						OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 				end
% 				BadIndices = ~Mask(OtherIndices);
% 				if(any(BadIndices))
% 					IdxOtherSubs{DimensionLookup(OtherDim)}(BadIndices) = IdxOtherSubs{DimensionLookup(OtherDim)}(BadIndices) - 1;
% 					%OtherIndices = sub2ind(size(Mask), IdxOtherSubs{:});
% 					switch(ND)
% 						case 1
% 							OtherIndices = IdxOtherSubs{1};
% 						case 2
% 							OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows;
% 						case 3
% 							OtherIndices(BadIndices) = IdxOtherSubs{1}(BadIndices) + (IdxOtherSubs{2}(BadIndices) - 1) * NumRows + (IdxOtherSubs{3}(BadIndices) - 1) * NumInSlice;
% 					end
% 				end
% 				%OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 				% now do minmod
% 				%CurGradMAG = CurGradMAG + ls_minmod(GradCentrePSI{OtherDim}, GradCentrePSI{OtherDim}(VariableImage(OtherIndices))).^2;
% 				CurGradMAG = CurGradMAG + ls_minmod_squared_c(GradCentrePSI{OtherDim}, GradCentrePSI{OtherDim}(VariableImage(OtherIndices)));
%             else
%                 CurGradMAG = CurGradMAG + GradPlusPSI{CurDim} .* GradPlusPSI{CurDim};
% 			end
% 			% check for boundary conditions
% 			
%         end
%         CurGradMAG = sqrt(CurGradMAG + 1);
%         %CurGradMAG(CurGradMAG == 0) = 1;
% 		NormGradPlusPSI = GradPlusPSI{CurDim} ./ CurGradMAG;
% 		%T = zeros(size(Mask));
% 		%T(Mask) = NormGradPlusPSI;
% 		%NormGradPlusPSI = T;
% 		%clear T;
% 		%disp('NormGradPlusPSI');
% 		%NormGradPlusPSI(PixelI - 2:PixelI + 2, PixelJ - 2:PixelJ + 2)
% 		IdxOtherSubs = IdxSubs;
% 		IdxOtherSubs{DimensionLookup(CurDim)} = max(IdxOtherSubs{DimensionLookup(CurDim)} - 1, 1);
% 		%OtherIndices = sub2ind(size(Mask), IdxOtherSubs{:});
% 		%OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 		%inline sub2ind
% 		switch(ND)
% 			case 1
% 				OtherIndices = IdxOtherSubs{1};
% 			case 2
% 				OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows;
% 			case 3
% 				OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 		end
% 		BadIndices = ~Mask(OtherIndices);
% 		if(any(BadIndices))
% 			IdxOtherSubs{DimensionLookup(CurDim)}(BadIndices) = IdxOtherSubs{DimensionLookup(CurDim)}(BadIndices) + 1;
% 		end
% 		%disp('Other subs');
% 		%5[IdxOtherSubs{1}(93), IdxOtherSubs{2}(93)]
% 		%OtherIndices = sub2ind(size(Mask), IdxOtherSubs{:});
% 		% inline sud2ind
% 		switch(ND)
% 			case 1
% 				OtherIndices = IdxOtherSubs{1};
% 			case 2
% 				OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows;
% 			case 3
% 				OtherIndices = IdxOtherSubs{1} + (IdxOtherSubs{2} - 1) * NumRows + (IdxOtherSubs{3} - 1) * NumInSlice;
% 		end
% 		% now we compute the central difference gradient using the
% 		CurvatureMinModPSI = CurvatureMinModPSI + NormGradPlusPSI - NormGradPlusPSI(VariableImage(OtherIndices));
% 		%NormGradPlusPSI(PixelI - 1:PixelI + 1, PixelJ - 1:PixelJ + 1)
% 		clear IdxOtherSubs BadIndices OtherIndices NormGradPlusPSI;
% 	end
% 	toc;
% 	tic;
	
	[CurvatureMinModPSI] = fls_evolution_illusory_contours_curvminmod_c(GradCentrePSI, GradPlusPSI, VariableImage, MaskIndices, IdxSubs);
% 	toc;
% 	keyboard;
	%clear GradCentrePSI;
	
	%if(~isempty(CurvatureWeight))
	%	CurvatureCoefficients = abs(DMask) + ParameterBETA .* CurvatureWeight(Mask);
	%else
		CurvatureCoefficients = abs(DMask) + ParameterBETA;
	%end
	
	%DeltaPSI = zeros(NumVoxelsInMask, 1);
	DeltaPSI = GradMAGMean .* CurvatureMinModPSI .* CurvatureCoefficients;
	%clear CurvatureMinModPSI GradMAGMean;
% 	GradCentrePSIX{1} = GradCentrePSIX{1}(Mask);
% 	GradCentrePSIX{2} = GradCentrePSIX{2}(Mask);
% 	GradCentrePSIX{3} = GradCentrePSIX{3}(Mask);
% 	
% 	GradCentrePSIY{1} = GradCentrePSIY{1}(Mask);
% 	GradCentrePSIY{2} = GrParameter13150BETA;adCentrePSIY{2}(Mask);
% 	GradCentrePSIY{3} = GradCentrePSIY{3}(Mask);
% 	
% 	GradCentrePSIZ{1} = GradCentrePSIZ{1}(Mask);
% 	GradCentrePSIZ{2} = GradCentrePSIZ{2}(Mask);
% 	GradCentrePSIZ{3} = GradCentrePSIZ{3}(Mask);
	
	
	% before we do the cutting off of the mask, we need to do the minmod
	% approximation for each dimension
	% see the explanation on the bottom of page 587 in "Euler's elastica
	% and curvature based inpainting"
	
	
	%for CurDim = 1:ND
        %$GradMinusPSI{CurDim} = GradMinusPSI{CurDim}(Mask);
        %GradPlusPSI{CurDim} = GradPlusPSI{CurDim}(Mask);
    %    GradCentrePSI{CurDim} = GradCentrePSI{CurDim}(Mask);
	%end
	%CurG = ls_select_gradients(GradMinusPSI, GradPlusPSI, GradCentreAbsDMaskA);
	
% 	tic;
%     Hamiltonians = cell(1, ND);
% 	GradAbsDDotGradPSI = zeros(NumInMask, 1);
%     
% 	for CurDim = 1:ND
%     
%         % do grad A dot grad psi, this using the forward and backward
%         % differences for psi and A
% 		%GradCentreAbsD
% 		TempMin = max(GradCentreAbsDMaskA{CurDim}, 0);
% 		TempMax = min(GradCentreAbsDMaskA{CurDim}, 0);
% 		%CurG = ls_select_gradients(GradMinusPSI{CurDim}, GradPlusPSI{CurDim}, GradCentreAbsDMaskA);
% 		
% 		Hamiltonians{CurDim} = GradMinusPSI{CurDim} .* TempMax + GradPlusPSI{CurDim} .* TempMin;
% 		GradAbsDDotGradPSI = GradAbsDDotGradPSI + Hamiltonians{CurDim};%CurG{CurDim} .* GradCentreAbsDMaskA{CurDim};%;GradMinusPSI{CurDim} .* TempMax + GradPlusPSI{CurDim} .* TempMin; 
% 		clear TempMin TempMax;
%             %GradMinusPSI{CurDim} .* mi
	
%n(GradCentreAbsD{CurDim}(Mask), 0) + ... 
% 		    %GradPlusPSI{CurDim} .* max(GradCentreAbsD{CurDim}(Mask), 0); 
%         
%         %GChosen{CurDim} = GChosen{CurDim}(Mask);
% 		%GradMAGMean = GradMAGMean + ...
%         %   (GradMinusPSI{CurDim} .* GradMinusPSI{CurDim} + GradPlusPSI{CurDim} .* GradPlusPSI{CurDim}) / 2;
% 		
% 		%GradMAGradPlusPSI = GradMAGradPlusPSI + GradPlusPSI{CurDim} .* GradPlusPSI{CurDim};
% 	end
% 	toc;
% 	tic;
	[Hamiltonians, GradAbsDDotGradPSI] = fls_evolution_illusory_contours_hamiltonians_c(GradMinusPSI, GradPlusPSI, GradCentreAbsDMaskA);
%	toc;
	
	DeltaPSI = DeltaPSI + GradAbsDDotGradPSI;
	%clear GradAbsDDotGradPSI;
% 	F = find(abs(GradAbsDDotGradPSI) > 100);
% 	if(~isempty(F))
% 
% 		G = zeros(length(F), 3);
% 		for CurDim = 1:ND
% 
% 			% do grad A dot grad psi, this using the forward and backward
% 			% differences for psi and A
% 			%GradCentreAbsD
% 			TempMin = min(GradCentreAbsDMaskA{CurDim}, 0);
% 			TempMax = max(GradCentreAbsDMaskA{CurDim}, 0);
% 			%CurG = ls_select_gradients(GradMinusPSI{CurDim}, GradPlusPSI{CurDim}, GradCentreAbsDMaskA);
% 			G(:, (CurDim - 1) * 4 + 1) = TempMin(F);
% 			G(:, (CurDim - 1) * 4 + 2) = TempMax(F);
% 			G(:, (CurDim - 1) * 4 + 3) = GradMinusPSI{CurDim}(F);
% 			G(:, (CurDim - 1) * 4 + 4) = GradPlusPSI{CurDim}(F);
% 			%GradMinusPSI{CurDim} .* TempMax + GradPlusPSI{CurDim} .* TempMin;
% 			%Hamiltonians{CurDim} = GradMinusPSI{CurDim} .* TempMax + GradPlusPSI{CurDim} .* TempMin;
% 			%GradAbsDDotGradPSI = GradAbsDDotGradPSI + Hamiltonians{CurDim};%CurG{CurDim} .* GradCentreAbsDMaskA{CurDim};%;GradMinusPSI{CurDim} .* TempMax + GradPlusPSI{CurDim} .* TempMin; 
% 			clear TempMin TempMax;
% 		end
% 		disp(G);
% 	end
	%clear CurG;
%	tic;
%     GodunovMinusPSI = zeros(NumInMask, 1);
% 	GodunovPlusPSI = zeros(NumInMask, 1);
% 	for CurDim = 1:ND
% 		GodunovMinusPSI = GodunovMinusPSI + max(min(GradMinusPSI{CurDim}, 0).^2, max(GradPlusPSI{CurDim}, 0).^2);
% 		GodunovPlusPSI = GodunovPlusPSI + max(max(GradMinusPSI{CurDim}, 0).^2, min(GradPlusPSI{CurDim}, 0).^2);
% 	end
% 	GodunovMinusPSI = sqrt(GodunovMinusPSI);
% 	GodunovPlusPSI = sqrt(GodunovPlusPSI);
% % \	if(ImageSpeedTermPresent)
% % 		HeavisideDelta = max(HeavisideDMask + ImageSpeedTermMask, 0) .* GodunovPlusPSI + min(HeavisideDMask + ImageSpeedTermMask, 0) .* GodunovMinusPSI;
% % 	elsefls_evolution_illusory_contours_heavisidedelta_c
% 	HeavisideDelta = max(HeavisideDMask, 0) .* GodunovPlusPSI + min(HeavisideDMask, 0) .* GodunovMinusPSI;
% 	
	
	HeavisideDelta = fls_evolution_illusory_contours_heavisidedelta_c(GradMinusPSI, GradPlusPSI, HeavisideDMask);
	%toc;
% 	keyboard;
	
	%Original, only prevents going into -D
	%DeltaPSI = DeltaPSI - ParameterLAMBDA .* HeavisideDelta);
	%CC blows internal boundaries outwards so that poorly initialised CC do not go to the wrong boundary
	DeltaPSI = DeltaPSI - ParameterLAMBDA .* (HeavisideDelta + D(Mask) * 5  + ExclusionWeight(Mask) * 15);
% 	clear GodunovMinusPSI GodunovPlusPSI HeavisideDelta;
% 	clear GradMinusPSI GradPlusPSI;
% 	end

	%for CurDim = 1:ND
	%	Hamiltonians{CurDim} = Hamiltonians{CurDim} - ParameterLAMBDA .* HeavisideDelta;
		%Hamiltonians{CurDim} = Hamiltonians{CurDim} + ImageSpeedTermMask;
	%end

	C = zeros(NumInMask, 1);
	AbsPSI = abs(LevelSetPSI(MaskIndices));
	C(AbsPSI <= BETA) = 1;
	M = AbsPSI > BETA & AbsPSI <= GAMMA;
	
	C(M) = (AbsPSI(M) - GAMMA).^2 .* (2 .* AbsPSI(M) + GAMMA - 3 .* BETA) ./ (GAMMA - BETA).^3;
	%clear AbsPSI M;
	DeltaPSI = DeltaPSI .* C;
	%clear C;
	%S = Hamiltonians;
	
% 	tic;
% 	for CurDim = 1:ND
% 		T = zeros(size(Mask));
% 		T(Mask) = Hamiltonians{CurDim};
% 		%$GradT = cell(1, ND);
% 		Hamiltonians{CurDim} = numeric_gradient_vector_return_c(T, Mask, RealH(1), RealH(2), RealH(3), CurDim);
% 		%Hamiltonians{CurDim} = GradT(Mask);
% 		clear T;
% 	end
% 	toc;
% 	tic;
% 	S2 = cell(1, ND);
	for CurDim = 1:ND
		%T = zeros(size(Mask));
		%T(Mask) = Hamiltonians{CurDim};
		Hamiltonians{CurDim} = numeric_gradient_vector_varimage_c(Hamiltonians{CurDim}, VariableImage, RealHCell{:}, CurDim);
		%Hamiltonians{CurDim} = GradT(Mask);
		%clear T;
	end
% 	toc;
% 	keyboard;
	% TimeStep = 1 ./ max(H1/DeltaX + H2/DeltaY + H3/DeltaZ +
	% 2b/(DeltaX.^2) + 2b/(DeltaY.^2) + 2b/(DeltaZ.^2))
	
	%$CurvatureCoefficients = ParameterBETA;%AbsoluteDistanceFunction(MaskIndices) + ParameterBETA;
	MaxHamiltonian = zeros(NumInMask, 1);
	for CurDim = 1:ND
		MaxHamiltonian = MaxHamiltonian + abs(Hamiltonians{CurDim}) ./ H(CurDim);
		Hamiltonians{CurDim} = [];
		MaxHamiltonian = MaxHamiltonian + 2 .* CurvatureCoefficients ./ H(CurDim) ./ H(CurDim);
	end
	%clear CurvatureCoefficients;
	%MaxHamiltonian = MaxHamiltonian .* C;
	%MaxHamiltonian = max(MaxHamiltonian);
% 	SE = strel('arbitrary', 25 * ones(repmat(1, 1, ndims(LevelSetPSI))));
% 	T = zeros(size(LevelSetPSI));
% 	T(MaskIndices) = MaxHamiltonian;
% 	T = imdilate(T, SE);
% 	MaxHamiltonian = T(MaskIndices);
	MaxHamiltonian = max(MaxHamiltonian);
	TimeStep = 0.25 ./ MaxHamiltonian;
    
	%FirstTerm = GradMAGMean .* CurvatureMinModPSI .* CurvatureCoefficients;

    clear AbsPSI M;
	
	%DeltaPSI = TimeStep .* C .* (GradMAGMean .* CurvatureMinModPSI .* CurvatureCoefficients - GradAbsDDotGradPSI - ParameterLAMBDA .* HeavisideDelta);
	%DeltaPSI = TimeStep .* C .* (GradMAGMean .* CurvatureMinModPSI .* CurvatureCoefficients - ParameterLAMBDA .* HeavisideDelta);
	DeltaPSI = DeltaPSI .* TimeStep;% .* (GradMAGMean .* CurvatureMinModPSI .* CurvatureCoefficients + GradAbsDDotGradPSI - ParameterLAMBDA .* HeavisideDelta);
	%clear C;
	%LastPSI = LevelSetPSI;
	%keyboard;
% 	clf;
% 	SR = 3;
% 	SC = 5;
%  	subplot(SR, SC, 1); T = zeros(size(OriginalPSI)); T(Mask) = GradMinusPSI{1}; imshow(T, []); title('GX-');
%  	subplot(SR, SC, 2); T = zeros(size(OriginalPSI)); T(Mask) = GradMinusPSI{2}; imshow(T, []); title('GY-');
%  	subplot(SR, SC, 3); T = zeros(size(OriginalPSI)); T(Mask) = GradPlusPSI{1}; imshow(T, []); title('GX+');
%  	subplot(SR, SC, 4); T = zeros(size(OriginalPSI)); T(Mask) = GradPlusPSI{2}; imshow(T, []); title('GY+');
% 	subplot(SR, SC, 5); T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSI{1}; imshow(T, []); title('GXCentre');
%  	subplot(SR, SC, 6); T = zeros(size(OriginalPSI)); T(Mask) = GradCentrePSI{2}; imshow(T, []); title('GYCentre');
% 	subplot(SR, SC, 7); T = zeros(size(OriginalPSI)); T(Mask) = CurvatureMinModPSI; imshow(T, []); title('CurvatureMinMod');
% 	subplot(SR, SC, 8); T = zeros(size(OriginalPSI)); T(Mask) = GradMAGMean; imshow(T, []); title('gradMAGmean');
% 	subplot(SR, SC, 9); hold off; T = zeros(size(OriginalPSI)); T(Mask) = DeltaPSI; imshow(T, []); title('DeltaPSI');
% 	hold on;
% 	[~, CC] = contour(LevelSetPSI, [0, 0]);
% 	set(CC, 'color', 'r');
%  	subplot(SR, SC, 10); imshow(Mask, []); title('Mask');
%  	subplot(SR, SC, 11); zeros(size(OriginalPSI)); T = LevelSetPSI; imshow(T, []); title('psi');
% 	subplot(SR, SC, 12); hold on; imshow(BackgroundImage, []); title('IMG');
% 	hold on;
% 	[~, CC] = contour(LevelSetPSI, [0, 0]);
% 	set(CC, 'color', 'r');
% 	subplot(SR, SC, 13); zeros(size(OriginalPSI)); T(Mask) = C; imshow(T, []); title('C');
% 	subplot(SR, SC, 14); zeros(size(OriginalPSI)); imshow(AccelerationFunction, []); title('A');
% 	subplot(SR, SC, 14); zeros(size(OriginalPSI)); imshow(D < 0, []); title('D');
% 	drawnow;
% 	disp(num2str(max(abs(DeltaPSI))));
	%keyboard;
	LevelSetPSI(MaskIndices) = LevelSetPSI(MaskIndices) + DeltaPSI;
	
	%clear DeltaPSI;
	%OriginalSeg = OriginalPSI > 0;
%	Slices = [100, 120, 140, 160];
%	for CurSlice = 1:length(Slices)
%		subplot(2, 2, CurSlice);
%		hold off;
%		imshow(D(:, :, Slices(CurSlice)) < 0, []);
%		hold on;
%		[C, CH] = contour(LevelSetPSI(:, :, Slices(CurSlice)), [0 0]);
%		set(CH, 'Color', 'r');
%		[C, CH] = contour(OriginalPSI(:, :, Slices(CurSlice)), [0 0]);
%		set(CH, 'Color', 'b');
%
%% 		BW = LevelSetPSI(:, :, Slices(CurSlice)) > 0;
%% 		RGB = zeros([size(BW) 3]);
%% 		RGB(:, :, 1) = double(BW);
%% 		image(RGB, 'AlphaData', double(BW) * 0.5);
%% 		RGB = zeros([size(BW) 3]);
%% 		RGB(:, :, 2) = double(OriginalSeg(:, :, Slices(CurSlice)));
%% 		image(RGB, 'AlphaData', double(OriginalSeg(:, :, Slices(CurSlice))) * 0.5);
%		drawnow;
%	end
 	%disp(num2str(z));
 	%if(mod(z, 50) == 0)
	%if(true)
% 		PVPSI = isosurface(LevelSetPSI, 0);
% 		clf;
% 		subplot 121;
% 		patch(PVD, 'FaceColor', 'r', 'EdgeColor', 'none');
% 		axis equal;
% 		hold on;
% 		patch(PVPSI, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% 		lighting gouraud; light;
% 		view(3);
% 		
% 		subplot 122;
% 		patch(PVPSI, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 1);
% 		lighting gouraud; light;
% 		view(3);
% 		axis equal;
% 		drawnow;
		
% 		T = zeros(size(D));
% 		T(Mask) = DeltaPSI;
% 		subplot 221;
% 		hold off;
% 		imshow(T, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% 		subplot 222;
% 		hold off;
% 		imshow(D > 0, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% % 		subplot 222;
% % 		hold off;
% % 		imshow(LevelSetPSI > 0, []);
% % 		hold on;
% % 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% % 		set(CH, 'Color', 'r');
% 		subplot 223;
% 		hold off;
% 		imshow(BackgroundImage, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% 		subplot 224;
% 		hold off;
% 		imshow(min(AccelerationFunction, 20), []);
% 		%imshow(A, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% 		
% 		%keyboard;
% 		%subplot 223;
% 		%hold off;
% % 		imshow(A, []);
% % 		hold on;
% % 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% % 		set(CH, 'Color', 'r');
% % 		drawnow;
% % 		subplot 224;
% % 		hold off;
% % 		imshow(T, []);
% %  		hold on;
% %  		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% %  		set(CH, 'Color', 'r');
%  		drawnow;
% %		keyboard;
%   	end
% 	[GradPSIX, GradPSIY, GradPSIZ] = numeric_gradient_c(LevelSetPSI, Mask, RealH(1), RealH(2), 1);
% 	GradPSIX = GradPSIX(Mask);
% 	GradPSIY = GradPSIY(Mask);
% 	GradPSIZ = GradPSIZ(Mask);
% 	GradMAG = sqrt(GradPSIX .* GradPSIX + GradPSIY .* GradPSIY + GradPSIZ .* GradPSIZ);
	
	%M = mean(GradMAG);
	%MaxError = max(abs(GradMAG - 1));
	%MaxG = max(GradMAG);
	
	%[MinDistance] = fls_distance_to_beta(GridPoints, LevelSetPSI, LastPSIBeforeReinit, H, BETA);
	%disp([num2str(M) '+' num2str(MaxError) '+' num2str(MaxG) ' + ' num2str(MinDistance)]);	
	%disp([num2str(M) '+' num2str(MaxError) '+' num2str(MaxG) ' + '
	%disp(num2str(MinDistance));	
	% the front is too close to the cutoff part of the speed function,
% 	% reinitialise
%  		if(mod(z, 100) == 0)
%  		clf;
%  		isosurface(LevelSetPSI, 0);
% %%
% 
% 		Slice = 94;
% 		SR = 2;
% 		SC = 3;
% 		%subplot(SR, SC, 1);
% 		hold off; T = zeros(size(Mask)); T = D > 0; imshow(T(:, :, Slice), []); title('LevelSetPSI'); colorbar;
% 		hold on;
% 		%delete(CH);
% 		%delete(CH2);
% 		[CCC, CH2] = contour(LevelSetPSI(:, :, Slice), [-BETA BETA]);
% 		set(CH2, 'Color', 'b');
% 		[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'r');
% 		[CCC, CH] = contour(LastPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'g');
% 		[CCC, CH] = contour(OriginalPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'y');
% %%
% % 		subplot(SR, SC, 2); T = zeros(size(Mask)); T(Mask) = GradMAGMean; imshow(T(:, :, Slice), []); colorbar;
% % 		subplot(SR, SC, 3); T = zeros(size(Mask)); T(Mask) = CurvatureMinModPSI; imshow(T(:, :, Slice), []); colorbar;
% % 		subplot(SR, SC, 4); T = zeros(size(Mask)); T(Mask) = CurvatureCoefficients; imshow(T(:, :, Slice), []); colorbar;
% % 		subplot(SR, SC, 5); T = zeros(size(Mask)); T(Mask) = GradAbsDDotGradPSI; imshow(T(:, :, Slice), []); colorbar;
% % 		subplot(SR, SC, 6); T = zeros(size(Mask)); T(Mask) = HeavisideDelta; imshow(T(:, :, Slice), []); colorbar;
% 		drawnow;
% 		keyboard;
% %%
	[NeedToReinit] = fls_check_for_reinit(LevelSetPSI, JustBelowBETAV, H, BETA);
	
	if(mod(z, 50) == 0)
% 		Incr = 10;
% 		MidSlice = floor(size(LevelSetPSI, 3) / 2);
%  		Slices = MidSlice - (Incr * 2):Incr:MidSlice + (Incr * 2);
% 		
% 		for CurSlice = 1:length(Slices)
% 			subplot(1, 5, CurSlice);
% 			hold off;
% 			if(isa(D, 'char'))
% 				load(D, 'FileD');
% 				imshow(FileD(:, :, Slices(CurSlice)) > 0, []);
% 				clear FileD;
% 			else
% 				imshow(D(:, :, Slices(CurSlice)) > 0, []);
% 			end
% 			
% 			%imshow(LevelSetPSI(:, :, Slices(CurSlice)), []);
% 			hold on;
% 			[CC, CCG] = contour(LevelSetPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'g');
% 			%[CC, CCG] = contour(LastPSI(:, :, Slices(CurSlice)), [0 0]);
% 			%set(CCG, 'Color', 'r');
% 			[CC, CCG] = contour(OriginalPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'b');
% 		end
% 		drawnow;
% 		
		NumChanged = sum(LastSeg(:) ~= (LevelSetPSI(:) > 0));
		
		%disp(['Voxels changed since last: ' num2str(NumChanged) ' of ' num2str(NumChanged / numel(LevelSetPSI))]);
		PercentageThreshold = 5e-4;
		if((NumChanged / numel(LevelSetPSI)) < PercentageThreshold && z > (NumIter / 2))
			%disp('Converged, terminating');
			break;
		end
		clear PercentageThreshold;
		LastSeg = (LevelSetPSI > 0);
	end
	if(NeedToReinit || mod(z, 20) == 0)

		%disp('Reinitialising');
		
		%LastPSIBeforeReinit = LevelSetPSI;
		%%%% inlined to save memory
		%LevelSetPSI = fls_reinitialise_psi(LevelSetPSI, H, UpwindScheme);
		%%%%%%%%%%%%%%% inlining start
		[BETA, GAMMA, DilationFactor] = fls_fastlocallevelset_beta_gamma(H, UpwindScheme);
		
		%disp('about to compute distance function');
		% make a mask that will speed things up
		SE = strel('arbitrary', ones(repmat(DilationFactor * 2 + 1, 1, ndims(OriginalPSI))));
		ReinitMask = (abs(LevelSetPSI) < 2 * GAMMA);
		ReinitMask = imdilate(ReinitMask, SE);
		%HigherNeighbours = imerode(OriginalPSI, SE) ~= OriginalPSI;
		%Mask = LowerNeighbours | HigherNeighbours;
		
		%OriginalPSI(~MaskBinary) = -1;%2 .* GAMMA .* sign(OriginalPSI(~Mask));
		ReinitMask = double(ReinitMask);
		% set the outside points to -1
		ReinitMask = ReinitMask - 1;
		%Mask(Mask == 0) = -1;

		% if OriginalPSI is the output of sign()
		% then initialise it to a signed distance function[D] =
		% 
		%LevelSetPSI = fls_signed_distance_crossing_time(OriginalPSI, H, Inf);
		%%
		%imshow(Mask(:, :, 154), []);
		%%

		S = sign(LevelSetPSI);
		
		if(ndims(LevelSetPSI) == 2)
			LevelSetPSI = sign(LevelSetPSI);
			LevelSetPSI(LevelSetPSI == 0) = -1;
			LevelSetPSI = computeDistanceFunction2d(LevelSetPSI, H, ReinitMask, 2);
		elseif(ndims(LevelSetPSI) == 3)
			%disp('Reinitialising: computing distance function');
			% constructing this using a non-signed distance function seems to stuff
			% up with lsmlib, so revert to using a signed version
			LevelSetPSI = sign(LevelSetPSI);
			LevelSetPSI(LevelSetPSI == 0) = -1;
			LevelSetPSI = computeDistanceFunction3d(LevelSetPSI, H, ReinitMask, 2);
		end
		
		%clear Mask;
		ReinitMask = (ReinitMask == 0);
		%LevelSetPSI(~ReinitMask) = 2 .* GAMMA .* sign(LevelSetPSI(~ReinitMask));
		LevelSetPSI(~ReinitMask) = 2 .* GAMMA .* S(~ReinitMask);
		LevelSetPSI = max(min(LevelSetPSI, 2 * GAMMA), -2 * GAMMA);
		
		clear SE ReinitMask S;
		%%%%%%%%%%%%%%%% inlining finished
		
		%Mask = abs(LevelSetPSI) <= GAMMA;
		if(length(H) == 2)
			CC = contourc(abs(LevelSetPSI), [BETA BETA]);
			CC = contour2cell(CC);
			JustBelowBETAV = cat(1, CC{:});
			clear CC;
		elseif(length(H) == 3)
			%[F, JustBelowBETAV] = isosurface(abs(LevelSetPSI), BETA);
			[~, JustBelowBETAV] = isosurface(abs(LevelSetPSI), BETA, 'noshared');
			JustBelowBETAV = unique(JustBelowBETAV, 'rows');
			%clear F;
		end
		%[Mask, MaskIndices, NumInMask, VariableImage, IdxSubs] = init_after_reinit(LevelSetPSI, D, H, GAMMA);
		%[Mask, MaskIndices, AMask, DMask, NumInMask, VariableImage, IdxSubs] = init_after_reinit(LevelSetPSI, D, H, ParameterMU, GAMMA);
		[Mask, MaskIndices, DMask, HeavisideDMask, GradCentreAbsDMaskA, NumInMask, VariableImage, IdxSubs, AccelerationFunction] = init_after_reinit(LevelSetPSI, D, ExclusionWeight, CurvatureWeight, H, ParameterMU, GAMMA);
		%AllMasks = AllMasks | Mask;
% 		Mask = abs(LevelSetPSI) <= GAMMA;
% 		MaskIndices = find(Mask);
% 		MaskIndices = int32(MaskIndices);
% 		NumInMask = length(MaskIndices);
		
		%[F, JustBelowBETAV] = isosurface(abs(LevelSetPSI), BETA);
		
		%SE = strel('arbitrary', ones(repmat(3, 1, ndims(LevelSetPSI))));
		%XC = (abs(LevelSetPSI) < BETA);
		%MaskJustBelowBETA = (imerode(XC, SE) == 0 & XC == 1);
		%MaskJustBelowBETA = MaskJustBelowBETA(Mask);
		
		%clear XC SE;
		%disp('reinitialise');
	end
	%disp('hello');
% 	if(mod(z, 10) == 0)
% 		clf;
% 		subplot 211;
% 		imshow(BackgroundImage, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% 		drawnow;
% 		subplot 212;
% 		imshow(D > 0, []);
% 		hold on;
% 		[CC, CH] = contour(LevelSetPSI, [0, 0]);
% 		set(CH, 'Color', 'r');
% 		drawnow;
% 		%keyboard;
% 	end
% 		clf;
		%NumChanged = sum(sign(LevelSetPSI(:)) ~= sign(LastPSI(:)));
% 		AX = zeros(1, 3);
% 		AX(1) = subplot(3, 3, 1);
% 		AX(2) = subplot(3, 3, 2);
% 		AX(3) = subplot(3, 3, 3);
% 		display_contours(AX, AVW, LevelSetPSI);
% 		AX = zeros(1, 3);
% 		AX(1) = subplot(3, 3, 4);
% 		AX(2) = subplot(3, 3, 5);
% 		AX(3) = subplot(3, 3, 6);
% 		display_contours(AX, A .* Mask, LevelSetPSI);
% 		drawnow;
		% want to display A, LevelSetPSI == 0, 
		
		%disp(num2str(NumChanged));
		%if(NumChanged < 5)
		%	break;
		%end
		%Slices = 1;
% 		AXSlice = floor(size(LevelSetPSI, 3) / 2);
% 		CORSlice = floor(size(LevelSetPSI, 1) / 2);
% 		SAGSlice = floor(size(LevelSetPSI, 2) / 2);
% 		
% 		subplot 131;
% 		hold off;
% 		imshow(D(:, :, AXSlice) > 0, []);
% 		hold on;
% 		[CC, CCG] = contour(LevelSetPSI(:, :, AXSlice), [0 0]);
% 		set(CCG, 'Color', 'g');
% 		[CC, CCG] = contour(LastPSI(:, :, AXSlice), [0 0]);
% 		set(CCG, 'Color', 'r');
% 		[CC, CCG] = contour(OriginalPSI(:, :, AXSlice), [0 0]);
% 		set(CCG, 'Color', 'b');
% 		
% 		subplot 132;
% 		hold off;
% 		imshow(squeeze(D(:, SAGSlice, :)) > 0, []);
% 		hold on;
% 		[CC, CCG] = contour(squeeze(LevelSetPSI(:, SAGSlice, :)), [0 0]);
% 		set(CCG, 'Color', 'g');
% 		[CC, CCG] = contour(squeeze(LastPSI(:, SAGSlice, :)), [0 0]);
% 		set(CCG, 'Color', 'r');
% 		[CC, CCG] = contour(squeeze(OriginalPSI(:, SAGSlice, :)), [0 0]);
% 		set(CCG, 'Color', 'b');
% 		subplot 133;
% 		hold off;
% 		imshow(squeeze(D(CORSlice, :, :)) > 0, []);
% 		hold on;
% 		[CC, CCG] = contour(squeeze(LevelSetPSI(CORSlice, :, :)), [0 0]);
% 		set(CCG, 'Color', 'g');
% 		[CC, CCG] = contour(squeeze(LastPSI(CORSlice, :, :)), [0 0]);
% 		set(CCG, 'Color', 'r');
% 		[CC, CCG] = contour(squeeze(OriginalPSI(CORSlice, :, :)), [0 0]);
% 		set(CCG, 'Color', 'b');
% 

		%Seed = [18, 39, 45];
		%LevelSetPSI(Seed(1) - 1:Seed(1) + 1, Seed(2) - 1:Seed(2) + 1, Seed(3) - 1:Seed(3) + 1)
		
		% need to display
		% GradCentrePSI
		% Grad
		% GodunovMinusPSI
		% GodunovPlusPSI
		% HeavisideDelta
		% DeltaPSI
		%keyboard;
		%%
% 		Slice = 20;
% 		
% 		SR = 4;
% 		SC = 3;
% 		subplot(SR, SC, 1); hold off; S = zeros(size(LevelSetPSI)); S(Mask) = GradCentrePSI{1}; imshow(S(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 2); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradCentrePSI{2}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 3); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradCentrePSI{3}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 4); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradMinusPSI{1}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 5); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradMinusPSI{2}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 6); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradMinusPSI{3}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 7); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradPlusPSI{1}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 8); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradPlusPSI{2}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 9); hold off; T = zeros(size(LevelSetPSI)); T(Mask) = GradPlusPSI{3}; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 10); hold off; T = zeros(size(LevelSetPSI)); T = LastPSI; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		subplot(SR, SC, 11); hold off; T = zeros(size(LevelSetPSI)); T = sign(LevelSetPSI); imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		[FX, FY, FZ] = gradient(LastPSI, H(1), H(2), H(3));
% 		subplot(SR, SC, 12); hold off;  T = zeros(size(LevelSetPSI)); T = FX; T(~Mask) = 0; T = T - S; imshow(T(:, :, Slice), []); colorbar; hold on; [CC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]); set(CH, 'Color', 'r');
% 		drawnow;
		%%
		%size(AVW)
% 		Incr = 10;
% 		MidSlice = floor(size(D, 3) / 2);
%  		Slices = MidSlice - (Incr * 2):Incr:MidSlice + (Incr * 2);
% 		
% 		for CurSlice = 1:length(Slices)
% 			subplot(1, 5, CurSlice);
% 			hold off;
% 			imshow(D(:, :, Slices(CurSlice)) > 0, []);
% 			%imshow(LevelSetPSI(:, :, Slices(CurSlice)), []);
% 			hold on;
% 			[CC, CCG] = contour(LevelSetPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'g');
% 			%[CC, CCG] = contour(LastPSI(:, :, Slices(CurSlice)), [0 0]);
% 			%set(CCG, 'Color', 'r');
% 			[CC, CCG] = contour(OriginalPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'b');
% 		
% 		end
% 		drawnow;
		%keyboard;
% 		
% 		for CurSlice = 1:length(Slices)
% 			subplot(2, 5, CurSlice);
% 			hold off;
% 			imshow(D(:, :, Slices(CurSlice)) > 0, []);
% 			hold on;
% 			[CC, CCG] = contour(LevelSetPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'g');
% 			[CC, CCG] = contour(LastPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'r');
% 			[CC, CCG] = contour(OriginalPSI(:, :, Slices(CurSlice)), [0 0]);
% 			set(CCG, 'Color', 'b');
% 		end
% 		drawnow;
		%LastPSI = LevelSetPSI;
% 		cla;
% 		isosurface(LevelSetPSI, 0);
% 		axis equal;
% 		xlim([1 size(LevelSetPSI, 2)]);
% 		ylim([1 size(LevelSetPSI, 1)]);
% 		zlim([1 size(LevelSetPSI, 3)]);
% 		lighting gouraud;
% 		light;
% 		view(3);
% 		drawnow;
		
		%ylim([GridPoints{2}(1) GridPoints{2}(end)]);
		%zlim([GridPoints{3}(1) GridPoints{3}(end)]);
	%end
	NumIter = z;
end

% 	
	%end
%%
	%if(z > 10)
	%	keyboard;
	%end
%%
	
	%XL = [285.6932  304.8182];
	%YL = [195.0718  212.6343];
	%XL = [0.5000   68.1866];
	%YL = [186.6982  248.8548];

% 	SC = 3;
% 	SR = 3;
% 	subplot(SR, SC, 1); T = zeros(size(Mask)); T(Mask) = FirstTerm; imshow(T(:, :, Slice), []); title('First Term');
% 	subplot(SR, SC, 2); T = zeros(size(Mask)); T = log(abs(D)); imshow(T(:, :, Slice), []); title('A');%T(Mask) = HeavisideDelta; imshow(T(:, :, Slice), []); title('HDelta');
% 	subplot(SR, SC, 3); T = zeros(size(Mask)); T(Mask) = GradAbsDDotGradPSI; imshow(T(:, :, Slice), []); title('GradAbsDDotGradPSI');
%  	hold on;
%  	[I, J] = find(Mask);
%  	H1 = quiver(J, I, GradCentreAbsDMaskA{1}, GradCentreAbsDMaskA{2});
% 	set(H1, 'Color', 'b');
% 	H2 = quiver(J, I, CurG{1}, CurG{2});
% 	set(H2, 'Color', 'g');
	
	%GradAbsDDotGradPSI = GradAbsDDotGradPSI + CurG{CurDim} .* GradCentreAbsDMaskA{CurDim};%;GradMinusPSI{CurDim} .* TempMin + GradPlusPSI{CurDim} .* TempMax; 
	
% 	subplot(SR, SC, 4); T = zeros(size(Mask)); T(Mask) = CurvatureMinModPSI; imshow(T(:, :, Slice), []); title('Curvature');
% 	subplot(SR, SC, 5); T = zeros(size(Mask)); T = CurvatureD; imshow(T(:, :, Slice), []); title('D');%T(Mask) = CurvatureCoefficients; imshow(T(:, :, Slice), []); title('Curvature Coefficients');
% 	subplot(SR, SC, 6); T = zeros(size(Mask)); T(Mask) = GradMAGMean; imshow(T(:, :, Slice), []); title('GradMAG');
% 	subplot(SR, SC, 7); T = zeros(size(Mask)); T(Mask) = C; imshow(T(:, :, Slice), []); title('C');
% 	subplot(SR, SC, 8); T = zeros(size(Mask)); T = LevelSetPSI; imshow(T(:, :, Slice), []); title('LevelSetPSI');
% 	hold on;
% 	%T = zeros(size(Mask)); T(Mask) = 1;
% 	%[CCC, CH] = contour(T(:, :, Slice), [0.5 0.5]);
% 	%set(CH, 'Color', 'b');
% 	[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [GAMMA GAMMA]);
% 	set(CH, 'Color', 'g');
% 	[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [-GAMMA -GAMMA]);
% 	set(CH, 'Color', 'g');
% 	[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]);
% 	set(CH, 'Color', 'r');
% 	[CCC, CH] = contour(OriginalPSI(:, :, Slice), [0 0]);
% 	set(CH, 'Color', 'y');
	%subplot(SR, SC, 9); 

	%streamline(XY);

	%XL = [127.3868  165.6368];
	%YL = [246.3750  281.5000];
	%for CurAxis = 1:9
	%	subplot(SR, SC, CurAxis);
		%xlim(XL);
		%ylim(YL);
	%	colorbar;
	%end
%	drawnow;
	% too close to BETA, reinitialise
% 	if(MinDistance < max(H))
% 		if(ndims(LevelSetPSI) == 2)
%  			LevelSetPSI = computeDistanceFunction2d(LevelSetPSI, H, [], 2);
%  		else
%  			LevelSetPSI = computeDistanceFunction3d(LevelSetPSI, H, [], 2);
% 		end
% 		LevelSetPSI = min(max(LevelSetPSI, -2 * GAMMA), 2 * GAMMA);
% 		LastPSIBeforeReinit = LevelSetPSI;
% 		disp('reinitialise');
% 	end
% 
	
% 	if(MaxError > GAMMA)
% 		%U = sign(LevelSetPSI);
% 		%U(U == 0) = -1;
% 		OldPSI = LevelSetPSI;
% 		%[LevelSetPSI] = fls_signed_distance_crossing_time(U, H, 2 * GAMMA);
% 		if(ndims(LevelSetPSI) == 2)
% 			LevelSetPSI = computeDistanceFunction2d(LevelSetPSI, H, [], 2);
% 		else
% 			LevelSetPSI = computeDistanceFunction3d(LevelSetPSI, H, [], 2);
% 		end
% 	
% 		%[LevelSetPSI] = fls_signed_distance_crossing_time(U, H, 2 * GAMMA);
% 		subplot(SR, SC, 9); hold off; T = HeavisideDelta; imshow(T(:, :, Slice), []); title('LevelSetPSI');
% 		hold on;
% 		%T = zeros(size(Mask)); T(Mask) = 1;
% 		%[CCC, CH] = contour(T(:, :, Slice), [0.5 0.5]);
% 		%set(CH, 'Color', 'b');
% 		[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [GAMMA GAMMA]);
% 		set(CH, 'Color', 'g');
% 		[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [-GAMMA -GAMMA]);
% 		set(CH, 'Color', 'g');
% 		[CCC, CH] = contour(LevelSetPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'r');
% 		[CCC, CH] = contour(OriginalPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'y');
% 		[CCC, CH] = contour(OldPSI(:, :, Slice), [0 0]);
% 		set(CH, 'Color', 'm');
% 		
% 		%keyboard;
% 		%LevelSetPSI = fls_solve_eikonal_equation(LevelSetPSI, H, 'ENO1', true(size(LevelSetPSI)), 1.5);
% 	end
% 	%keyboard;
%%

% 	FV = isosurface(X, Y, Z, LevelSetPSI, 0);
% 	if(~isempty(P))
% 		delete(P);
% 	end
% 	P = patch(FV, 'EdgeColor', 'none', 'FaceColor', 'g');
% 	axis equal; lighting gouraud; light; view(3); drawnow;
% 	
	%if mod(z, 10) == 0
	%	LevelSetPSI = fls_solve_eikonal_equation(LevelSetPSI, H, 'ENO1', true(size(LevelSetPSI)), 1.5);
 	%end
% 	SubplotCols = 4;
% 	SubplotRows = 2;
% 	if mod(z, 100) == 0
% 		disp([num2str(z) ': ' num2str(TimeStep)]);
% 		AX(1) = subplot(SubplotRows, SubplotCols, 1);
% 		imshow(LevelSetPSI, []); colorbar;
% 		title('LevelSetPSI');
% 		subplot(SubplotRows, SubplotCols, 4);
% 		T = zeros(size(Mask));
% 		T(Mask) = FirstTerm;
% 		imshow(T, []); colorbar;
% 		title('First Term');
% 		AX(2) = subplot(SubplotRows, SubplotCols, 8);
% 		T = zeros(size(Mask));
% 		T(Mask) = CurvatureMinModPSI;%FirstTerm;
% 		imshow(T, []); colorbar;
% 		title('Curvature');
% 		linkaxes(AX, 'xy');
% 		subplot(SubplotRows, SubplotCols, [2 3 6 7]);
% 		hold off;
% 		imshow(D < 0, []); colorbar;
% 		title('Zero Level Set');
% 		hold on;
% 		[C, CH] = contour(LevelSetPSI, [0 0]);
% 		set(CH, 'EdgeColor', 'b');
% 		T = zeros(size(Mask));
% 		T(Mask) = GradMAGMean;
% 		subplot(SubplotRows, SubplotCols, 5);
% 		imshow(T, []); colorbar;
% 	%  	subplot(SubplotRows, SubplotCols, 4:6);
% 	%  	plot([O, LevelSetPSI(Mask)]);
% 		drawnow;
% 	%	keyboard;
% 		%NumString = num2str(z / 10);
% 		%NumString = [repmat('0', 1, 5 - length(NumString)), NumString];
% 		%print(fullfile('movie', ['movie-' NumString]), '-dpng');
% 		%keyboard;
% 	end

	%z = z + 1;
% 	if mod(z, 1000) == 0
% 		return;
% %%
% 		C = 45;
% 		clf;
% 		hold off;
% 		imshow(D(:, :, C) > 0, []);
% 		hold on;
% 		[C1, H1] = contour(OriginalPSI(:, :, C), [0 0]);
% 		set(H1, 'Color', 'r');
% 		[C2, H2] = contour(LevelSetPSI(:, :, C), [0 0]);
% 		set(H2, 'Color', 'g');
% 		
% %%
% 		FVOriginal = isosurface(OriginalPSI, 0);
% 		FVPSI = isosurface(LevelSetPSI, 0);
% %%
% 		clf;
% 		P1 = patch(FVOriginal, 'FaceColor', 'r', 'EdgeColor', 'none');
% 		P2 = patch(FVPSI, 'FaceColor', 'g', 'EdgeColor', 'none');
% 		axis equal;
% 		lighting gouraud;
% 		light;
% 		
% %%
% 	end
	%disp(num2str(z));
	%keyboard;
%%

function [ABR] = initialise_acceleration_function(D, ExclusionWeight, CurvatureWeight, H, ParameterMU)

if(ndims(D) == 3)
	CurvatureD = ls_curvature_gaussian_3d(D, H);	
else
	CurvatureD = ls_curvature(D, H, true(size(D)), 'gaussian');
end


if(any(ExclusionWeight(:)))
% 	SE = strel('disk', 11);
% 	BR = imdilate(ExclusionWeight, SE);
% 
% 	GaussianFilter = gaussian_filter_max_1d(8);
% 	BR = imfilter(double(BR), GaussianFilter, 'same', 'conv', 'replicate');
% 	BR = imfilter(BR, GaussianFilter', 'same', 'conv', 'replicate');
% 	
% 	GaussianFilter = gaussian_filter_max_1d(1);
% 	BR = imfilter(BR, reshape(GaussianFilter, [1 1 length(GaussianFilter)]), 'same', 'conv', 'replicate');
	BR = ExclusionWeight;
	M = max(BR(:));
	
	if(M == 0)
		M = 1;
	end
	
	BR = BR ./ M;
	BR = (1 - BR);
else
	BR = ones(size(D));
end

%%

%ACutoff = 1;
% A = 1 + ParameterMU .* max(CurvatureD, 0);
% 
% %A = gaussian_smooth_3d(A, 2);
% 
% Filter = gaussian_filter_max_1d(2);
% for CurDim = 1:ndims(A)
% 	IDX = ones(1, ndims(A));
% 	IDX(CurDim) = length(Filter);
% 	A = imfilter(A, reshape(Filter, IDX), 'same', 'conv', 'replicate');
% end
% clear Filter IDX;
% 
%ACutoff = 1;
%keyboard;
%if(~isempty(CurvatureWeight))
%	T = CurvatureWeight;
%else
	T = 1;
%end
%keyboard;
ABR = BR .* (1 + T .* ParameterMU .* max(CurvatureD, 0));
clear CurvatureD;

%A = gaussian_smooth_3d(A, 2);

Filter = gaussian_filter_max_1d(2);
for CurDim = 1:ndims(ABR)
	IDX = ones(1, ndims(ABR));
	IDX(CurDim) = length(Filter);
	ABR = imfilter(ABR, reshape(Filter, IDX), 'same', 'conv', 'replicate');
end
clear Filter IDX;

% use this to display the acceleration modifier
% C = 228;
% subplot 221;
% imshow(rot90(squeeze(ExclusionWeight(:, C, :))), []);
% subplot 222;
% imshow(rot90(squeeze(BR(:, C, :))), []);
% subplot 223;
% imshow(rot90(squeeze(A(:, C, :))), []);
% subplot 224;
% imshow(rot90(squeeze(ABR(:, C, :))), []);
% 

%keyboard;

%function [Mask, MaskIndices, AMask, DMask, HeavisideDMask, GradCentreAbsD, GradCentreAbsDMaskA, NumInMask, VariableImage, IdxSubs] = init_after_reinit(LevelSetPSI, D, H, ParameterMU, GAMMA)
function [Mask, MaskIndices, DMask, HeavisideDMask, GradCentreAbsDMaskA, NumInMask, VariableImage, IdxSubs, A] = init_after_reinit(LevelSetPSI, D, ExclusionWeight, CurvatureWeight, H, ParameterMU, GAMMA)

ND = ndims(LevelSetPSI);
if(ND == 2 && any(size(LevelSetPSI) == 1))
	ND = 1;
end
Mask = abs(LevelSetPSI) <= GAMMA;
MaskIndices = find(Mask);
MaskIndices = int32(MaskIndices);
NumInMask = length(MaskIndices);

VariableImage = zeros(size(Mask), 'int32');
VariableImage(Mask) = int32(1:length(MaskIndices));
IdxSubs = cell(1, ND);
%[IdxSubs{1:ND}] = ind2sub(size(Mask), MaskIndices);
[IdxSubs{:}] = ind2sub_addo(int32(size(Mask)), MaskIndices);

if ND == 2
	RealH = [H(:); 1];
else
	RealH = H(:);
end

RealHCell = num2cell(RealH);

GradCentreAbsD = cell(1, 3);
%HCell = num2cell(H(:));
if(isa(D, 'char'))
	load(D, 'FileD');
	A = initialise_acceleration_function(FileD, ExclusionWeight, CurvatureWeight, H, ParameterMU);
	FileD = abs(FileD);
	%disp(num2str(size(FileD)));
	%disp(num2str(size(Mask)));

	[GradCentreAbsD{:}] = numeric_gradient_vector_return_c(FileD, Mask, RealHCell{:});
	GradCentreAbsD = GradCentreAbsD(1:ND);
	clear FileD;
	load(D, 'FileD');
	DMask = FileD(Mask);
	clear FileD;
else	
	A = initialise_acceleration_function(D, ExclusionWeight, CurvatureWeight, H, ParameterMU);
	DMask = D(Mask);
	[GradCentreAbsD{:}] = numeric_gradient_vector_return_c(abs(D), Mask, RealHCell{:});
	GradCentreAbsD = GradCentreAbsD(1:ND);
end
AMask = A(Mask);
HeavisideDMask = 0.5 + atan(50 * DMask) / pi;

GradCentreAbsDMaskA = cell(1, ND);
%GradPlusAbsDMaskA = cell(1, ND);
for CurDim = 1:ND
	%GradMinusAbsDMaskA{CurDim} = GradMinusAbsD{CurDim}(MaskIndices) .* AMask;
%	keyboard;
	GradCentreAbsDMaskA{CurDim} = GradCentreAbsD{CurDim} .* AMask;
	%GradPlusAbsDMaskA{CurDim} = GradPlusAbsD{CurDim}(MaskIndices) .* AMask;
end


%clear A;

% function [varargout] = display_contours(AX, AVW, LevelSetPSI)
% 
% 
% DisplaySlices = [76 floor(size(AVW, 2) / 2) 39];
% 
% axes(AX(1));
% imshow(rot90(squeeze(AVW(DisplaySlices(1), :, :)), 1), []);
% hold on;
% [A, B] = contour(rot90(squeeze(LevelSetPSI(DisplaySlices(1), :, :)), 1), [0 0]);
% set(B, 'Color', 'b');
% 
% axes(AX(2));
% imshow(AVW(:, :, DisplaySlices(3)), []);
% hold on;
% 
% [A, B] = contour(LevelSetPSI(:, :, DisplaySlices(3)), [0 0]);
% set(B, 'Color', 'b');
% 
% axes(AX(3));
% imshow(rot90(squeeze(AVW(:, DisplaySlices(2), :)), 1), []);
% hold on;
% [A, B] = contour(rot90(squeeze(LevelSetPSI(:, DisplaySlices(2), :)), 1), [0 0]);
% set(B, 'Color', 'b');
