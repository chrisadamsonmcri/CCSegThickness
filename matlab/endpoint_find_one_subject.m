function [Contours, SolvedImage, InitialPoints, ArcLengthsLeft, ArcLengthsRight, LeftOffsetIDXStart, CurLeftOffsetIDX, RightOffsetIDXStart, CurRightOffsetIDX] = endpoint_find_one_subject(BW, OutputFile)

% endpoint finding routine for the logical mask BW
% returns Contours which is a structure with the fields
%		xi [N, 1]: the X coordinates of the inner contour
%		yi [N, 1]: the Y coordinates of the inner contour
%		xo [M, 1]: the X coordinates of the outer contour
%		yo [M, 1]: the Y coordinates of the outer contour


% find initial guesses based on matlab's extrema in regionprops

if(nargin == 1)
	OutputFile = [];
end

[~, NumRegions] = bwlabel(BW);
if(NumRegions > 1)
	error('More than one region');
end

% cut the CC into anterior and posterior segments by separating halfway along the bounding box
R = regionprops(BW, 'BoundingBox');
MidCol = floor(R.BoundingBox(1) + R.BoundingBox(3) / 2);
BWLeft = BW;
BWRight = BW;
BWLeft(:, MidCol + 1:end) = 0;
BWRight(:, 1:MidCol) = 0;

% calculate extrema for each half
RLeft = regionprops(BWLeft, 'Extrema');
RRight = regionprops(BWRight, 'Extrema');

% the desirable extrema on the left half is bottom-right (5)
% the desirable extrema on the right half is bottom-left (6)
% see regionprops
LeftStart = RLeft.Extrema(5, :);
RightStart = RRight.Extrema(6, :);

% get a boundary contour with bwboundaries
BWBoundaryContour = bwboundaries(BW);
BWBoundaryContour = BWBoundaryContour{1};

% find the index on the boundary contour that most closely corresponds to the extrema points on the left half found previously
XC = bsxfun(@minus, LeftStart([2 1]), BWBoundaryContour);
[~, LeftStartIDX] = min(sqrt(sum(XC .* XC, 2)));

% shift the boundary contour so it starts at the left extrema point
BWBoundaryContour = circshift(BWBoundaryContour, [-(LeftStartIDX - 1), 0]);

% find the index on the boundary contour that most closely corresponds to the extrema points on the right half found previously
XC = bsxfun(@minus, RightStart([2 1]), BWBoundaryContour);
[~, RightStartIDX] = min(sqrt(sum(XC .* XC, 2)));
InitialPoints = [BWBoundaryContour(1, :); BWBoundaryContour(RightStartIDX, :)];
% find the search region prevent the left from crossing the midpoint of the CC

W = 20;
% find the first X point that exceeds the midpoint of the CC
LeftLimitPositive = find(BWBoundaryContour(:, 2) > MidCol, 1, 'first') - W - 1;
% find the last X point that exceeds the midpoint of the CC
LeftLimitNegative = -find(BWBoundaryContour([1, end:-1:2], 2) > MidCol, 1, 'first') + W + 1;

% find the limits for the right half
% move the start of the right contour to the first element
BT = circshift(BWBoundaryContour, [-(RightStartIDX - 1), 0]);
% find the first X point that is left of the midpoint of the CC
RightLimitPositive = find(BT(:, 2) < MidCol, 1, 'first') - W - 1;
% find the last X point that is left of the midpoint of the CC
RightLimitNegative = -find(BT([1, end:-1:2], 2) < MidCol, 1, 'first') + W + 1;

clear BT;

% set up vectors containing the search region
LeftOffsets = LeftLimitNegative:LeftLimitPositive;
RightOffsets = RightLimitNegative:RightLimitPositive;

% perform optimization
CurLeftOffset = 0;
CurRightOffset = 0;
CurLeftOffsetIDX = find(LeftOffsets == CurLeftOffset);
CurRightOffsetIDX = find(RightOffsets == CurRightOffset);
LeftOffsetIDXStart = find(LeftOffsets == CurLeftOffset);
RightOffsetIDXStart = find(RightOffsets == CurRightOffset);

% store the objective function values as we compute them
ArcLengthsLeft = zeros(size(LeftOffsets));
ArcLengthsRight = zeros(size(RightOffsets));

% 
SearchSize = 10;
% perform optimisation on the left offset

% algorithm
% for each direction perform gradient ascent until maxima is reached
% compute means for arc lengths on the negative and positive gradients within SearchSize
% am I at a maxima (current objective function greater than both means and greater than all arc lengths in the neighbourhood)
%	yes, terminate
%	no
%		compute gradients either side by looking at SearchSize points either side
%		move in the direction where the mean is greater
%keyboard;

% use this to precompute the closed mask for the Laplace method, for speed
MaskClosed = [];

LastFourMoves = zeros(1, 4);
while(1)
	% get the current object function
	if(ArcLengthsLeft(CurLeftOffsetIDX) == 0)
		[ArcLengthsLeft(CurLeftOffsetIDX), MaskClosed] = endpoint_find_objective_function(BWBoundaryContour, CurLeftOffset, CurRightOffset, RightStartIDX, MaskClosed);
	end
	CurObjectiveFunction = ArcLengthsLeft(CurLeftOffsetIDX);
	
	% find the mean of the objective function if we were to move the offset negative
	NegativeObjectiveMean = 0;
	NegativeIDX = [];
	if(CurLeftOffsetIDX > 1)
		NegativeIDX = CurLeftOffsetIDX - 1:-1:max(1, CurLeftOffsetIDX - SearchSize);
		for z = NegativeIDX
			if(ArcLengthsLeft(z) == 0)
				[ArcLengthsLeft(z), MaskClosed] = endpoint_find_objective_function(BWBoundaryContour, LeftOffsets(z), CurRightOffset, RightStartIDX, MaskClosed);
			end
		end
		NegativeObjectiveMean = mean(ArcLengthsLeft(NegativeIDX));
	end
	
	% find the mean of the objective function if we were to move the offset positive
	PositiveObjectiveMean = 0;
	PositiveIDX = [];
	if(CurLeftOffsetIDX < length(LeftOffsets))
		PositiveIDX = CurLeftOffsetIDX + 1:min(length(LeftOffsets), CurLeftOffsetIDX + SearchSize);
		for z = PositiveIDX
			if(ArcLengthsLeft(z) == 0)
				[ArcLengthsLeft(z), MaskClosed] = endpoint_find_objective_function(BWBoundaryContour, LeftOffsets(z), CurRightOffset, RightStartIDX, MaskClosed);
			end
		end
		PositiveObjectiveMean = mean(ArcLengthsLeft(PositiveIDX));
	end

	% am I at a maxima?
	if(~isempty(NegativeIDX))
		AllGreaterThanNegative = all(CurObjectiveFunction >= ArcLengthsLeft(NegativeIDX));
	else
		AllGreaterThanNegative = true;
	end
	if(~isempty(PositiveIDX))
		AllGreaterThanPositive = all(CurObjectiveFunction >= ArcLengthsLeft(PositiveIDX));
	else
		AllGreaterThanPositive = true;
	end
	
	if((CurObjectiveFunction > PositiveObjectiveMean && CurObjectiveFunction > NegativeObjectiveMean) || (AllGreaterThanNegative && AllGreaterThanPositive))
		%disp('At a maxima');
		break;
	else
		LastFourMoves(1:3) = LastFourMoves(2:4);
		if(PositiveObjectiveMean < NegativeObjectiveMean)
			CurLeftOffsetIDX = CurLeftOffsetIDX - 1;
			LastFourMoves(4) = -1;
			%disp('Moving negative');
		else
			CurLeftOffsetIDX = CurLeftOffsetIDX + 1;
			LastFourMoves(4) = 1;
			%disp('Moving positive');
		end
		CurLeftOffset = LeftOffsets(CurLeftOffsetIDX);
		
		% check for oscillation
		if(isequal(LastFourMoves, [-1, 1, -1, 1]) || isequal(LastFourMoves, [1, -1, 1, -1]))
			disp([mfilename '.m: Oscillating']);
			break;
		end
	end
end

% local optimisation, prone to wrong results
% while(1)
% 	% get the current object function
% 	if(ArcLengthsRight(CurRightOffsetIDX) == 0)
% 		ArcLengthsRight(CurRightOffsetIDX) = endpoint_find_objective_function(BWBoundaryContour, CurLeftOffset, CurRightOffset, RightStartIDX);
% 	end
% 	CurObjectiveFunction = ArcLengthsRight(CurRightOffsetIDX);
% 	
% 	% find the mean of the objective function if we were to move the offset negative
% 	NegativeObjectiveMean = 0;
% 	NegativeIDX = [];
% 	if(CurRightOffsetIDX > 1)
% 		NegativeIDX = CurRightOffsetIDX - 1:-1:max(1, CurRightOffsetIDX - SearchSize);
% 		for z = NegativeIDX
% 			if(ArcLengthsRight(z) == 0)
% 				ArcLengthsRight(z) = endpoint_find_objective_function(BWBoundaryContour, CurLeftOffset, RightOffsets(z), RightStartIDX);
% 			end
% 		end
% 		NegativeObjectiveMean = mean(ArcLengthsRight(NegativeIDX));
% 	end
% 	
% 	% find the mean of the objective function if we were to move the offset positive
% 	PositiveObjectiveMean = 0;
% 	PositiveIDX = [];
% 	if(CurRightOffsetIDX < length(RightOffsets))
% 		PositiveIDX = CurRightOffsetIDX + 1:min(length(RightOffsets), CurRightOffsetIDX + SearchSize);
% 		for z = PositiveIDX
% 			if(ArcLengthsRight(z) == 0)
% 				ArcLengthsRight(z) = endpoint_find_objective_function(BWBoundaryContour, CurLeftOffset, RightOffsets(z), RightStartIDX);
% 			end
% 		end
% 		PositiveObjectiveMean = mean(ArcLengthsRight(PositiveIDX));
% 	end
% 	
% 	% am I at a maxima?
% 	if(~isempty(NegativeIDX))
% 		AllGreaterThanNegative = all(CurObjectiveFunction > ArcLengthsRight(NegativeIDX));
% 	else
% 		AllGreaterThanNegative = true;
% 	end
% 	if(~isempty(PositiveIDX))
% 		AllGreaterThanPositive = all(CurObjectiveFunction > ArcLengthsRight(PositiveIDX));
% 	else
% 		AllGreaterThanPositive = true;
% 	end
% 	
% 	if((CurObjectiveFunction > PositiveObjectiveMean && CurObjectiveFunction > NegativeObjectiveMean) || (AllGreaterThanNegative && AllGreaterThanPositive))
% 		%disp('At a maxima');
% 		break;
% 	else
% 		if(PositiveObjectiveMean < NegativeObjectiveMean)
% 			CurRightOffsetIDX = CurRightOffsetIDX - 1;
% 			%disp('Moving negative');
% 		else
% 			CurRightOffsetIDX = CurRightOffsetIDX + 1;
% 			%disp('Moving positive');
% 		end
% 		CurRightOffset = RightOffsets(CurRightOffsetIDX);
% 	end
% end

% alternative method, due to the way the start point is selected, the right point always needs to go positive
MaxSoFar = 0;
for z = max(1, CurRightOffsetIDX - SearchSize):length(RightOffsets)
	if(ArcLengthsRight(z) == 0)
		[ArcLengthsRight(z), MaskClosed] = endpoint_find_objective_function(BWBoundaryContour, CurLeftOffset, RightOffsets(z), RightStartIDX, MaskClosed);
		if(ArcLengthsRight(z) < MaxSoFar / 2)
			break;
		else
			MaxSoFar = max(MaxSoFar, ArcLengthsRight(z));
		end
	end
end
% find the last maxima before the arc length starts to decrease a lot
ArcLengthsRightDone = find(ArcLengthsRight > 0);
T = ArcLengthsRight(ArcLengthsRightDone);
DilateSize = floor(SearchSize);
if(mod(DilateSize, 2) == 0)
	DilateSize = DilateSize + 1;
end

DilateMinusT = imdilate(T(:)', ones(1, DilateSize)) - T;
I = find(DilateMinusT == 0, 1, 'last');
CurRightOffsetIDX = ArcLengthsRightDone(I);
CurRightOffset = RightOffsets(CurRightOffsetIDX);

BT = circshift(BWBoundaryContour, [-CurLeftOffset, 0]);
Contours.xi = BT(1:RightStartIDX - CurLeftOffset + CurRightOffset, 2);
Contours.yi = BT(1:RightStartIDX - CurLeftOffset + CurRightOffset, 1);
Contours.xo = BT(RightStartIDX - CurLeftOffset + CurRightOffset + 1:size(BT, 1), 2);
Contours.yo = BT(RightStartIDX - CurLeftOffset + CurRightOffset + 1:size(BT, 1), 1);
%clear BT;

Contours.xo = flipud(Contours.xo);
Contours.yo = flipud(Contours.yo);
[~, ~, ~, SolvedImage] = laplace_get_points_2d_auto_mw(Contours.xi, Contours.yi, Contours.xo(end:-1:1), Contours.yo(end:-1:1), 1, 1, 1, 0, MaskClosed);

if(~isempty(OutputFile))
	clf;
	
	BoundingBox = regionprops(BW, 'BoundingBox');
	BoundingBox = BoundingBox.BoundingBox;
	
	subplot(2, 2, [1, 2]);
	T = double(BW);
	T(BWLeft) = 1;
	T(BWRight) = 2;
	imshow(imcrop(T, BoundingBox), []);
	hold on;
	
	%U = BWBoundaryContour
	I = mod(1:size(BWBoundaryContour, 1), size(BWBoundaryContour, 1)) + 1;
	U = BWBoundaryContour(I, :) - BWBoundaryContour;
	
	Q = quiver(BWBoundaryContour(:, 2) - BoundingBox(1) + 0.5, BWBoundaryContour(:, 1) - BoundingBox(2) + 0.5, U(:, 2), U(:, 1));
	
	%BWBoundaryContour
	plot([LeftStart(1), RightStart(1)] - BoundingBox(1), [LeftStart(2), RightStart(2)] - BoundingBox(2), 'r*', 'MarkerSize', 15);
	plot([BT(1, 2), BT(RightStartIDX - CurLeftOffset + CurRightOffset + 1, 2)] - BoundingBox(1), [BT(1, 1), BT(RightStartIDX - CurLeftOffset + CurRightOffset + 1, 1)] - BoundingBox(2), 'g*', 'MarkerSize', 15);
	set(gca, 'XTick', [], 'YTick', []);
	subplot(2, 2, 3);
	plot(LeftOffsets, ArcLengthsLeft);
	subplot(2, 2, 4);
	plot(RightOffsets, ArcLengthsRight);
	
	
 	FigPos = fullscreen_fig_pos;
 	set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
	%keyboard;
	cc_seg_save_figure_paper_png(OutputFile);
% 	exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 	IMG = imread(OutputFile);
% 	IMG = imautocropwhite(IMG);
% 	imwrite(IMG, OutputFile);

end

%keyboard;
