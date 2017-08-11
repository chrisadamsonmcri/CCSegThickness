function [V] = cc_seg_parc_to_values(ParcValues, ParcScheme, NumPoints)

%
% this projects values from a parcellation, Witelson, Hofer and Frahm,
% Emsell to n points at their correct locations so that 

if nargin < 3
	NumPoints = 100;
end

ParcSchemeNumValues.witelson = 5;
ParcSchemeNumValues.hoferfrahm = 5;
ParcSchemeNumValues.emsell = 8;

if(~isfield(ParcSchemeNumValues, lower(ParcScheme)))
	error(['Parcellation scheme ' ParcScheme ' not supported'])
else
	if numel(ParcValues) ~= ParcSchemeNumValues.(lower(ParcScheme))
		error([ParcScheme ' parcellation scheme requires ' ParcSchemeNumValues.(lower(ParcScheme)) ' values'])
	end
end
% 		% check the number of 
% switch(lower(ParcScheme))
% 	case 'witelson'
% 		if numel(ParcValues) ~= 5
% 			error('Witelson parcellation scheme requires 5 values')
% 		end
% 	case 'hoferfrahm'
% 		if numel(ParcValues) ~= 7
% 			error('Hofer/Frahm parcellation scheme requires 7 values')
% 		end
% 	case 'emsell'
% 		if numel(ParcValues) ~= 8
% 			error('Emsell parcellation scheme requires 8 values')
% 		end
% 	otherwise
% 		error(['Parcellation scheme ' ParcScheme ' not supported'])
% end
% 	
load ideal_callosum_interp_points_curved;
%keyboard;
InterpPoints = imresize(InterpPoints, 0.25, 'nearest');
M = InterpPoints > 0;
InterpRange = [min(InterpPoints(M)); max(InterpPoints(M))];

InterpPoints(M) = 1 - (InterpPoints(M) - InterpRange(1)) ./ (InterpRange(2) - InterpRange(1));
% %InterpPoints(M) = InterpPoints(M) .* (length(Values) - 1) + 1;
% 
%  subplot(2, 2, 1);
%  hold off;
%  imshow(InterpPoints, []);
%  hold on;
[MY, MX] = find(M);
%keyboard;
LeftEndpointX = min(MX);
LeftEndpointY = mean(MY(MX == LeftEndpointX));
RightEndpointX = max(MX);
RightEndpointY = mean(MY(MX == RightEndpointX));
% plot(LeftEndpointX, LeftEndpointY, 'r*');
% plot(RightEndpointX, RightEndpointY, 'b*');

EndpointVector = -[RightEndpointX - LeftEndpointX, RightEndpointY - LeftEndpointY];
EndpointVectorMAG = sqrt(sum(EndpointVector .* EndpointVector));
EndpointVectorNorm = EndpointVector ./ EndpointVectorMAG;

% this is the scalar resolute of the vector onto the endpoint vector
MDotEndpointVector = sum([MX - RightEndpointX, MY - LeftEndpointY] .* EndpointVectorNorm, 2);
% this is the proportion that said point is along that vector
MDotEndpointVectorProportion = MDotEndpointVector ./ EndpointVectorMAG;

% this is the equation for the vector B(1) * X + B(2)
EndpointVectorBeta = [EndpointVector(2) / EndpointVector(1), ... % slope
	LeftEndpointY - (EndpointVector(2) / EndpointVector(1) * LeftEndpointX)]; % intercept

X = [1, size(InterpPoints, 2)]';
MAboveEndpointVector = (MY < polyval(EndpointVectorBeta, MX));

Y = polyval(EndpointVectorBeta, X);
%keyboard;

% I = zeros(size(M));
% I(M) = MDotEndpointVectorProportion;
% subplot(2, 2, 2);
% imshow(I, []);
% 
% I = zeros(size(M));
% I(M) = MAboveEndpointVector;
% subplot(2, 2, 3);
% hold off;
% imshow(I, []);
% hold on;
% plot(X, Y, '*-');
%keyboard;

% now do the subdivisions
SubProportions.hoferfrahm = [1 / 6, 1 / 2, 2 / 3, 3 / 4];
SubProportions.witelson = [1 / 3, 1 / 2, 2 / 3, 4 / 5];
SubProportions.emsell = [1 / 6, 1 / 3, 1 / 2, 2 / 3, 4 / 5];

ParcIMG = zeros(size(MX), 'uint8');

CurSubProps = SubProportions.(lower(ParcScheme));

switch(lower(ParcScheme))
	case {'hoferfrahm', 'witelson'}
		ParcIMG(MDotEndpointVectorProportion < CurSubProps(1)) = 1;
		disp(num2str(max(ParcIMG)));
		for z = 2:length(CurSubProps)
			ParcIMG((MDotEndpointVectorProportion < CurSubProps(z)) & (ParcIMG == 0)) = z;
			disp(num2str(max(ParcIMG)));
		end
		ParcIMG(ParcIMG == 0) = length(CurSubProps) + 2;
		disp(num2str(max(ParcIMG)));
		ParcIMG(~MAboveEndpointVector) = 0;
	case 'emsell'
		ParcIMG(MDotEndpointVectorProportion < CurSubProps(1)) = 2;
		for z = 2:length(CurSubProps)
			ParcIMG((MDotEndpointVectorProportion < CurSubProps(z)) & (ParcIMG == 0) & MAboveEndpointVector) = z + 1;
			%disp(num2str(max(ParcIMG)));
		end
		ParcIMG(ParcIMG == 0 & MAboveEndpointVector) = length(CurSubProps) + 2;
		ParcIMG(ParcIMG == 0 & MDotEndpointVectorProportion > 0.5) = 8;
		%ParcIMG(~MAboveEndpointVector) = 9;
		PointsAntNearEndpointVector = (MAboveEndpointVector & MDotEndpointVectorProportion > 0.5 & (polyval(EndpointVectorBeta, MX)) - MY < 2);
		PointsAntNearEndpointVectorIDX = find(PointsAntNearEndpointVector);
		
		[~, AnchorPointPropIDX] = max(MX(PointsAntNearEndpointVector));
		
		AnchorPointProp = MDotEndpointVectorProportion(PointsAntNearEndpointVectorIDX(AnchorPointPropIDX));
		ParcIMG(MDotEndpointVectorProportion < AnchorPointProp & MDotEndpointVectorProportion > 1 / 2 & ~MAboveEndpointVector) = 1;
end

T = zeros(size(M), 'uint8');
T(M) = ParcIMG;
% %T(M) = PointsAntNearEndpointVector;
ParcIMG = T;
% %clear T;
% subplot(2, 2, 4);
% hold off;
% %RGB = label2rgb(ParcIMG, 'lines', 'k', 'shuffle');
% imshow(ParcIMG, []);
% hold on;
%plot(X, Y, '*-');
%keyboard;

InterpPoints(~M) = NaN;

% generate the levels for contour generation
Levels = linspace(0, 1, NumPoints + 2);
Levels = Levels(2:end - 1);

C = contourc(InterpPoints, Levels);
%[C, CC] = contour(InterpPoints, Levels);
%set(CC, 'color', 'r');
C = contour2cell(C);

ParcInterpPoints = zeros(NumPoints, 2);
for z = 1:length(C)
	C{z} = curve_resample_arc_length(C{z}, 25);
	
	ParcInterpPoints(z, :) = C{z}(12, :);
	%plot(C{z}(12, 1), C{z}(12, 2), 'r*');
end

ParcLabels = interp2(double(ParcIMG), ParcInterpPoints(:, 1), ParcInterpPoints(:, 2), 'nearest');

if(any(ParcLabels == 0))
	error('This shouldnt happen');
end
%V =
V = ParcValues(ParcLabels);
%keyboard;