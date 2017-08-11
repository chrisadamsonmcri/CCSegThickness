function [RegisteredProfiles, BestOffsets, CostFunction] = registerSubjectProfilesToMean(Profiles, MeanProfile)

DeMeanProfile = MeanProfile - mean(MeanProfile);

if(numel(MeanProfile) ~= size(Profiles, 1))
	error('Profiles do not match in size to the mean');
end
NumNodes = numel(MeanProfile);

TemplateLeftPeak = round(NumNodes * 0.12);
TemplateRightPeak = round(NumNodes * 0.90);

leftPeakValues = 5:round(NumNodes / 4);
rightPeakValues = round(NumNodes * 3 / 4):NumNodes - 5;

CostFunction = zeros(numel(leftPeakValues), numel(rightPeakValues));

for curLeftIDX = 1:numel(leftPeakValues)
	curLeftPeak = leftPeakValues(curLeftIDX);
	for curRightIDX = 1:numel(rightPeakValues)
		curRightPeak = rightPeakValues(curRightIDX);
		R = resampleProfiles(Profiles, curLeftPeak, curRightPeak, TemplateLeftPeak, TemplateRightPeak);
		T = bsxfun(@times, bsxfun(@minus, R, mean(R)), DeMeanProfile);
		CostFunction(curLeftIDX, curRightIDX) = sum(T(:));
		
		%keyboard;
		% resample all profiles as if 
	end
end
%keyboard;
[~, I] = max(CostFunction(:));

[ID, JD] = ind2sub(size(CostFunction), I);
BestOffsets = [leftPeakValues(ID), rightPeakValues(JD)];
RegisteredProfiles = resampleProfiles(Profiles, BestOffsets(1), BestOffsets(2), TemplateLeftPeak, TemplateRightPeak);


function [R] = resampleProfiles(Profiles, LeftPeak, RightPeak, TemplateLeftPeak, TemplateRightPeak)

% hardcoded
NumNodes = size(Profiles, 1);

ReparamX = zeros(NumNodes, 1);
ReparamX(1:LeftPeak) = linspace(1, TemplateLeftPeak, LeftPeak);
ReparamX(LeftPeak:RightPeak) = linspace(TemplateLeftPeak, TemplateRightPeak, RightPeak - LeftPeak + 1);
ReparamX(RightPeak:end) = linspace(TemplateRightPeak, NumNodes, NumNodes - RightPeak + 1);

R = zeros(size(Profiles));

for z = 1:size(Profiles, 2)
	R(:, z) = interp1(ReparamX, Profiles(:, z), 1:NumNodes);
end