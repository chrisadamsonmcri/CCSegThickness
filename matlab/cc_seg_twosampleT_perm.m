function [PValue, ObservedT, omniP, ObservedP] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms)

% [P_at_Node, omniP] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms)
%
% DESCRIPTION
%	Two sample T-test between thickness profile arrays ThicknessA and ThicknessB
%	Each column of ThicknessA and ThicknessB is the set of thicknesses for one callosum
%	The number of columns is the number of subjects
%	P-values are computed with a permutation test, NumPerms permutations
%
%This was made change based on Callosum_Paired.m to make it works
%Jian Chen 15/06/10
% Rewritten by Chris Adamson 2013
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%

if(~ismatrix(ThicknessA) || ~isnumeric(ThicknessA) || ~ismatrix(ThicknessB) || ~isnumeric(ThicknessB))
	error('The thickness matrices are not matrices or not numeric');
end

if(~isequal(size(ThicknessA, 1), size(ThicknessB, 1)))
	error('ThicknessA and ThicknessB are different sizes');
end

if(nargin < 3)
	NumPerms = 100000;
end

if(NumPerms < 1 || floor(NumPerms) ~= NumPerms)
	error('NumPerms needs to be a non-zero, positive integer');
end

ThicknessAll = cat(2, ThicknessA, ThicknessB);

[~, NumSubjectsA] = size(ThicknessA);
[NumVariables, NumSubjectsB] = size(ThicknessB);

[ObservedT, ObservedP] = cc_seg_twosampleT_perm_calcT(ThicknessA, ThicknessB);

ObservedTSign = sign(ObservedT);
ObservedT = abs(ObservedT);
%[ObservedP, ObservedT]

[ObservedTSorted, ObservedTSortedIDX] = sort(ObservedT);

NumSubjectsTotal = NumSubjectsA + NumSubjectsB;

%keyboard;
% generate randomisation vectors
%[~, RandomisationIDX] = sort(rand([NumPerms, NumSubjectsTotal], 'single'), 2);
% sorting the columns takes half the time
[~, RandomisationIDX] = sort(rand([NumSubjectsTotal, NumPerms], 'single'));

if(NumSubjectsTotal < 256)
	RandomisationIDX = uint8(RandomisationIDX);
elseif(NumSubjectsTotal < 65536)
	RandomisationIDX = uint16(RandomisationIDX);
else
	RandomisationIDX = uint32(RandomisationIDX);
end

RandomisationIDX = RandomisationIDX';

GroupAIDX = 1:NumSubjectsA;
GroupBIDX = NumSubjectsA + 1:NumSubjectsTotal;

RandomisationIDX(:, GroupAIDX) = sort(RandomisationIDX(:, GroupAIDX), 2);
RandomisationIDX(:, GroupBIDX) = sort(RandomisationIDX(:, GroupBIDX), 2);
RandomisationIDX = unique(RandomisationIDX, 'rows');

%RandomisationIDX = setdiff(RandomisationIDX, 1:NumSubjectsTotal, 'rows');

% take out the randomisation that doesn't swap any subjects,
% will always be the first row
if(all(RandomisationIDX(1, :) == (1:NumSubjectsTotal)))
	RandomisationIDX = RandomisationIDX(2:end, :);
end

% transpose so that each permutation is a column, makes it faster to extract each permutation
RandomisationIDX = RandomisationIDX';

NumActualRandomisations = size(RandomisationIDX, 2);
MaxRandomisedT = zeros(NumActualRandomisations, 1);

% this algorithm is the step down maxT adjusted p-value algorithm from Westfall and Young 1993
% Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment
% see the Permutation_test_algorithm.png image in the matlab/cc_seg directory

CountRandGreater = zeros(NumVariables, 1);
%OneOverNumSubjectsMinusOne = 1 ./ (NumSubjects - 1);
% perform the randomisations
for z = 1:NumActualRandomisations
	
	RandomisedThicknessA = ThicknessAll(:, RandomisationIDX(GroupAIDX, z));
	RandomisedThicknessB = ThicknessAll(:, RandomisationIDX(GroupBIDX, z));
	RandomisedT = cc_seg_twosampleT_perm_calcT(RandomisedThicknessA, RandomisedThicknessB);
	RandomisedT = abs(RandomisedT);
	MaxRandomisedT(z) = max(RandomisedT);
	
	% step down test
	% http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_multtest_sect014.htm
	% pesudocode
	% go down the order of the observed T, going from the 
	% StepDownP(1) is the randomised T value of the lowest T value from the observed
	% then 
	% StepDownP(i) is the maximum of StepDownP(i - 1) and the randomised T value of the (i'th) highest T value from the observed
	% 
	% so if you observe a high randomised T for a low observed T, then all the other T values are given this T value
	% so the randomised T values are boosted if any previous randomised T has exceeded the observed T until the randomised T is greater
	% than the last randomised T
	
	%StepDownT = zeros(NumVariables, 1);
% 	StepDownT = RandomisedT(ObservedTSortedIDX);
% 	for CurT = NumVariables - 1:-1:1
% 		StepDownT(CurT) = max(StepDownT(CurT + 1), StepDownT(CurT));
% 	end
	StepDownT = cummax(RandomisedT(ObservedTSortedIDX));
	
	% then the we add one to the count if the Randomised T value is greater than the observed
	%keyboard;
	CountRandGreater = CountRandGreater + double(StepDownT >= ObservedTSorted);
	%clear StepDownT CurT SignSwap M XC S;
end
% so the PValue is the proportion of times that the randomised T value was greater
PValue = CountRandGreater ./ NumActualRandomisations;

%T = PValue;
% the PValues are sorted according to the observed T values, reverse this ordering to restore the nodes
for z = NumVariables - 1:-1:1
	%PValue(z) = max(PValue(z + 1), PValue(z));
	if PValue(z) < PValue(z + 1)
		PValue(z) = PValue(z + 1);
	end
end

R(ObservedTSortedIDX) = 1:length(ObservedTSortedIDX);
PValue(ObservedTSortedIDX) = PValue;
%PValue = [PValue(:), ObservedTSortedIDX(:)];

omniP = sum(MaxRandomisedT >= max(ObservedT)) ./ NumActualRandomisations;
ObservedT = ObservedT .* ObservedTSign;

function [T, P] = cc_seg_twosampleT_perm_calcT(ThicknessA, ThicknessB)

NA = size(ThicknessA, 2);
NB = size(ThicknessB, 2);

MeanA = sum(ThicknessA, 2) ./ NA;
MeanB = sum(ThicknessB, 2) ./ NB;

AXC = bsxfun(@minus, ThicknessA, MeanA);
BXC = bsxfun(@minus, ThicknessB, MeanB);

VarA = sum(AXC .* AXC, 2) ./ (NA - 1);
VarB = sum(BXC .* BXC, 2) ./ (NB - 1);

Den = sqrt(VarA ./ NA + VarB ./ NB);

T = (MeanA - MeanB) ./ Den;

if nargout > 1
	% calculate degrees of freedom, Welch's t-test, https://en.wikipedia.org/wiki/Welch%27s_t_test
	DF = ((VarA ./ NA + VarB ./ NB) .* (VarA ./ NA + VarB ./ NB)) ./ ((VarA .* VarA ./ NA ./ NA) ./ (NA - 1) + (VarB .* VarB ./ NB ./ NB) ./ (NB - 1));
	%DF = Numerator ./ Den;
	% two-tailed
	P = 2 * tcdf(-abs(T), DF);
	% right-tailed
	%P = tcdf(-T, DF);
	% left-tailed
	%P = tcdf(T, DF);
	
end
