function [PValue, ObservedTSigned, omniP, ObservedP] = cc_seg_pairedT_perm(ThicknessA, ThicknessB, NumPerms, TestType)

% [P_at_Node, omniP] = cc_seg_pairedT_perm(ThicknessA, ThicknessB, NumPerms)
%
% DESCRIPTION
%	Paired T-test between thickness profile arrays ThicknessA and ThicknessB
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
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
%

if(~ismatrix(ThicknessA) || ~isnumeric(ThicknessA) || ~ismatrix(ThicknessB) || ~isnumeric(ThicknessB))
	error('The thickness/data matrices are not matrices or not numeric');
end

if(~isequal(size(ThicknessA), size(ThicknessB)))
	error('ThicknessA and ThicknessB are different sizes');
end

if(nargin < 3)
	NumPerms = 100000;
end

if(nargin < 4)
	TestType = 'right-sided';
end

if(NumPerms < 0 || floor(NumPerms) ~= NumPerms)
	error('NumPerms needs to be a non-zero, positive integer');
end

[NumVariables, NumSubjects] = size(ThicknessA);

SQRTNumSubjects = sqrt(NumSubjects);
ObservedDifference = ThicknessA - ThicknessB;
ObservedT = mean(ObservedDifference, 2) ./ (std(ObservedDifference, 0, 2) ./ SQRTNumSubjects);
%DF = NumSubjects - 1;
% right-sided test
switch(TestType)
	case 'left-sided'
		ObservedP = tcdf(ObservedT, NumSubjects - 1);
	case 'right-sided'
		ObservedP = tcdf(-ObservedT, NumSubjects - 1);
	case 'two-sided'
		ObservedP = 2 * tcdf(-abs(ObservedT), NumSubjects - 1);
end

ObservedTSigned = ObservedT;
ObservedT = abs(ObservedT);
[ObservedTSorted, ObservedTSortedIDX] = sort(ObservedT);

if(NumVariables < 256)
	ObservedTSortedIDX = uint8(ObservedTSortedIDX);
elseif(NumVariables < 65536)
	ObservedTSortedIDX = uint16(ObservedTSortedIDX);
else
	ObservedTSortedIDX = uint32(ObservedTSortedIDX);
end

if NumPerms == 0
	PValue = [];
	omniP = [];
	return;
end
%keyboard;
% generate randomisation vectors
RandomisationSwap = (rand([NumPerms, NumSubjects], 'single') >= rand([NumPerms, NumSubjects], 'single'));
RandomisationSwap = unique(RandomisationSwap, 'rows');

% transpose so that each permutation is a column, makes it faster to extract each permutation
RandomisationSwap = RandomisationSwap';

NumActualRandomisations = size(RandomisationSwap, 2);
% compute the swaps from the previous randomisation
% the RandomisedDifference is initialised to the observed difference
RandomisedDifference = ObservedDifference;

% and each row of RandomisationSwapFromLast contains only the signs that are different from the previous randomisation
% the first row is the first randomisation verbatim since
RandomisationSwapFromLast = cat(2, RandomisationSwap(:, 1), ...
	xor(RandomisationSwap(:, 1:NumActualRandomisations - 1), RandomisationSwap(:, 2:NumActualRandomisations)));

MaxRandomisedT = zeros(NumActualRandomisations, 1);

% this algorithm is the step down maxT adjusted p-value algorithm from Westfall and Young 1993
% Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment
% see the Permutation_test_algorithm.png image in the matlab/cc_seg directory

CountRandGreater = zeros(NumVariables, 1);
OneOverNumSubjectsMinusOne = 1 ./ (NumSubjects - 1);

% perform the randomisations
for z = 1:NumActualRandomisations
	%	slow version, construct the randomisation from scratch using multiplication by signed indicator functions
	% 	SignSwap = zeros(1, NumSubjects);
	% 	SignSwap(RandomisationSwap(:, z)) = -1;
	% 	SignSwap(~RandomisationSwap(:, z)) = 1;
	% 	RandomisedDifferenceOrig = bsxfun(@times, ObservedDifference, SignSwap);
	% and only the signs that different than the previous randomisation are swapped
	RandomisedDifference(:, RandomisationSwapFromLast(:, z)) = -RandomisedDifference(:, RandomisationSwapFromLast(:, z));
	M = sum(RandomisedDifference, 2) ./ NumSubjects;
	XC = bsxfun(@minus, RandomisedDifference, M);
	S = sqrt(OneOverNumSubjectsMinusOne .* sum(XC .* XC, 2));
	
	RandomisedT = abs(M ./ (S ./ SQRTNumSubjects));
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

% the PValues are sorted according to the observed T values, reverse this ordering to restore the nodes
for z = NumVariables - 1:-1:1
	PValue(z) = max(PValue(z + 1), PValue(z));
end

%R(ObservedTSortedIDX) = 1:length(ObservedTSortedIDX);
PValue(ObservedTSortedIDX) = PValue;
%PValue = [PValue(:), ObservedTSortedIDX(:)];

omniP = sum(MaxRandomisedT >= max(ObservedT)) ./ NumActualRandomisations;
%keyboard;
% %RANDOMISATION DISTRIBUTIOObservedTSortedIDXN
% tmax=[];
% Ci=zeros(size(T));
% for count=1:size(RandStore,1)
%   RandGp1=Gp1;
%   RandGp2=Gp2;
%   RandGp1(:,find(RandStore(count,:)==0))=Gp2(:,find(RandStore(count,:)==0));
%   RandGp2(:,find(RandStore(count,:)==0))=Gp1(:,find(RandStore(count,:)==0));
%   Diff=(RandGp1-RandGp2);
%   t=mean(Diff)./std(Diff);
%   tmax=[tmax;max(max(t))];
% 
%   %STEP DOWN TEST
%   vi(1)=t(iT(1));
%   for count2=2:length(t)
%     vi(count2)=max([vi(count2-1) t(iT(count2))]);
%   end;
%   Ci=Ci+(vi>=Tsort);
% end;
% P=Ci./size(RandStore,1);
% P=fliplr(P);
% for count=2:length(P)
%   P(count)=max([P(count-1) P(count)]);
% end;
% P=fliplr(P);
% Ptrue(iT)=P;
% 
% P_at_Node = [iT' Ptrue'];
% 
% 
% %OMNIBUS TEST
% omniP=sum(tmax>=max(T))./size(RandStore,1);
