function [varargout] = cc_seg_stattest(InputDir, OutputDir, Groups, TestType, OutputFile, varargin)

% cc_seg_stattest(InputDir, OutputDir, Groups, TestType, ..., param1, val1, param2, val2, ...)
%
% DESCRIPTION
%	Performs group-wise statistical testing of thickness profiles generated
%	by cc_seg_process
% PARAMETERS
%	InputDir (string): directory of compressed NIFTI files (*.nii.gz) for
%	processing
%	OutputDir (string): directory that contains the output of
%	cc_seg_process
%	Groups [N]: vector of group labels (1 and 2) one for each subject
%		can be a MATLAB vector or as a string to point to a file that
%		contains a newline delimited list of group labels
%	TestType (string):
%		'twosample': for Welch's 2-sample t-test (Group 1 vs. Group 2)
%		'paired': for Student's 1-sample t-test (Group 1 - Group 2)
%	OutputFile (string): mat file where results are saved, the following
%	fields are written
%		PermPValues: permutation test p-values
%		ObservedTValues: T statistics from observed distribution 
%		ObservedPValues: p-values from observed distribution
%		ObservedPFDRValues: FDR corrected p-values from observed
%		distribution
%
% PARAMETERS
%	'CullPercent' [1]: percentage of nodes to cull from the left and right
%	boundaries, defaults to 10% from either boundary
%	'NumPerms' [1]: number of permutations for testing, 100000 by default
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%


if nargin < 5
	error('Not enough arguments specified');
end

TestTypes = {'twosample', 'paired'};
if(~ismember(TestType, TestTypes))
	error(['TestType must be one of: ' deblank(sprintf('''%s'' ', TestTypes{:}))]);
end

D = dir(fullfile(InputDir, '*.nii.gz'));

if(ischar(Groups))
	if(exist(Groups, 'file') == 2)
		FID = fopen(Groups, 'r');
		GroupsReal = textscan(FID, '%d');
		GroupsReal = GroupsReal{1};
		fclose(FID);
	else
		error('Groups file does not exist');
	end
else
	GroupsReal = Groups;
end

if(numel(GroupsReal) ~= length(D) || ~isnumeric(GroupsReal))
	error('Groups error');
end

if(~all(ismember(GroupsReal(:), [1 2]')))
	error('Groups must only be 1 or 2');
end

CullPercent = 10;
NumPerms = 100000;
DisplayResults = false;

for z = 1:2:length(varargin)
	if(length(varargin) >= z + 1)
		switch(lower(varargin{z}))
			case 'cullpercent'
				CullPercent = varargin{z + 1};
			case 'numperms'
				NumPerms = varargin{z + 1};
			case 'displayresults'
				DisplayResults = varargin{z + 1};
		end
	end
end

if(~isnumeric(CullPercent) || ~isscalar(CullPercent))
	error('CullPercent should be a numeric scalar');
end
if(~isnumeric(NumPerms) || ~isscalar(NumPerms))
	error('NumPerms should be a numeric scalar');
end
if(~islogical(DisplayResults) || ~isscalar(DisplayResults))
	error('DisplayResults should be a logical scalar');
end

if(NumPerms < 1 || floor(NumPerms) ~= NumPerms)
	error('NumPerms must be a positive integer');
end

[OutputBases{1:length(D)}] = deal(D.name);
OutputBases = strrep(OutputBases, '.nii.gz', '');

for z = 1:length(D)
	OutputBases{z} = [sprintf('%04d', z) '-' OutputBases{z}];
end

% load all thickness profiles
ThicknessProfiles = cell(length(D), 1);

for z = 1:length(D)
	OutputBase = fullfile(OutputDir, OutputBases{z});
	T = load([OutputBase '_thickness'], 'Thickness');
	ThicknessProfiles{z} = T.Thickness;
	clear T OutputBase;
end

ThicknessProfiles = cat(2, ThicknessProfiles{:});

GroupNumbers = histc(GroupsReal, 1:2);

if(any(GroupNumbers < 2))
	error('Insufficient subjects in at least one group');
end

ThicknessA = ThicknessProfiles(:, GroupsReal == 1);
ThicknessB = ThicknessProfiles(:, GroupsReal == 2);

NumNodes = size(ThicknessA, 1);

CutoffNumber = round(CullPercent * NumNodes);
IDX = CutoffNumber:NumNodes - CutoffNumber;

switch(lower(TestType))
	case 'twosample'
		TestFunc = @cc_seg_twosampleT_perm;
	case 'paired'
		TestFunc = @cc_seg_pairedT_perm;
end

[PermPValues, ObservedTValues, ~, ObservedPValues] = feval(TestFunc, ThicknessA(IDX, :), ThicknessB(IDX, :), NumPerms);
[~, ~, ObservedPFDRValues] = fdr(ObservedPValues, 0.05);

NonIDX = setdiff(1:NumNodes, IDX);

T = zeros(NumNodes, 1); T(IDX) = PermPValues; T(NonIDX) = NaN; PermPValues = T;
T = zeros(NumNodes, 1); T(IDX) = ObservedTValues; T(NonIDX) = NaN; ObservedTValues = T;
T = zeros(NumNodes, 1); T(IDX) = ObservedPValues; T(NonIDX) = NaN; ObservedPValues = T;
T = zeros(NumNodes, 1); T(IDX) = ObservedPFDRValues; T(NonIDX) = NaN; ObservedPFDRValues = T;

save(OutputFile, 'PermPValues', 'ObservedTValues', 'ObservedPValues', 'ObservedPFDRValues');
