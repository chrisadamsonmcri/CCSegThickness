function [varargout] = cc_seg_stattest_display_results(OutputFile, PType, varargin)

%   cc_seg_stattest_display_results(OutputFile, PType, param1, val1, param2, val2, ...)
%  
%   DESCRIPTION
%     Displays the results of groupwise statistical testing from CCSegStatTest
%   PARAMETERS
%     OutputFile (string): hdf5 file where results of CCSegStatTest were saved, the following
%     fields are required
%         permP: permutation test p-values
%         observedT: T statistics from observed distribution 
%         observedP: p-values from observed distribution
%         observedPFDR: FDR corrected p-values from observed distribution
%     PType (string): determines the type of p-values to display, the
%     following types are supported:
%         'observed': the observed p-values
%         'FDR': the FDR-corrected p-values
%         'permutation': the permutation-test p-values
%   OPTIONAL KEY/VALUE PAIR PARAMETERS
%     'GroupLabels' (cell array of strings) [2]: a two-element cell array
%     with the labels of the groups. {'Group 1', 'Group 2'} by default
%	  'OutputFile' (string): the figure will be saved to this file, format
%	  will be a png file

% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
%

if nargin < 1
	error('Not enough arguments');
end

if(exist(OutputFile, 'file') ~= 2)
	error('Output file does not exist');
end

GroupLabels = {'Group 1', 'Group 2'};
PNGFile = [];

for z = 1:2:length(varargin)
	if(length(varargin) >= z + 1)
		switch(lower(varargin{z}))
			case 'grouplabels'
				GroupLabels = varargin{z + 1};
			case 'outputfile'
				PNGFile = varargin{z + 1};
		end
	end
end

if(~iscellstr(GroupLabels) || numel(GroupLabels) ~= 2)
	error('GroupLabels must be a 2-element cell array of strings');
end

[PATHSTR, NAME, EXT] = fileparts(OutputFile);

switch(lower(EXT))
    case '.hdf5'
		try
			I = h5info(OutputFile);
			[N{1:length(I.Datasets)}] = deal(I.Datasets.Name);
			for z = 1:length(N)
				S.(N{z}) = h5read(OutputFile, ['/' N{z}]);
			end
		catch e
			I = hdf5info(OutputFile);
			%keyboard;
			[N{1:length(I.GroupHierarchy.Datasets)}] = deal(I.GroupHierarchy.Datasets.Name);
			% strip the leading slashes
			N = strrep(N, '/', '');
			for z = 1:length(N)
				%disp(N{z});
				S.(N{z}) = hdf5read(OutputFile, ['/' N{z}]);
				if isa(S.(N{z}), 'hdf5.h5string')
					T = cell(1, length(S.(N{z})));
					for k = 1:length(T)
						T{k} = S.(N{z})(k).Data;
					end
					S.(N{z}) = T;
					clear T k;
				end
			end
		end
    case '.mat'
        S = load(OutputFile);
end
%keyboard;
switch(lower(PType))
	case 'observed'
		P = S.observedP;
	case 'fdr'
		P = S.observedPFDR;
	case 'permutation'
		P = S.permP;
end

P(isnan(P)) = 1;

cc_seg_tube_p_display(P, S.observedT, GroupLabels);

if ~isempty(PNGFile)
	FigPos = get(gcf, 'Position');
	set(gcf, 'PaperPosition', FigPos, 'InvertHardCopy', 'off', 'PaperUnits', 'points');

	exportfig(gcf, PNGFile, 'Format', 'png', 'Color', 'rgb', 'Width', FigPos(3), 'Height', FigPos(4));
	delete(gcf);
end