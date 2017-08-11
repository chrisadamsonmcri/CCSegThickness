function [MainAX, LegAX] = cc_seg_tube_p_display(P, PSign, GroupLabels)

% cc_seg_tube_p_display(P, PSign, GroupLabels)
%
% DESCRIPTION
%	Utility script for cc_seg_stattest_display_results
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%

if(nargin < 3)
	GroupLabels = {'Group 1', 'Group 2'};
end

if(~isempty(PSign))
	T = sign(PSign(:));
	T(T == 0) = 1;
	[V, F, CData] = cc_seg_tube_load_cdata(P(:) .* T, 'nearest');
else
	[V, F, CData] = cc_seg_tube_load_cdata(P(:), 'nearest');
end
%keyboard;

CMAPSize = 256;
NonSignificantColour = [0.7, 0.7, 0.7];

if(~isempty(PSign))
	NegCMAP = cool(CMAPSize);
	PosCMAP = autumn(CMAPSize);
	CMAPX = [-1, -0.05 - eps, linspace(-0.05, -1 / CMAPSize, CMAPSize), 0, linspace(1 / CMAPSize, 0.05, CMAPSize), 0.05 + eps, 1];
	CMAP = [NonSignificantColour; NonSignificantColour; NegCMAP; [0 0 0]; PosCMAP; NonSignificantColour; NonSignificantColour];
else
	CMAPX = [0, linspace(1 / CMAPSize, 0.05, CMAPSize), 0.05 + eps, 1];
	CMAP = [jet(CMAPSize + 1); NonSignificantColour; NonSignificantColour];
end

[MainAX, LegAX] = cc_seg_tube_axes_display(V, F, CData, CMAPX, CMAP, GroupLabels, '\itp',  'PValuePlot', true);