function [MainAX, LegAX] = cc_seg_tube_confintdiff_display(ConfIntDiff, LegendLabel, GroupLabels, CLims)

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
if(nargin < 4)
	CLims = [];
end
[V, F, CData] = cc_seg_tube_load_cdata(ConfIntDiff, 'nearest');
%keyboard;
CMAPSize = 256;

%NonSignificantColour = [0.7, 0.7, 0.7];

%CMAPX = [-1, -0.05 - eps, linspace(-0.05, -1 / CMAPSize, CMAPSize), 0, linspace(1 / CMAPSize, 0.05, CMAPSize), 0.05 + eps, 1];
%CMAP = [NonSignificantColour; NonSignificantColour; NegCMAP; [0 0 0]; PosCMAP; NonSignificantColour; NonSignificantColour];

if isempty(CLims)
	CMAPX = linspace(min(ConfIntDiff(ConfIntDiff > 0)), max(ConfIntDiff), CMAPSize);
else
	CMAPX = linspace(CLims(1), CLims(2), CMAPSize);
end
CMAP = jet(CMAPSize);
%CMAPX

[MainAX, LegAX] = cc_seg_tube_axes_display(V, F, CData, CMAPX, CMAP, GroupLabels, LegendLabel, 'PValuePlot', false);