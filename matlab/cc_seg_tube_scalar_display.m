function [MainAX, LegAX] = cc_seg_tube_scalar_display(Values, GroupLabels, LegendLabel, ValueLimits)

if(nargin < 2)
	GroupLabels = {'Group 1', 'Group 2'};
end

if(nargin < 3)
	LegendLabel = 'value';
end

if(nargin < 4)
	ValueLimits = [];
end

[V, F, CData] = cc_seg_tube_load_cdata(Values);

if(~isempty(ValueLimits))
	n = 256;
	CMAP = jet(n);
	CMAPX = linspace(ValueLimits(1), ValueLimits(2), n);
else
	
	if(all(CData(:) == CData(1)))
		CMAP = [1 1 1; 1 1 1];
		CMAPX = [CData(1) - 1; CData(1) + 1];
	else
		[CMAP, ~, CMAPX] = bluewhitered_image(517, CData);
	end
end

CMAP = squeeze(CMAP);

[MainAX, LegAX] = cc_seg_tube_axes_display(V, F, CData, CMAPX, CMAP, GroupLabels, LegendLabel, false); % last false is PValuePlot