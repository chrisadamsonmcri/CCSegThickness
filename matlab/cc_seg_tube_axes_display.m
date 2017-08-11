function [MainAX, LegAX] = cc_seg_tube_axes_display(V, F, CData, CMAPX, CMAP, GroupLabels, LegendLabel, varargin)

PValuePlot = false;

for z = 1:2:length(varargin)
	if length(varargin) > z
		switch(lower(varargin{z}))
			case 'pvalueplot'
				PValuePlot = varargin{z + 1};
		end
	end
end

%keyboard;
%size(CMAP)
T = zeros(size(CData, 1), 3);
%keyboard;
if(~PValuePlot)
	
	for z = 1:3
		T(:, z) = interp1(CMAPX, CMAP(:, z), CData, 'linear');
	end
	NonSignificantColour = [0.7, 0.7, 0.7];
	
	M = any(isnan(T), 2);
	
	if(any(M))
		T(M, :) = repmat(NonSignificantColour, sum(M), 1);
	end
	%keyboard;
else
	% it is a p-value plot, we need to perform a positive and a negative
	% interpolation to get the colours right
	if(min(CMAPX) < 0)

		NegCMAPX = CMAPX(CMAPX <= 0);
		NegCMAP = CMAP(CMAPX <= 0, :);
		% replace the last colour, which is black, with the second last colour,
		% which is the lowest negative p-value colour
		NegCMAP(end, :) = NegCMAP(end - 1, :);

		PosCMAPX = CMAPX(CMAPX >= 0);
		PosCMAP = CMAP(CMAPX >= 0, :);
		% replace the first colour, which is black, with the second colour,
		% which is the lowest positive p-value colour
		PosCMAP(1, :) = PosCMAP(2, :);

		M = (CData < 0);
		for z = 1:3
			T(M, z) = interp1(NegCMAPX, NegCMAP(:, z), CData(M), 'linear');
		end
		M = (CData >= 0);
		for z = 1:3
			T(M, z) = interp1(PosCMAPX, PosCMAP(:, z), CData(M), 'linear');
		end
	else
		for z = 1:3
			T(:, z) = interp1(CMAPX, CMAP(:, z), CData, 'linear');
		end
		NonSignificantColour = [0.7, 0.7, 0.7];

		M = any(isnan(T), 2);

		if(any(M))
			T(M, :) = repmat(NonSignificantColour, sum(M), 1);
		end
	end
	%T = CData;
end
%size(T)
FigPos = [6 40 1270 878];
set(gcf, 'Position', FigPos);
clf;
MainAX = axes('Position', [0.1300 0.2100 0.5750 0.6150]);
%keyboard;
patch('Vertices', V, 'Faces', F, 'FaceVertexCData', T, 'FaceColor', 'interp', 'EdgeColor', 'none', 'AmbientStrength', 0.8);
%keyboard;
axis equal ij off;
lighting gouraud;
view(-32, 10);
set(gca, 'Clipping', 'off');
%keyboard;
VERS = version('-release');
VERS = str2double(VERS(1:4));
if VERS < 2015
	zoom(20);
else
	zoom(20);
end
%keyboard;
light('Position', [0 0 -1], 'Style', 'infinite')
set(gcf, 'Color', 'k');
%view(164, 16);  % for ISMRM

MinV = min(V);
MaxV = max(V);

% the CC is in the X-Z plane, X is AP, Z is IS
m = 256;
TextY = (MinV(2) + MaxV(2)) / 2;
TextZ = MinV(3) - (MaxV(3) - MinV(3)) / 10;

TextProps = {'Color', 'w', 'FontSize', 20, 'HorizontalAlignment', 'center'};

text(MinV(1), TextY, TextZ, 'A', TextProps{:});
text(MaxV(1), TextY, TextZ, 'P', TextProps{:});
% 

%keyboard;
%[MinV(1) + (MinV(1) + MaxV(1)) / 10, TextY, TextZ]
%[MinV(1) + 9 * (MinV(1) + MaxV(1)) / 10, TextY, TextZ]

% H = arrow([MinV(1) + (MinV(1) + MaxV(1)) / 10, TextY, TextZ], ...
% 	[MinV(1) + 9 * (MinV(1) + MaxV(1)) / 10, TextY, TextZ], ...
% 	5, ...
% 	'EdgeColor', 'w', ...
% 	'FaceColor', 'w', ...
% 	'Ends', 'both', ...
% 	'LineWidth', 2.5);
% %%
% set(H, 'AmbientStrength', 1);
% %%
% arrow fixlimits;

%LegAX = [];
% V = get(H, 'Vertices');
% hold on;
% plot3(V(:, 1), V(:, 2), V(:, 3), '*');
% text(V(:, 1), V(:, 2), V(:, 3), cellstr(num2str((1:size(V, 1))')), 'Color', 'r');
% keyboard;
%return;
IMG = repmat(reshape(CMAP, [size(CMAP, 1), 1, 3]), [1, 20, 1]);

AXPos = get(gca, 'Position');
LegAX = axes('Position', [AXPos(1) + AXPos(3) + 0.025, AXPos(2), 0.05, AXPos(4)]);
imagesc(IMG);
axis equal tight on;

% if we are displaying p-values then make the labels positive since we
% don't want -0.05, so display 0.05, 0, 0.05
if(PValuePlot)
	if(min(CMAPX) < 0)
		CMAPXDisplay = [-0.05, 0, 0.05];
		YTickIDX = interp1(CMAPX, 1:length(CMAPX), CMAPXDisplay);
		YTickLabels = cellstr(num2str(abs(CMAPXDisplay(:)), '%.2f'));
	else
		CMAPXDisplay = [0, 0.05];
		YTickIDX = interp1(CMAPX, 1:length(CMAPX), CMAPXDisplay);
		YTickLabels = cellstr(num2str(abs(CMAPXDisplay(:)), '%.2f'));
	end
else
	%CMAPX
	if all(CMAPX < 0) || all(CMAPX > 0)
		CMAPXDisplay = [min(CMAPX), max(CMAPX)];
	else
		CMAPXDisplay = unique([min(CMAPX), 0, max(CMAPX)]);
	end
	
	YTickIDX = interp1(CMAPX, 1:length(CMAPX), CMAPXDisplay);
	YTickLabels = cellstr(num2str(CMAPXDisplay(:), '%.2f'));
end
%IDX = [2  ceil(size(CMAP, 1) / 2) size(CMAP, 1) - 1];

%YTickIDX
%keyboard;
set(gca, 'XTick', [], 'YTick', YTickIDX, 'YTickLabel', YTickLabels, 'YColor', 'w', 'FontSize', 20);

% find a serif font
if(PValuePlot)
	FontName = 'Times';
else
	FontName = 'Helvetica';
end
if(~isempty(LegendLabel))
	ylabel(LegendLabel, 'Color', 'w', 'FontSize', 24, 'FontName', FontName);
end

axis xy;

TextProps = {'Color', 'w', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized'};

if(~isempty(GroupLabels))
	text(0.5, -0.1, [GroupLabels{1} ' < ' GroupLabels{2}], TextProps{:});
	text(0.5, 1.1, [GroupLabels{1} ' > ' GroupLabels{2}], TextProps{:});
end

% FigPos = fullscreen_fig_pos;
% set(gcf, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off');
% exportfig(gcf, 'cc_seg_tube_example', 'Format', 'eps', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
