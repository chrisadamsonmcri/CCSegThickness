clear;

load cc_seg_atlases;

R = [163.7248  168.6860  201.4264  125.0233];
MatchedFilter = imcrop(TemplateAVW, R);
MatchedFilterProb = imcrop(ProbAVW, R);
MatchedFilterFornixProb = imcrop(FornixProbAVW, R);

MaxProb = max(MatchedFilterProb(:));
MatchedFilterProb = MatchedFilterProb ./ MaxProb;

MaxFornixProb = max(MatchedFilterFornixProb(:));
MatchedFilterFornixProb = MatchedFilterFornixProb ./ MaxFornixProb;

FigPos =  [6 134 1563 784];
clf;
iptsetpref('ImShowAxesVisible', 'off');
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off', 'Color', 'w');
SR = 1;
SC = 3;
AX = zeros(2, 3);
AX(1, 1) = subplot(SR, SC, 1);
imshow(TemplateAVW, []);
RectProps = {'Position', R, 'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle', '--'};
rectangle(RectProps{:});
TitleProps = {'FontSize', 20};
title(roman_label(1), TitleProps{:});
AX(1, 2) = subplot(SR, SC, 2);
imshow(ProbAVW, []);
rectangle(RectProps{:});
title(roman_label(2), TitleProps{:});
AX(1, 3) = subplot(SR, SC, 3);
imshow(FornixProbAVW, []);
rectangle(RectProps{:});
title(roman_label(3), TitleProps{:});

BotAXPositions = zeros(3, 4);
AXPos = get(AX(1, 1), 'Position');
BotAXPositions(1, :) = [AXPos(1) + AXPos(3) * 3 / 5, AXPos(2) - AXPos(4) / 100, AXPos(3) / 2, AXPos(4) / 2];
AX(2, 1) = axes('Position', BotAXPositions(1, :));
imshow(MatchedFilter, []);
fix_figure_aspect(gcf, gca, size(MatchedFilter));

AXPos = get(AX(1, 2), 'Position');
BotAXPositions(2, :) = [AXPos(1) + AXPos(3) * 3 / 5, AXPos(2) - AXPos(4) / 100, AXPos(3) / 2, AXPos(4) / 2];
AX(2, 2) = axes('Position', BotAXPositions(2, :));
T = MatchedFilterProb(3:end - 2, 3:end - 2);
T = padarray(T, [1, 1], 1, 'both');
T = padarray(T, [1, 1], 0, 'both');
imshow(T, []);
fix_figure_aspect(gcf, gca, size(T));

AXPos = get(AX(1, 3), 'Position');
BotAXPositions(3, :) = [AXPos(1) + AXPos(3) * 3 / 5, AXPos(2) - AXPos(4) / 100, AXPos(3) / 2, AXPos(4) / 2];
AX(2, 3) = axes('Position', BotAXPositions(3, :));
T = MatchedFilterFornixProb(3:end - 2, 3:end - 2);
T = padarray(T, [1, 1], 1, 'both');
T = padarray(T, [1, 1], 0, 'both');
imshow(T, []);
fix_figure_aspect(gcf, gca, size(T));
AXPos = get(AX(2, 3), 'Position');

clear T;

TopAXPositions = zeros(3, 4);
for z = 1:3
	BotAXPositions(z, :) = get(AX(2, z), 'Position');
	
	fix_figure_aspect(gcf, AX(1, z), size(TemplateAVW));
	TopAXPositions(z, :) = get(AX(1, z), 'Position');
	axes(AX(2, z));
end

for z = 1:3
	TopRightCorner = [TopAXPositions(z, 1) + (R(1) + R(3)) / size(TemplateAVW, 2) * TopAXPositions(z, 3), TopAXPositions(z, 2) + TopAXPositions(z, 4) - (R(2)) / size(TemplateAVW, 1) * TopAXPositions(z, 4)];
	BottomLeftCorner = [TopAXPositions(z, 1) + R(1) / size(TemplateAVW, 2) * TopAXPositions(z, 3), TopAXPositions(z, 2) + TopAXPositions(z, 4) - (R(2) + R(4)) / size(TemplateAVW, 1) * TopAXPositions(z, 4)];
	LineProps = {'Color', repmat(0.5, 1, 3), 'LineWidth', 2, 'LineStyle', '--'};
	annotation('line', [BottomLeftCorner(1), BotAXPositions(z, 1)], [BottomLeftCorner(2), BotAXPositions(z, 2)], LineProps{:});
	annotation('line', [TopRightCorner(1), BotAXPositions(z, 1) + BotAXPositions(z, 3)], [TopRightCorner(2), BotAXPositions(z, 2) + BotAXPositions(z, 4)], LineProps{:});
	%annotation('rectangle', [BottomLeftCorner, TopRightCorner - BottomLeftCorner], 'EdgeColor', 'r');
	%annotation('rectangle', TopAXPositions(z, :), 'EdgeColor', 'g');
end

LegAX = zeros(1, 2);
LegAX(1) = axes('Position', [BotAXPositions(2, 1), BotAXPositions(2, 2) - BotAXPositions(2, 4) / 4, BotAXPositions(2, 3), BotAXPositions(2, 4) / 4.5]);
CMAPX = linspace(0, MaxProb, 256);
CMAP = repmat(linspace(0, 1, 256)', 1, 3);
CMAPIMG = repmat(reshape(CMAP, [1, size(CMAP, 1), 3]), [25, 1, 1]);
imshow(CMAPIMG, [], 'XData', CMAPX, 'YData', linspace(0, MaxProb / 10, size(CMAPIMG, 1)));
axis on;
set(gca, 'YTick', []);

LegAX(2) = axes('Position', [BotAXPositions(3, 1), BotAXPositions(3, 2) - BotAXPositions(3, 4) / 4, BotAXPositions(3, 3), BotAXPositions(3, 4) / 4.5]);
CMAPX = linspace(0, MaxFornixProb, 256);
CMAP = repmat(linspace(0, 1, 256)', 1, 3);
CMAPIMG = repmat(reshape(CMAP, [1, size(CMAP, 1), 3]), [25, 1, 1]);
imshow(CMAPIMG, [], 'XData', CMAPX, 'YData', linspace(0, MaxFornixProb / 10, size(CMAPIMG, 1)));
axis on;
set(gca, 'YTick', []);
%annotation('rectangle', AXPos, 'EdgeColor', 'r');

%imshow(MatchedFilterFornixProb, []);

TitleProps = {'FontSize', 16, 'BackgroundColor', 'w', 'Interpreter', 'latex'};
title(AX(2, 1), '$$T$$', TitleProps{:});
title(AX(2, 2), '$$PCC$$', TitleProps{:});
title(AX(2, 3), '$$PFornix$$', TitleProps{:});
OutputFile = fullfile('cc_seg_paper_figures', 'oasis_atlases');
exportfig(gcf, OutputFile, 'Color', 'rgb', 'Format', 'eps', 'Width', FigPos(3), 'Height', FigPos(4));
% IMG = imread([F '_annotated.png']);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, [F '_annotated.png']);
delete(gcf);