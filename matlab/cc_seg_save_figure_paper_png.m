function [varargout] = cc_seg_save_figure_paper_png(OutputFile, FigPos)

if nargin < 2
	FigPos = fullscreen_fig_pos;
end
[filepath, name] = fileparts(OutputFile);

[~, ~, ~] = mkdir(filepath);
%FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w', 'InvertHardcopy', 'off');
exportfig(gcf, fullfile(filepath, name), 'Format', 'png', 'Width', FigPos(3) * 2, 'Height', FigPos(4) * 2, 'Color', 'rgb');

IMG = imread(fullfile(filepath, [name '.png']));
IMG = imautocropwhite(IMG, 20);
imwrite(IMG, fullfile(filepath, [name '.png']));
