function [varargout] = cc_seg_save_figure_paper_eps(OutputFile, FigPos)

if nargin < 2
	FigPos = fullscreen_fig_pos;
end
[filepath, name] = fileparts(OutputFile);

[~, ~, ~] = mkdir(filepath);
%FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w', 'InvertHardcopy', 'off');
exportfig(gcf, fullfile(filepath, name), 'Format', 'eps', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

gzip(fullfile(filepath, [name '.eps']));
delete(fullfile(filepath, [name '.eps']));
% IMG = imread(OutputFile);
% IMG = imautocropwhite(IMG);
% imwrite(IMG, OutputFile);
