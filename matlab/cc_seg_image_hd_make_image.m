function [varargout] = cc_seg_image_hd_make_image(NumNodes, InputP, InputT, IDX, GroupLabels, OutputPNG)

P = ones(NumNodes, 1);
T = P;
P(IDX) = InputP;
T(IDX) = sign(InputT);
%delete(gcf);
cc_seg_tube_p_display(min(max(P, 0), 1), T, GroupLabels);
%OutputPNG = fullfile(OutputDir, ['perm_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);

FigPos = get(gcf, 'Position');
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
%delete(gcf);