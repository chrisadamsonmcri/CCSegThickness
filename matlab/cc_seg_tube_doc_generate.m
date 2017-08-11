clear;

P = rand([40, 1]) * 0.1;
PSign = sign(randn([40, 1]));
PSign = median_filter_1d(PSign, 5);

cc_seg_tube_p_display(P, PSign, {'CTL', 'PAT'});
%cc_seg_tube_scalar_display(P .* PSign, [], '\itp');
return;
OutputPNG = 'cc_seg_tube_example.png';
FigPos = get(gcf, 'Position');
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardcopy', 'off');
exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3) * 2, 'Height', FigPos(4) * 2, 'Color', 'rgb');

IMG = imread(OutputPNG);
IMG = 255 - IMG;
IMG = imautocropwhite(IMG, 20);
IMG = 255 - IMG;
imwrite(IMG, OutputPNG);
delete(gcf);