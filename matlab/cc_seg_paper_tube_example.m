clear;

% N = 60;
% ANGLES = linspace(0, 2 * pi, N + 1);
% ANGLES = ANGLES(1:N)';
% 
% XYZ = [cos(ANGLES), sin(ANGLES), zeros(size(ANGLES))];
% 
% [V, F, CData] = tube_surface(XYZ, randn(size(XYZ, 1), 1), 0.1);
% clf;
% patch('Vertices', V, 'Faces', F, 'FaceVertexCData', CData, 'FaceColor', 'interp', 'EdgeColor', 'none');
% axis equal;
% lighting gouraud;
% light;

P = ones(40, 1);

P(25) = 0.0003;
P(26) = 0.0003;
P(24) = 0.0022;
P(27) = 0.0033;
P(20) = 0.0073;
P(23) = 0.0095;
P(21) = 0.0101;
P(1) = 0.0109;
P(19) = 0.004;
P(18) = 0.0056;
P(31) = 0.006;
P(33) = 0.006;
P(22) = 0.0063;
P(30) = 0.0068;
P(28) = 0.0068;
P(15) = 0.0073;
P(29) = 0.0078;
P(13) = 0.0106;
P(17) = 0.0125;
P(16) = 0.0125;
P(14) = 0.0135;
P(35) = 0.0152;
P(32) = 0.0165;
P(34) = 0.0216;
P(2) = 0.0281;
P(5) = 0.0347;
P(12) = 0.0349;
P(6) = 0.0352;

load ideal_callosum_interp_points_curved;

M = InterpPoints > 0;
InterpRange = [min(InterpPoints(M)); max(InterpPoints(M))];

InterpPoints(M) = (InterpPoints(M) - InterpRange(1)) ./ (InterpRange(2) - InterpRange(1)) .* (length(P) - 1) + 1;

B = bwboundaries(InterpPoints > 0);
B{1} = B{1}(1:10:end, :);
I = sub2ind(size(InterpPoints), B{1}(:, 1), B{1}(:, 2));

PAtB = interp1(1:length(P), P, InterpPoints(I));

[V, F, CData] = tube_surface([B{1}(:, 2), B{1}(:, 1), zeros(size(B{1}, 1), 1)], PAtB, 50);
CMAP = spharm_p_cmap(CData);

V(:, 2) = -V(:, 2);
V = V(:, [1 3 2]);

clf;
%subplot 111;
MainAX = axes('Position', [0.2300    0.2100    0.5750    0.6150]);
PT = patch('Vertices', V, 'Faces', F, 'FaceVertexCData', CMAP, 'FaceColor', 'interp', 'EdgeColor', 'none', 'AmbientStrength', 0.8);
axis equal ij off;
lighting gouraud;
view(-32, 10);
zoom(20);
light('Position',[0 0 -1],'Style','infinite')
set(gcf, 'Color', 'k');
%set(gca, 'Color', 'w');

CMAPX = [linspace(0, 0.05, 256), linspace(0.05, 1, 20)];
RGBCMAP = spharm_p_cmap(CMAPX);
IMG = repmat(reshape(RGBCMAP, [size(RGBCMAP, 1), 1, 3]), 1, 20);

% RGBCMAPX = [RGBCMAPX, L];
% RGBCMAP = [RGBCMAP; repmat([0, 0, 1], 20, 1)];
% IMG = repmat(reshape(RGBCMAP, [size(RGBCMAP, 1), 1, 3]), 1, 20);
H = 0.2;
AXPos = get(gca, 'Position');
LEGAX = axes('Position', [AXPos(1) + AXPos(3) + 0.01, AXPos(2) + H * AXPos(4), 0.05, (1 - 2 * H) * AXPos(4)]);
iptsetpref('imshowaxesvisible', 'on');
imshow(IMG);
axis xy;
IDX = [1:50:255, numel(CMAPX)];
set(gca, 'XTick', [], 'YTick', IDX, 'YTickLabel', cellstr(num2str(CMAPX(IDX)', '%.2f')), 'YColor', 'w', 'FontSize', 24);
ylabel('{\itp}-value', 'Color', 'w', 'FontSize', 24);

FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off');
%%
keyboard;
OutputFile = 'cc_seg_tube_example.png';
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

IMG = imread(OutputFile);
IMG = imcomplement(IMG);
%clf;
%imshow(IMG);
IMG = imautocropwhite(IMG, 10);
IMG = imcomplement(IMG);
imwrite(IMG, OutputFile, 'png');
%gzip('cc_seg_tube_example.eps');
%delete('cc_seg_tube_example.eps');
delete(gcf);