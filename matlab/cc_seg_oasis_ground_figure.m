clear;

D = dir(fullfile('oasis_database', 'flipped', 'OAS*_msp.nii.gz'));
FigureDir = 'oasis_lkccseg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

ConvexHullSubjects = [2,14,27,40,61,67,69,75,77,79,95,106,115,116,120,144,152,183,199,209,234,239];

ConvexHullUsed = false(1, length(D));

%Output
IDX = 1:length(D);

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	S(z) = load([OutputBase '_seg']);
	T(z) = load([OutputBase '_thickness']);
end

disp(['First join activated: ' num2str(sum([S.FirstJoinActivated]))]);
disp(['Second join activated: ' num2str(sum([S.SecondJoinActivated]))]);
disp(['Vein removal activated: ' num2str(sum([S.VeinRemovalActivated]))]);
%disp([': ' num2str(max([S.NumAddedByReconstruct]))]);
%disp(['Second join activated: ' num2str(sum([S.SecondJoinActivated]))]);
%return;
clf;
FigPos = [43 99 1216 637];
set(gcf, 'Position', FigPos);
SR = 1;
SC = 3;
AX = zeros(1, 3);
BH = cell(1, 2);
TextProps = {'FontSize', 16};
AX(1) = subplot(SR, SC, 1);
BH{1} = boxplot([[S.InitialDice]; [S.FinalDice]]');
ylabel('Dice coefficient', TextProps{:});
set(gca, 'XTick', [1 2], 'XTickLabel', {'Initial', 'Final'});
set(gca, TextProps{:});


AX(2) = subplot(SR, SC, 2);
BH{2} = boxplot([[S.InitialHaussdorf]; [S.FinalHaussdorf]]');
ylabel('Hausdorff distance (mm)', TextProps{:});
set(gca, 'XTick', [1 2], 'XTickLabel', {'Initial', 'Final'});
AX(3) = subplot(SR, SC, 3);
EstThickness = cat(2, T.Thickness);
GNDThickness = zeros(size(EstThickness));

IDX = [1, 20:20:100];
for z = 1:length(T)
	GNDThickness(:, z) = T(z).GND.Thickness;
end

EstGNDCorr = zeros(1, length(T));
for z = 1:length(T)
	EstGNDCorr(z) = corr(EstThickness(:, z), GNDThickness(:, z));
end

MeanEstThickness = mean(EstThickness, 2);
MeanGNDThickness = mean(GNDThickness, 2);
SEEstThickness = std(EstThickness, 0, 2) ./ sqrt(size(EstThickness, 2));
SEGNDThickness = std(GNDThickness, 0, 2) ./ sqrt(size(GNDThickness, 2));
VarEstThickness = var(EstThickness, 0, 2);
VarGNDThickness = var(GNDThickness, 0, 2);

plot([MeanEstThickness, MeanGNDThickness], 'LineWidth', 5);
hold on;
for z = 1:size(EstThickness, 1) - 1
	patch('Vertices', cat(2, [z, z, z + 1, z + 1]', ...
		[MeanEstThickness(z) - VarEstThickness(z), ...
		 MeanEstThickness(z) + VarEstThickness(z), ...
		 MeanEstThickness(z + 1) + VarEstThickness(z + 1), ...
		 MeanEstThickness(z + 1) - VarEstThickness(z + 1)]'), ...
		 'Faces', [1 2 3 4], 'FaceColor', repmat(0.75, [1 3]), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
	patch('Vertices', cat(2, [z, z, z + 1, z + 1]', ...
		[MeanGNDThickness(z) - VarGNDThickness(z), ...
		 MeanGNDThickness(z) + VarGNDThickness(z), ...
		 MeanGNDThickness(z + 1) + VarGNDThickness(z + 1), ...
		 MeanGNDThickness(z + 1) - VarGNDThickness(z + 1)]'), ...
		 'Faces', [1 2 3 4], 'FaceColor', repmat(0.5, [1 3]), 'EdgeColor', 'none', 'FaceAlpha', 0.8)
end
set(gca, 'XTick', IDX);
YL = ylim;
YL(1) = 0;
ylim(YL);
xlabel('Node', TextProps{:});
ylabel('Thickness (mm)', TextProps{:});

legend({'Estimated', 'Ground Truth'}, 'Location', 'North');

Nudge = 0.02;

for z = 1:3
	AXPos = get(AX(z), 'Position');
	AXPos(1) = AXPos(1) - Nudge * (z - 1) - 0.1 - 0.035 * double(z >= 2);
	set(AX(z), 'Position', AXPos);
	title(AX(z), roman_label(z), 'FontSize', 20);
end
AXPos = get(AX(3), 'Position');
AXPos(3) = AXPos(3) + Nudge * 5;
set(AX(3), 'Position', AXPos);
	
AX(4) = axes('Position', [AXPos(1) + AXPos(3) + 0.025, AXPos(2), 0.125, AXPos(4)]);

BH{3} = boxplot(EstGNDCorr, 'Labels', {[]});
ylabel('{\it\rho}', TextProps{:});
title(AX(4), roman_label(4), 'FontSize', 20);
set(AX, TextProps{:}, 'LineWidth', 7);

for z = 1:3
	set(BH{z}([1, 2, 3, 4, 5, 6], :), 'LineWidth', 5);
	set(BH{z}(7, :), 'MarkerSize', 20, 'LineWidth', 5);
end

FigPos = [43 99 1216 637];
set(gcf, 'Position', FigPos);
keyboard;
FigPos(3:4) = FigPos(3:4) * 2;
OutputFile = fullfile('oasis_lkccseg', 'paper_figures', 'oasis_ground_truth_comparison');
cc_seg_save_figure_paper_png(OutputFile, FigPos);
cc_seg_save_figure_paper_eps(OutputFile, FigPos);
% %cc_seg_save_figure_paper('mindboggle_ground_truth_comparison', FigPos);
% 
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w', 'InvertHardcopy', 'off');
% exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3) * 2, 'Height', FigPos(4) * 2, 'Color', 'rgb');
% 
% IMG = imread([OutputFile, '.png']);
% IMG = imautocropwhite(IMG, 10);
% imwrite(IMG, [OutputFile, '.png']);
delete(gcf);
return;
ThicknessA = cat(2, T.Thickness);
ThicknessB = zeros(size(ThicknessA));

for z = 1:length(T)
	ThicknessB(:, z) = T(z).GND.Thickness(:);
end

W = 10;
IDX = W:size(ThicknessA, 1) - W;

ThicknessC = ThicknessB(IDX, :);
ThicknessD = ThicknessC;
BB = [30, 60];
ThicknessD(BB, :) = ThicknessD(BB, :) - 3;
ThicknessD = imfilter(ThicknessD, fspecial('gaussian', [11, 1], 1), 'same', 'conv', 'replicate');
ThicknessD = ThicknessD + randn(size(ThicknessD)) * 0.2;
disp('Permutation testing');
[PValue, ObservedTSign, omniP, ObservedP] = cc_seg_pairedT_perm(ThicknessC, ThicknessD, 100000);
[~, ~, FDRP] = fdr(ObservedP, 0.05);
%%
load ideal_callosum_interp_points_curved;

P = FDRP;
M = InterpPoints > 0;
InterpRange = [min(InterpPoints(M)); max(InterpPoints(M))];

InterpPoints(M) = (InterpPoints(M) - InterpRange(1)) ./ (InterpRange(2) - InterpRange(1)) .* (length(P) - 1) + 1;

B = bwboundaries(InterpPoints > 0);
B{1} = B{1}(1:10:end, :);
I = sub2ind(size(InterpPoints), B{1}(:, 1), B{1}(:, 2));

PAtB = interp1(1:length(P), P, InterpPoints(I));

[V, F, CData] = tube_surface([B{1}(:, 2), B{1}(:, 1), zeros(size(B{1}, 1), 1)], PAtB, 50);
CMAP = spharm_p_cmap(CData);

clf;
%subplot 111;
MainAX = axes('Position', [0.2300    0.2100    0.5750    0.6150]);
PT = patch('Vertices', V, 'Faces', F, 'FaceVertexCData', CMAP, 'FaceColor', 'interp', 'EdgeColor', 'none', 'AmbientStrength', 0.8);
axis equal ij;
lighting gouraud;
view(2);
light('Position', [0 0 -1], 'Style', 'infinite')
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');

CMAPX = [linspace(0, 0.05, 256), linspace(0.05, 1, 20)];
RGBCMAP = spharm_p_cmap(CMAPX);
IMG = repmat(reshape(RGBCMAP, [size(RGBCMAP, 1), 1, 3]), 1, 20);

% RGBCMAPX = [RGBCMAPX, L];
% RGBCMAP = [RGBCMAP; repmat([0, 0, 1], 20, 1)];
% IMG = repmat(reshape(RGBCMAP, [size(RGBCMAP, 1), 1, 3]), 1, 20);
AXPos = get(gca, 'Position');
LEGAX = axes('Position', [AXPos(1) + AXPos(3) + 0.1, AXPos(2), 0.05, AXPos(4)]);
iptsetpref('imshowaxesvisible', 'on');
imshow(IMG);
axis xy;
IDX = [1:50:255, numel(CMAPX)];
set(gca, 'XTick', [], 'YTick', IDX, 'YTickLabel', cellstr(num2str(CMAPX(IDX)', '%.2f')), 'YColor', 'w');
ylabel('{\itp} value', 'Color', 'W');

%FigPos = get(gcf, 'Position');
%set(gcf, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off');
%exportfig(gcf, 'first_pass', 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
