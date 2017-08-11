clear;
D = dir(fullfile('datasets', 'mindboggle101', '*.nii.gz'));
FigureDir = 'mindboggle_presentation';
OutputDir = FigureDir;
D = D([1, 3:end]);
for z = 1:length(D)
	BaseName = strrep(D(z).name, '.nii.gz', '');
	CCName = strrep(D(z).name, 'msp', 'cc');
	OutputMindboggle(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
	OutputMindboggleThickness(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_thickness']));
end

D = dir(fullfile('oasis_database', 'flipped', 'OAS*_msp.nii.gz'));
FigureDir = 'oasis_presentation';
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');
%Output
%IDX = 1:length(D);
%IDX = 183;
try
	for z = 1:length(D)
		BaseName = strrep(D(z).name, '.nii.gz', '');
		CCName = strrep(D(z).name, 'msp', 'cc');
		OutputOASIS(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
		OutputOASISThickness(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_thickness']));
	end
catch e
	disp(num2str(z));
end

clf;
SR = 2;
SC = 3;
subplot(SR, SC, 1);
T = [[OutputOASIS.InitialDice]'; [OutputOASIS.FinalDice]'];
G = char(repmat('Initial Segmentation', length(OutputOASIS), 1), repmat('Final Segmentation', length(OutputOASIS), 1));
boxplot(T, G);
clear T G;
title('Dice Coefficients');

subplot(SR, SC, 2);
T = [[OutputOASIS.InitialHaussdorf]'; [OutputOASIS.FinalHaussdorf]'];
G = char(repmat('Initial Segmentation', length(OutputOASIS), 1), repmat('Final Segmentation', length(OutputOASIS), 1));

boxplot(T, G);
ylabel('Distance (mm)');
title('Hausdorff distances');
clear T G;
[~, I] = max([OutputOASIS.FinalHaussdorf]);
subplot(SR, SC, 3);
imshow(OutputOASIS(I).IMG, []);
hold on;
[~, CC] = contour(OutputOASIS(I).FinalSeg, [0.5, 0.5]);
set(CC, 'Color', 'r');
[~, CC] = contour(OutputOASIS(I).GroundSeg, [0.5, 0.5]);
set(CC, 'Color', 'g');
title({'Worst segmentation according to Hausdorff', 'Estimate (Red), Green (Ground truth)'});
subplot(SR, SC, 4);
NumNodes = 100;
GroundThicknesses = zeros(NumNodes, length(OutputOASISThickness));
for z = 1:length(OutputOASISThickness)
	GroundThicknesses(:, z) = OutputOASISThickness(z).GND.Thickness;
end

EstimatedThicknesses = cat(2, [OutputOASISThickness.Thickness]);
plot([mean(EstimatedThicknesses, 2), mean(GroundThicknesses, 2)]);
title({'Mean thickness profiles', 'Estimates (Blue), Ground Truth (Green)'});
xlabel('Node');
ylabel('Thickness (mm)');

%get correlation coefficients
subplot(SR, SC, 5);
C = zeros(1, length(OutputOASISThickness));
for z = 1:length(OutputOASISThickness)
	T = corrcoef(EstimatedThicknesses(:, z), GroundThicknesses(:, z));
	C(z) = T(1, 2);
	clear T;
end
boxplot(C);
title({'Correlation coefficients between', 'estimated and ground truth thickness profiles'});
subplot(SR, SC, 6);
[~, I] = min(C);
imshow(OutputOASIS(I).IMG, []);
hold on;
[~, CC] = contour(OutputOASIS(I).FinalSeg, [0.5, 0.5]);
set(CC, 'Color', 'r');
[~, CC] = contour(OutputOASIS(I).GroundSeg, [0.5, 0.5]);
set(CC, 'Color', 'g');
title({'Worst segmentation according to low correlation of thickness profiles', 'Estimate (Red), Green (Ground truth)'});
% W = 0.4;
% annotation('textbox', 'Position', [0.5 - W / 2, 0.9, W, 0.1], 'String', 'OASIS dataset results', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 24, 'FontWeight', 'bold');
OutputFile = 'results_oasis_overlap_thickness.png';
FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w', 'InvertHardcopy', 'off');
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
%delete(gcf);
IMG = imread(OutputFile);
IMG = imautocropwhite(IMG, 20);
imwrite(IMG, OutputFile);
		
clf;
plot(mean(cat(2, OutputMindboggleThickness.Thickness), 2));
title('Mean thickness profiles');
xlabel('Node');
ylabel('Thickness (mm)');
OutputFile = 'results_mindboggle_thickness.png';
FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w', 'InvertHardcopy', 'off');
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

IMG = imread(OutputFile);
IMG = imautocropwhite(IMG, 20);
imwrite(IMG, OutputFile);

clf;

%IDX = 1:25;
SR = 5;
SC = 5;

AX = zeros(SR, SC);

CurIDX = 1;
for CurRow = 1:SR
	for CurCol = 1:SC
		AX(CurRow, CurCol) = axes('Position', [(CurCol - 1) / SC, (CurRow - 1) / SR, 1 / SC, 1 / SR]);
		imshow(OutputMindboggle(CurIDX).IMG, []);
		hold on;
		[~, CC] = contour(OutputMindboggle(CurIDX).FinalSegArtefactsRemoved, [0.5, 0.5]);
		set(CC, 'Color', 'r');
		CurIDX = CurIDX + 1;
	end
end
OutputFile = 'results_mindboggle_seg.png';
FigPos = fullscreen_fig_pos;
set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

IMG = imread(OutputFile);
IMG = imautocropwhite(IMG, 20);
imwrite(IMG, OutputFile);

delete(gcf);
