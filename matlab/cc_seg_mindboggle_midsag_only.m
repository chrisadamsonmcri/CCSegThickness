clear;

InputDir = 'RawData';
GroundTruthDir = 'ground_truth';
D = dir(fullfile(InputDir, '*.nii.gz'));
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

Mode = 'skip';

switch(lower(Mode))
	case 'seg'
		DoSeg = true;
		DoThickness = false;
		DoManEditOnly = false;
	case 'seg_and_thickness'
		DoSeg = true;
		DoThickness = true;
		DoManEditOnly = false;
	case 'thickness'
		DoSeg = false;
		DoThickness = true;
		DoManEditOnly = false;
	case 'manedit_thickness'
		DoSeg = false;
		DoThickness = true;
		DoManEditOnly = true;
	otherwise
		DoSeg = false;
		DoThickness = false;
		DoManEditOnly = false;
end

IDX = 1:length(D);
%IDX = 5;

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, BaseName);
	GroundTruthFile = fullfile(GroundTruthDir, [BaseName '_cc.nii.gz']);
	%keyboard;
	%disp(GroundTruthFile);
	if(DoSeg)
		disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
		SegReturnCode = cc_seg_one_subject_pipe_seg_figures(fullfile(InputDir, D(z).name), ...
		GroundTruthFile, ...
		OutputDir, ...
		BaseName, 'acpcdetect', false);
	end
	
	if(DoThickness)
		if(~DoManEditOnly || exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
			ThicknessReturnCode = cc_seg_one_subject_pipe_thickness_figures(OutputBase);
		end
	end
	%delete(gcf);
end

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, BaseName);
	
	S(z) = load([OutputBase '_seg']);
	T(z) = load([OutputBase '_thickness']);
end

disp(['First join activated: ' num2str(sum([S.FirstJoinActivated]))]);
disp(['Second join activated: ' num2str(sum([S.SecondJoinActivated]))]);
disp(['Vein removal activated: ' num2str(sum([S.VeinRemovalActivated]))]);
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
	AXPos(1) = AXPos(1) - Nudge * (z - 1) - 0.1 - 0.03 * double(z >= 2);
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


FigPos(3:4) = FigPos(3:4) * 2;
set(gcf, 'Position', FigPos);
%keyboard;
OutputFile = 'mindboggle_ground_truth_comparison';
%OutputFile = fullfile('oasis_lkccseg', 'paper_figures', 'mindboggle_ground_truth_comparison');

cc_seg_save_figure_paper_png(OutputFile, FigPos);
%cc_seg_save_figure_paper_eps(OutputFile, FigPos);

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

[PValue, ObservedTSign, omniP, ObservedP] = cc_seg_pairedT_perm(ThicknessA, ThicknessB, 100000);
[~, ~, FDRP] = fdr(ObservedP);


