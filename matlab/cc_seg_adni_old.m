clear;

InputDir = 'CCSegInput';
D = dir(fullfile(InputDir, '*.nii.gz'));

[BaseNames{1:length(D)}] = deal(D.name);
BaseNames = strrep(BaseNames, '.nii.gz', '');
SubjectIDs = strrep(BaseNames, '_sc', '');

%FigureDir = 'cc_seg_flirt';
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

%Mode = 'seg_and_thickness';
Mode = 'skip';
%Mode = 'thickness';

IDX = 1:length(D);
% IDX = 762:length(D);
% IDX = 761;
% % flirt
% IDX = [226, 254, 340, 519, 530, 582, 592, 665, 674, 711, 756];
% IDX2 = [226, 519, 756];
% IDX = setdiff(IDX, IDX2);
% IDX = 592;
%IDX = 672;
%IDX = 1;

IDs = cell(1, length(D));
TimePoints = zeros(1, length(D));
CurSubject = 1;
%BaseName = strrep(D(z).name, '.nii.gz', '');

for z = IDX

	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	%OutputBase = fullfile(OutputDir, BaseName);
	switch(lower(Mode))
		case 'midsag'
			cc_seg_one_subject_pipe_midsag(fullfile(InputDir, D(z).name), ...
				OutputDir, ...
				BaseName, 'flirt');
			
		case 'seg_and_thickness'
		% 	if(exist([OutputBase '_seg.mat'], 'file') ~= 2)
				disp(['Subject ' num2str(CurSubject) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName], 'flirt');
				if(SegReturnCode == 0)
					ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
				end
		%	end
			if(exist([OutputBase '_seg.mat'], 'file') ~= 2 || exist([OutputBase '_thickness.mat'], 'file') ~= 2 )
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName], 'flirt');
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			end
		case 'seg'
			SegReturnCode = cc_seg_one_subject_pipe_seg(OutputBase, [], true);
			%ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
		case 'thickness'
			disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
			ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
		case 'manedit_thickness'
			if(exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			end
	end
	CurSubject = CurSubject + 1;
	%delete(gcf);
end

NumNodes = 100;

ThicknessProfiles = zeros(NumNodes, length(IDX));

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	T = load([OutputBase '_thickness.mat']);
	
	S = T.Thickness;
	S(~T.ValidStreamlines) = NaN;
	ThicknessProfiles(:, z) = S;
	clear S T;
end

FID = fopen('thickness_profiles.csv', 'w');
fprintf(FID, '%s', BaseNames{1});
fprintf(FID, ',%s', BaseNames{2:end});
fprintf(FID, '\n');

fprintf(FID, ['%f' repmat(',%f', 1, size(ThicknessProfiles, 2) - 1) '\n'], ThicknessProfiles');
fclose(FID);

% read the 
FID = fopen('ADNI_CorticalThickness_ALL_INFO.csv', 'r');
tline = fgetl(FID);
Headers = textscan(tline, '%s', 'Delimiter', ',');
clear tline;
Headers = Headers{1};

A = textscan(FID, ['%s%d%d%c%d%s' repmat('%f', 1, length(Headers) - 6)], 'Delimiter', ',', 'ReturnOnError', false);
fclose(FID);

Headers = strrep(Headers, ' ', '_');
Headers = strrep(Headers, '-', '_');

for z = 1:length(Headers)
	ADNIData.(genvarname(Headers{z})) = A{z};
end
clear A;

% sort the thickness profiles according to the locations in the ADNI data

[~, LOC] = ismember(SubjectIDs, ADNIData.PTID);

ThicknessProfiles(:, LOC) = ThicknessProfiles;
clear LOC;

[GroupIDX, GroupLabels] = grp2idx(ADNIData.dx);

NumNodes = size(ThicknessProfiles, 1);
W = 10;
IDX = W:(NumNodes - W);
%IDX
NumPerms = 100000;
NumGroups = length(GroupLabels);

% PermPValues = cell(NumGroups, NumGroups);
% PermTSignValues = cell(NumGroups, NumGroups);
% ObservedPValues = cell(NumGroups, NumGroups);
% ObservedPFDRValues = cell(NumGroups, NumGroups);
% 
% OutputDir = 'twosample_group';
% [~, ~, ~] = mkdir(OutputDir);
% 
% for GroupOne = 1:NumGroups - 1
% 	for GroupTwo = GroupOne + 1:NumGroups
% 		disp(['Permutation testing, ' GroupLabels{GroupOne} ' vs. ' GroupLabels{GroupTwo}]);
% 		ThicknessA = ThicknessProfiles(IDX, GroupIDX == GroupOne);
% 		ThicknessB = ThicknessProfiles(IDX, GroupIDX == GroupTwo);
% 		ThicknessA = ThicknessA(:, all(~isnan(ThicknessA), 1));
% 		ThicknessB = ThicknessB(:, all(~isnan(ThicknessB), 1));
% 		[PermPValues{GroupOne, GroupTwo}, ...
% 			PermTSignValues{GroupOne, GroupTwo}, ...
% 			~, ...
% 			ObservedPValues{GroupOne, GroupTwo}] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms);
% 		[~, ~, ObservedPFDRValues{GroupOne, GroupTwo}] = fdr(ObservedPValues{GroupOne, GroupTwo}, 0.05);
% 			%return;
% %		OutputPNG = fullfile(OutputDir, ['perm_group_' GroupLabels{GroupOne} '_' GroupLabels{GroupTwo} '.png']);
% %		cc_seg_image_hd_make_image(NumNodes, ...
% %			PermPValues{GroupOne, GroupTwo}, ...
% %			PermTSignValues{GroupOne, GroupTwo}, ...
% %			IDX, GroupLabels([GroupOne, GroupTwo]), OutputPNG);
% %		OutputPNG = fullfile(OutputDir, ['fdr_group_' GroupLabels{GroupOne} '_' GroupLabels{GroupTwo} '.png']);
% %		cc_seg_image_hd_make_image(NumNodes, ...
% %			ObservedPFDRValues{GroupOne, GroupTwo}, ...
% %			PermTSignValues{GroupOne, GroupTwo}, ...
% %			IDX, GroupLabels([GroupOne, GroupTwo]), OutputPNG);
% %		OutputPNG = fullfile(OutputDir, ['observed_group_' GroupLabels{GroupOne} '_' GroupLabels{GroupTwo} '.png']);
% %		cc_seg_image_hd_make_image(NumNodes, ...
% %			ObservedPValues{GroupOne, GroupTwo}, ...
% %			PermTSignValues{GroupOne, GroupTwo}, ...
% %			IDX, GroupLabels([GroupOne, GroupTwo]), OutputPNG);
% 	end
% end
% %delete(gcf);
% %keyboard;
% 
% CurI = 1;
% T = ~cellfun(@isempty, ObservedPValues);
% NumNonZero = sum(T(:));
% %keyboard;
% % ObservedP, ObservedT, PermP, ObservedPFDRValues
% 
% Headers = cell(NumNonZero, 4);
% T = cell(NumNonZero, 4);
% 
% for GroupOne = 1:NumGroups - 1
% 	for GroupTwo = GroupOne + 1:NumGroups
% 		Headers{CurI, 1} = [GroupLabels{GroupOne} 'vs' GroupLabels{GroupTwo} ' Observed P'];
% 		Headers{CurI, 2} = [GroupLabels{GroupOne} 'vs' GroupLabels{GroupTwo} ' Observed T'];
% 		Headers{CurI, 3} = [GroupLabels{GroupOne} 'vs' GroupLabels{GroupTwo} ' Perm P'];
% 		Headers{CurI, 4} = [GroupLabels{GroupOne} 'vs' GroupLabels{GroupTwo} ' FDR P'];
% 		
% 		T{CurI, 1} = ones(NumNodes, 1);
% 		T{CurI, 2} = zeros(NumNodes, 1);
% 		T{CurI, 3} = ones(NumNodes, 1);
% 		T{CurI, 4} = ones(NumNodes, 1);
% 		
% 		T{CurI, 1}(IDX) = ObservedPValues{GroupOne, GroupTwo};
% 		T{CurI, 2}(IDX) = PermTSignValues{GroupOne, GroupTwo};
% 		T{CurI, 3}(IDX) = PermPValues{GroupOne, GroupTwo};
% 		T{CurI, 4}(IDX) = ObservedPFDRValues{GroupOne, GroupTwo};
% 		CurI = CurI + 1;
% 	end
% end
% %%
% FID = fopen('twosample_groups_stats.csv', 'w');
% 
% for z = 1:numel(T)
% 	fprintf(FID, '"%s"', Headers{z});
% 	fprintf(FID, ',%f', T{z});
% 	fprintf(FID, '\n');
% end
% 
% fclose(FID);
%%
%keyboard;
% 1- We can stick with the main groups (AD, MCI, Control).
% 
% 2- These are the lobar cortical thickness measures we can use in correlations with CC:
% ThickAvg_LOBES_lh_cingulate
% ThickAvg_LOBES_lh_frontal
% ThickAvg_LOBES_lh_occipital
% ThickAvg_LOBES_lh_parietal
% ThickAvg_LOBES_lh_temporal
% ThickAvg_LOBES_rh_cingulate
% ThickAvg_LOBES_rh_frontal
% ThickAvg_LOBES_rh_occipital
% ThickAvg_LOBES_rh_parietal
% ThickAvg_LOBES_rh_temporal
% 
% 3- These are the cognitive measures we can use:
% MMSCORE_sc
% MMSEDifference_12mo_sc
% 
% 4- Yes this is baseline, cross-sectional data. We already have longitudinal measures for TBM (MEASURE_2_m12) and for ventricle volumes (VolumeDifference_12mo_sc_L and VolumeDifference_12mo_sc_R). I can look for the longitudinal freesurfer measures.
% 
% 
LobeMeasures = {'ThickAvg_LOBES_lh_cingulate', ...
	'ThickAvg_LOBES_lh_frontal', ...
	'ThickAvg_LOBES_lh_occipital', ...
	'ThickAvg_LOBES_lh_parietal', ...
	'ThickAvg_LOBES_lh_temporal', ...
	'ThickAvg_LOBES_rh_cingulate', ...
	'ThickAvg_LOBES_rh_frontal', ...
	'ThickAvg_LOBES_rh_occipital', ...
	'ThickAvg_LOBES_rh_parietal', ...
	'ThickAvg_LOBES_rh_temporal'};

GroupI = cell(length(GroupLabels), 1);

for CurGroup = 1:length(GroupLabels)
	GroupI{CurGroup} = find(GroupIDX == CurGroup);
end

NumSubjectsInGroup = cellfun('length', GroupI);
 
ANProps = {'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none', 'FontSize', 20, 'Color', 'w'};

OutputDir = 'regressions_lobar_thickness';
[~, ~, ~] = mkdir(OutputDir);

for CurMeasure = 1:length(LobeMeasures)
	RSQValues = zeros(NumNodes, length(GroupLabels));
	PValues = ones(NumNodes, length(GroupLabels));
	PFDRValues = ones(NumNodes, length(GroupLabels));
	FValues = zeros(NumNodes, length(GroupLabels));
	BValues = zeros(NumNodes, length(GroupLabels));
	
	for CurGroup = 1:length(GroupLabels)
		for CurNode = IDX
			[B,~,~,~,STATS] = regress(ADNIData.(LobeMeasures{CurMeasure})(GroupI{CurGroup}), [ThicknessProfiles(CurNode, GroupI{CurGroup})', ones(NumSubjectsInGroup(CurGroup), 1)]);
			
			%the R2 statistic, the F statistic and its p value, and an estimate of the error variance
			RSQValues(CurNode, CurGroup) = STATS(1);
			FValues(CurNode, CurGroup) = STATS(2);
			PValues(CurNode, CurGroup) = STATS(3);
			BValues(CurNode, CurGroup) = B(1);
			clear B STATS;
		end
		[~, ~, PFDRValues(IDX, CurGroup)] = fdr(PValues(IDX, CurGroup), 0.05);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(RSQValues(:, CurGroup), [], '{\itr}^2');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' LobeMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [LobeMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_rsq.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%keyboard;
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(BValues(:, CurGroup), [], '{\it\beta}');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' LobeMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [LobeMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_beta.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' LobeMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [LobeMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_p.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PFDRValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' LobeMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [LobeMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_pfdr.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		%keyboard;
	end
end

CurI = 1;
%T = ~cellfun(@isempty, PValues);
NumNonZero = size(PValues, 2);
%keyboard;
% ObservedP, ObservedT, PermP, ObservedPFDRValues

Headers = cell(NumNonZero, 5);
T = cell(NumNonZero, 5);

for GroupOne = 1:NumGroups
	Headers{CurI, 1} = [GroupLabels{GroupOne} ' BETA'];
	Headers{CurI, 2} = [GroupLabels{GroupOne} ' RSQ'];
	Headers{CurI, 3} = [GroupLabels{GroupOne} ' F'];
	Headers{CurI, 4} = [GroupLabels{GroupOne} ' P'];
	Headers{CurI, 5} = [GroupLabels{GroupOne} ' P FDR'];
		
	T{CurI, 1} = zeros(NumNodes, 1);
	T{CurI, 2} = zeros(NumNodes, 1);
	T{CurI, 3} = zeros(NumNodes, 1);
	T{CurI, 4} = ones(NumNodes, 1);
	T{CurI, 5} = ones(NumNodes, 1);
		
	T{CurI, 1} = BValues(:, GroupOne);
	T{CurI, 2} = RSQValues(:, GroupOne);
	T{CurI, 3} = FValues(:, GroupOne);
	T{CurI, 4} = PValues(:, GroupOne);
	T{CurI, 5} = PFDRValues(:, GroupOne);
	CurI = CurI + 1;
end
%%
FID = fopen('regressions_lobar_thickness.csv', 'w');

for z = 1:numel(T)
	fprintf(FID, '"%s"', Headers{z});
	fprintf(FID, ',%f', T{z});
	fprintf(FID, '\n');
end

fclose(FID);

%keyboard;
%%

CogMeasures = {'MMSCORE_sc', 'MMSEDifference_12mo_sc'};
 
OutputDir = 'regressions_cogmeasures';
[~, ~, ~] = mkdir(OutputDir);
for CurMeasure = 1:length(CogMeasures)
	RSQValues = zeros(NumNodes, length(GroupLabels));
	PValues = ones(NumNodes, length(GroupLabels));
	PFDRValues = ones(NumNodes, length(GroupLabels));
	FValues = zeros(NumNodes, length(GroupLabels));
	BValues = zeros(NumNodes, length(GroupLabels));
	
	for CurGroup = 1:length(GroupLabels)
		for CurNode = IDX
			[B,~,~,~,STATS] = regress(ADNIData.(CogMeasures{CurMeasure})(GroupI{CurGroup}), [ThicknessProfiles(CurNode, GroupI{CurGroup})', ones(NumSubjectsInGroup(CurGroup), 1)]);
			
			%the R2 statistic, the F statistic and its p value, and an estimate of the error variance
			RSQValues(CurNode, CurGroup) = STATS(1);
			FValues(CurNode, CurGroup) = STATS(2);
			PValues(CurNode, CurGroup) = STATS(3);
			BValues(CurNode, CurGroup) = B(1);
			clear B STATS;
		end
		[~, ~, PFDRValues(IDX, CurGroup)] = fdr(PValues(IDX, CurGroup), 0.05);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(RSQValues(:, CurGroup), [], '{\itr}^2');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_rsq.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(BValues(:, CurGroup), [], '{\it\beta}');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_beta.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_p.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PFDRValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_pfdr.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		%keyboard;
	end
end

CurI = 1;
%T = ~cellfun(@isempty, PValues);
NumNonZero = size(PValues, 2);
%keyboard;
% ObservedP, ObservedT, PermP, ObservedPFDRValues

Headers = cell(NumNonZero, 5);
T = cell(NumNonZero, 5);

for GroupOne = 1:NumGroups
	Headers{CurI, 1} = [GroupLabels{GroupOne} ' BETA'];
	Headers{CurI, 2} = [GroupLabels{GroupOne} ' RSQ'];
	Headers{CurI, 3} = [GroupLabels{GroupOne} ' F'];
	Headers{CurI, 4} = [GroupLabels{GroupOne} ' P'];
	Headers{CurI, 5} = [GroupLabels{GroupOne} ' P FDR'];
		
	T{CurI, 1} = zeros(NumNodes, 1);
	T{CurI, 2} = zeros(NumNodes, 1);
	T{CurI, 3} = zeros(NumNodes, 1);
	T{CurI, 4} = ones(NumNodes, 1);
	T{CurI, 5} = ones(NumNodes, 1);
		
	T{CurI, 1} = BValues(:, GroupOne);
	T{CurI, 2} = RSQValues(:, GroupOne);
	T{CurI, 3} = FValues(:, GroupOne);
	T{CurI, 4} = PValues(:, GroupOne);
	T{CurI, 5} = PFDRValues(:, GroupOne);
	CurI = CurI + 1;
end
%%
FID = fopen('regressions_cogmeasures_thickness.csv', 'w');

for z = 1:numel(T)
	fprintf(FID, '"%s"', Headers{z});
	fprintf(FID, ',%f', T{z});
	fprintf(FID, '\n');
end

fclose(FID);

CogMeasures = {'CLOCKSCOR_24moChange', ...
'COPYSCOR_24moChange', ...
'LIMMTOTAL_24moChange', ...
'AVTOT6-1_24moChange', ...
'DSPANBLTH_24moChange', ...
'CATANIMSC_24moChange', ...
'TRAASCOR_24moChange', ...
'TRABSCOR_24moChange', ...
'LDELTOTAL_24moChange', ...
'BNTTOTAL_24moChange', ...
'WAISR_Score_24mo_change', ...
'TTAU_m12Change', ...
'ABETA142_m12Change', ...
'ABETA_m24Change', ...
'TAU_m24Change', ...
'PTAU_m24Change', ...
'TAUAB_m24Change', ...
'PTAUAB_m24Change'};

CogMeasures = strrep(CogMeasures, ' ', '_');
CogMeasures = strrep(CogMeasures, '-', '_');

OutputDir = 'regressions_cogmeasures2';
[~, ~, ~] = mkdir(OutputDir);

FID = fopen('regressions_cogmeasures2_groupsizes.txt', 'w');

for CurMeasure = 1:length(CogMeasures)
	
	if(~isfield(ADNIData, CogMeasures{CurMeasure}))
		continue;
	end
	ValidSubjects = ~isnan(ADNIData.(CogMeasures{CurMeasure}));
	
	RSQValues = zeros(NumNodes, length(GroupLabels));
	PValues = ones(NumNodes, length(GroupLabels));
	PFDRValues = ones(NumNodes, length(GroupLabels));
	FValues = zeros(NumNodes, length(GroupLabels));
	BValues = zeros(NumNodes, length(GroupLabels));
	
	CurADNIData = ADNIData;
	F = fieldnames(ADNIData);
	
	for z = 1:length(F)
		CurADNIData.(F{z}) = CurADNIData.(F{z})(ValidSubjects);
	end
	
	CurThicknessProfiles = ThicknessProfiles(:, ValidSubjects);
	CurGroupIDX = GroupIDX(ValidSubjects);
	CurGroupI = GroupI;
	for z = 1:length(CurGroupI)
		CurGroupI{z} = find(CurGroupIDX == z);
	end
	
	NumSubjectsInGroup = cellfun('length', CurGroupI);
	
	T = [num2cell(NumSubjectsInGroup), GroupLabels]';
	
	fprintf(FID, [CogMeasures{CurMeasure} ': ' sprintf('%d %s, ', T{:}) '\n']);
	clear T;
	%keyboard;

	for CurGroup = 1:length(GroupLabels)
		for CurNode = IDX
			[B,~,~,~,STATS] = regress(CurADNIData.(CogMeasures{CurMeasure})(CurGroupI{CurGroup}), [ThicknessProfiles(CurNode, CurGroupI{CurGroup})', ones(NumSubjectsInGroup(CurGroup), 1)]);
			
			%the R2 statistic, the F statistic and its p value, and an estimate of the error variance
			RSQValues(CurNode, CurGroup) = STATS(1);
			FValues(CurNode, CurGroup) = STATS(2);
			PValues(CurNode, CurGroup) = STATS(3);
			BValues(CurNode, CurGroup) = B(1);
			clear B STATS;
		end
		[~, ~, PFDRValues(IDX, CurGroup)] = fdr(PValues(IDX, CurGroup), 0.05);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(RSQValues(:, CurGroup), [], '{\itr}^2');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_rsq.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_scalar_display(BValues(:, CurGroup), [], '{\it\beta}');
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_beta.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_p.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		[MainAX, LegAX] = cc_seg_tube_p_display(PFDRValues(:, CurGroup), ones(NumNodes, 1), []);
		AXPos = get(MainAX, 'Position'); 
		annotation('textbox', 'Position', [AXPos(1), AXPos(2) + AXPos(4), AXPos(3), 1 - (AXPos(2) + AXPos(4))], 'String', [GroupLabels{CurGroup} ': ' CogMeasures{CurMeasure}], ANProps{:});
		OutputPNG = fullfile(OutputDir, [CogMeasures{CurMeasure} '_' GroupLabels{CurGroup} '_pfdr.png']);
		FigPos = get(gcf, 'Position');
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
		%delete(gcf)
		IMG = imread(OutputPNG); IMG = 255 - IMG; IMG = imautocropwhite(IMG, 10); IMG = 255 - IMG; imwrite(IMG, OutputPNG);
		
		%keyboard;
	end
end

fclose(FID);

CurI = 1;
%T = ~cellfun(@isempty, PValues);
NumNonZero = size(PValues, 2);
%keyboard;
% ObservedP, ObservedT, PermP, ObservedPFDRValues

Headers = cell(NumNonZero, 5);
T = cell(NumNonZero, 5);

for GroupOne = 1:NumGroups
	Headers{CurI, 1} = [GroupLabels{GroupOne} ' BETA'];
	Headers{CurI, 2} = [GroupLabels{GroupOne} ' RSQ'];
	Headers{CurI, 3} = [GroupLabels{GroupOne} ' F'];
	Headers{CurI, 4} = [GroupLabels{GroupOne} ' P'];
	Headers{CurI, 5} = [GroupLabels{GroupOne} ' P FDR'];
		
	T{CurI, 1} = zeros(NumNodes, 1);
	T{CurI, 2} = zeros(NumNodes, 1);
	T{CurI, 3} = zeros(NumNodes, 1);
	T{CurI, 4} = ones(NumNodes, 1);
	T{CurI, 5} = ones(NumNodes, 1);
		
	T{CurI, 1} = BValues(:, GroupOne);
	T{CurI, 2} = RSQValues(:, GroupOne);
	T{CurI, 3} = FValues(:, GroupOne);
	T{CurI, 4} = PValues(:, GroupOne);
	T{CurI, 5} = PFDRValues(:, GroupOne);
	CurI = CurI + 1;
end
%%
FID = fopen('regressions_cogmeasures2_thickness.csv', 'w');

for z = 1:numel(T)
	fprintf(FID, '"%s"', Headers{z});
	fprintf(FID, ',%f', T{z});
	fprintf(FID, '\n');
end

fclose(FID);
delete(gcf);