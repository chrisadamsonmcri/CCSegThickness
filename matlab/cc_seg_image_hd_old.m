clear;

InputDir = 'ConvertedNII';
D = dir(fullfile(InputDir, '*.nii.gz'));
%FigureDir = 'cc_seg_flirt';
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

%Mode = 'seg_and_thickness';
Mode = 'skip';
IDX = 1:length(D);
%IDX = 173;
IDs = cell(1, length(D));
TimePoints = zeros(1, length(D));

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	switch(lower(Mode))
		case 'seg_and_thickness'
		% 	if(exist([OutputBase '_seg.mat'], 'file') ~= 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
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
		case 'thickness'
			disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
			ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
		case 'manedit_thickness'
			if(exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			end
	end
	
	[tokens] = regexp(BaseName, '^\d+-\d+-(\d+-\d+)-\d+-\d+-\d+-SRT([123])$', 'tokens');
	%keyboard;
	IDs{z} = tokens{1}{1};
	TimePoints(z) = str2double(tokens{1}{2});
	clear tokens;
	%Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
	
end
%delete(gcf);
for z = 1:length(IDs)
	if(length(IDs{z}) == 4)
 		IDs{z} = [IDs{z}(1:3) '0' IDs{z}(4)];
	end
end

[UniqueIDs, ~, J] = unique(IDs);
UniqueIDs = strrep(UniqueIDs, '-', '.');

ThicknessProfiles = cell(length(UniqueIDs), 3);

for z = 1:length(UniqueIDs)
	CurIDX = find(J == z);
	CurTimePoints = TimePoints(CurIDX);
	[SortedCurTimePoints, I] = sort(CurTimePoints);
	CurIDX = CurIDX(I);
	for k = 1:length(SortedCurTimePoints)
		BaseName = strrep(D(CurIDX(k)).name, '.nii.gz', '');
		OutputBase = fullfile(OutputDir, [sprintf('%04d', CurIDX(k)) '-' BaseName]);
		T = load([OutputBase '_thickness'], 'Thickness');
		ThicknessProfiles{z, SortedCurTimePoints(k)} = T.Thickness;
		clear T;
	end
end

ISE = cellfun(@isempty, ThicknessProfiles);
I = find(~ISE(:), 1, 'first');
SZ = size(ThicknessProfiles{I});
[ThicknessProfiles{ISE}] = deal(NaN(SZ));

ThicknessProfiles = reshape(ThicknessProfiles, [1, size(ThicknessProfiles, 1), 3]);
ThicknessProfiles = cell2mat(ThicknessProfiles);
%keyboard;
% read in the groups
% FID = fopen('image-hd-ids-groups.csv', 'r');
% A = textscan(FID, '%s%d', 'Delimiter', ',');
% GroupFileIDs = A{1};
% GroupFileGroups = A{2};
% 
% %clear A;
% fclose(FID);
% % pad the IDs that are XX.Y as XX.0Y
for z = 1:length(UniqueIDs)
	if(length(UniqueIDs{z}) == 4)
 		UniqueIDs{z} = [UniqueIDs{z}(1:3) '0' UniqueIDs{z}(4)];
	end
end

FID = fopen('image-hd-ids-groups-all.csv', 'r');
LineOne = fgetl(FID);
LineTwo = fgetl(FID);
LineOneHeaders = textscan(LineOne, '%s', 'Delimiter', ',');
LineOneHeaders = LineOneHeaders{1};
for z = 2:length(LineOneHeaders)
	if(isempty(LineOneHeaders{z}))
		LineOneHeaders{z} = LineOneHeaders{z - 1};
	end
end
LineOneHeaders{length(LineOneHeaders) + 1} = LineOneHeaders{length(LineOneHeaders)};
LineTwoHeaders = textscan(LineTwo, '%s', 'Delimiter', ',');
LineTwoHeaders = LineTwoHeaders{1};

LineOneHeaders = strrep(LineOneHeaders, '#', '');
LineOneHeaders = strrep(LineOneHeaders, ' ', '_');

F = strfind(LineTwoHeaders, '(');
%F = regexp(LineTwoHeaders, '^([^\(]+)\(?', 'tokens');

for z = 1:length(F)
	if(~isempty(F{z}))
		LineTwoHeaders{z} = LineTwoHeaders{z}(1:F{z}(1) - 1);
	end
end

clear F;
LineTwoHeaders = deblank(LineTwoHeaders);
LineTwoHeaders = strrep(LineTwoHeaders, ' ', '_');

FormatString = ['%d%s%d%s%s%s' repmat('%f', 1, length(LineOneHeaders) - 6)];

A = textscan(FID, FormatString, 'Delimiter', ',');

fclose(FID);
%clear LineOne LineTwo;

for z = 1:length(LineOneHeaders)
	if(isempty(LineTwoHeaders{z}))
		NeuropsychData.(LineOneHeaders{z}) = A{z};
	else
		NeuropsychData.(LineOneHeaders{z}).(LineTwoHeaders{z}) = A{z};
	end
end

clear A;

for z = 1:length(NeuropsychData.ID)
	if(length(NeuropsychData.ID{z}) == 4)
 		NeuropsychData.ID{z} = [NeuropsychData.ID{z} '0'];
	end
end

%keyboard;
GroupNames = {'CTL', 'PRE', 'SYM'};

% assign the groups to the IDs
TF = ismember(UniqueIDs, NeuropsychData.ID);

ThicknessProfiles = ThicknessProfiles(:, TF, :);
UniqueIDs = UniqueIDs(TF);

TF = ismember(NeuropsychData.ID, UniqueIDs);

%LOC = LOC(LOC > 0);
%Groups = NeuropsychData.Group(LOC);

TopLevelFields = fieldnames(NeuropsychData);

for TopLevelFieldIDX = 1:length(TopLevelFields)
	if(~isstruct(NeuropsychData.(TopLevelFields{TopLevelFieldIDX})))
		NeuropsychData.(TopLevelFields{TopLevelFieldIDX}) = NeuropsychData.(TopLevelFields{TopLevelFieldIDX})(TF);
	else
		SecondLevelFields = fieldnames(NeuropsychData.(TopLevelFields{TopLevelFieldIDX}));
		for SecondLevelFieldIDX = 1:length(SecondLevelFields)
			NeuropsychData.(TopLevelFields{TopLevelFieldIDX}).(SecondLevelFields{SecondLevelFieldIDX}) = NeuropsychData.(TopLevelFields{TopLevelFieldIDX}).(SecondLevelFields{SecondLevelFieldIDX})(TF);
		end
	end
end
%keyboard;

% so instead of having all the data in a vector, organise into a
% [Subject, TimePoint] matrix

[UniqueNeuropsychIDs, ~, J] = unique(NeuropsychData.ID);
NeuropsychData.ID = UniqueNeuropsychIDs;

NumTimePoints = 3;
% make linear indices for the resultant arrays
SZ = [length(UniqueNeuropsychIDs), NumTimePoints];

SubI = sub2ind(SZ, J, double(NeuropsychData.TimePoint));
clear J UniqueNeuropsychIDs;

TopLevelFields = setdiff(fieldnames(NeuropsychData), {'TimePoint', 'ID'});

for TopLevelFieldIDX = 1:length(TopLevelFields)
	if(~isstruct(NeuropsychData.(TopLevelFields{TopLevelFieldIDX})))
		if(isnumeric(NeuropsychData.(TopLevelFields{TopLevelFieldIDX})))
			T = NaN(SZ);
		else
			T = cell(SZ);
		end
		T(SubI) = NeuropsychData.(TopLevelFields{TopLevelFieldIDX});
		NeuropsychData.(TopLevelFields{TopLevelFieldIDX}) = T;
		clear T;
	else
		SecondLevelFields = fieldnames(NeuropsychData.(TopLevelFields{TopLevelFieldIDX}));
		for SecondLevelFieldIDX = 1:length(SecondLevelFields)
			if(isnumeric(NeuropsychData.(TopLevelFields{TopLevelFieldIDX}).(SecondLevelFields{SecondLevelFieldIDX})))
				T = NaN(SZ);
			else
				T = cell(SZ);
			end
			T(SubI) = NeuropsychData.(TopLevelFields{TopLevelFieldIDX}).(SecondLevelFields{SecondLevelFieldIDX});
			NeuropsychData.(TopLevelFields{TopLevelFieldIDX}).(SecondLevelFields{SecondLevelFieldIDX}) = T;
			clear T;
		end
	end
end

% now deal with group and gender which is special, only one point per time
% point is required
% remove the NaN values from the group labels by using max (ignores NaN)

NeuropsychData.Group = max(NeuropsychData.Group, [], 2);

ISE = ~cellfun(@isempty, NeuropsychData.Gender);

[~, J] = max(ISE, [], 2);
%[I, J] = find(ISE);

%keyboard;
SubI = sub2ind(size(ISE), (1:size(ISE, 1))', J);
NeuropsychData.Gender = NeuropsychData.Gender(SubI);
NeuropsychData.Gender(strcmp(NeuropsychData.Gender, '1')) = {'M'};
NeuropsychData.Gender(strcmp(NeuropsychData.Gender, '2')) = {'F'};
clear J SubI ISE;
NeuropsychData.Gender = char(upper(NeuropsychData.Gender));

NeuropsychData = rmfield(NeuropsychData, 'TimePoint');

% correction for the third timepoint
NeuropsychData.Brain_Volumes.total_Intracranial_vol(:, 3) = NeuropsychData.Brain_Volumes.total_Intracranial_vol(:, 3) / 1000;

NumGroups = length(unique(NeuropsychData.Group));
%GroupFileGroups = GroupFileGroups(TF);
%clear TF LOC k z;
%keyboard;

[NumNodes, NumSubjects, NumTimePoints] = size(ThicknessProfiles);

% compare each group's mean profile over time

% AX = zeros(1, 3);
% SurfHandles = zeros(1, 3);
% ThicknessDiffs = cell(1, 3);
% clf;
% for CurGroup = 1:3
% 	T = zeros(size(ThicknessProfiles, 1), 3);
% 	for CurTimePoint = 1:3
% 		T(:, CurTimePoint) = nanmean(ThicknessProfiles(:, Groups == CurGroup, CurTimePoint), 2);
% 	end
% 	%ThicknessDiffs{CurGroup} = cat(2, zeros(size(ThicknessProfiles, 1), 1), diff(T, 1, 2));
% 	ThicknessDiffs{CurGroup} = cat(2, zeros(size(ThicknessProfiles, 1), 1), T(:, 2) - T(:, 1), zeros(size(ThicknessProfiles, 1), 1), T(:, 3) - T(:, 2), zeros(size(ThicknessProfiles, 1), 1));
% 	%AX(CurGroup) = subplot(2, 2, CurGroup);
% 	AX(CurGroup) = axes('Position', [(CurGroup - 1) / 3, 0, 1 / 3, 1]);
% 	%SurfHandles(CurGroup) = surf(T);
% 	
% 	SurfHandles(CurGroup) = surf(repmat(1:0.5:3, size(T, 1), 1), repmat((1:size(ThicknessProfiles, 1))', 1, size(ThicknessDiffs{CurGroup}, 2)), [T(:, 1), (T(:, 1) + T(:, 2)) / 2, T(:, 2), (T(:, 2) + T(:, 3)) / 2, T(:, 3)]);
% 	%set(SurfHandles(z), 'CData', D, 'FaceColor', 'interp');%keyboard;
% 	axis equal;
% 	view(-23, 28);
% 	%colorbar;
% 	%title(['Group ' num2str(CurGroup)]);
% 	%get(gca, 'DataAspect');
% 	%set(gca, 'DataAspect', [2 1 1]);
% end
% 
% [newmap, newimg, cmapx] = bluewhitered_image(256, cat(2, ThicknessDiffs{:}));
% %keyboard;
% newimg = mat2cell(newimg, size(ThicknessProfiles, 1), repmat(size(ThicknessDiffs{1}, 2), 1, length(ThicknessDiffs)), 3);
% for z = 1:3
% 	set(SurfHandles(z), 'CData', newimg{z}, 'FaceColor', 'interp');
% end
% %%
% set(AX, 'DataAspect', [0.3 1 1]);
% %%
% Nudge = 0.05;
% for z = 1:3
% 	AXPos = get(AX(z), 'Position');
% 	AXPos(1) = AXPos(1) + Nudge;
% 	set(AX(z), 'Position', AXPos);
% 	text(0.01, 1.05, GroupNames{z}, 'Units', 'normalized', 'Parent', AX(z));
% end
% 
% Nudge = 0.1;
% 
% for z = 2:3
% 	AXPos = get(AX(z), 'Position');
% 	AXPos(1) = AXPos(1) - Nudge * (z - 1);
% 	set(AX(z), 'Position', AXPos);
% end
% 
% AXPos = get(AX(3), 'Position');
% LegAX = axes('Position', [AXPos(1) + AXPos(3) / 1.2, AXPos(2) + AXPos(4) / 1.9, 0.025, AXPos(4) / 5]);
% C = repmat(reshape(newmap, [256, 1, 3]), [1, 10, 1]);
% imagesc([0, 1], [cmapx(1), cmapx(end)], C);
% axis xy;
% set(gca, 'XTick', []);
% ylabel({'Difference in thickness (mm)', 'from previous timepoint'});
% %keyboard;
% FigPos = fullscreen_fig_pos;
% OutputPNG = 'mean_thicknesses_by_group_surf.png';
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'w');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% delete(gcf);
%set(gcf, 'Color', 'w');


% Chris, the other variables they are keen on (as per the email chain below) are:
% 
% ITIPTAP-slow and �fast
% UPSIT
% SDMT
% Stroop
% CAG repeats
% Disease burden score
% 
% They are in the attached database.
% 
% Nellie below suggests we just do correlations. The first four measures are in all groups;
% the last two just in groups 2 & 3. I think for all variables, we could run them in groups
% 2 & 3 combined, and then group 1 only and 3 only. What do you think?

%M = isnan(NeuropsychData.Brain_Volumes.total_Intracranial_vol);
%NeuropsychData.Brain_Volumes.total_Intracranial_vol(M) = mean(NeuropsychData.Brain_Volumes.total_Intracranial_vol(~M));

LaurenDataFile = 'HD_Data.csv';

FID = fopen(LaurenDataFile, 'r');
% the first line is trash
tline = fgetl(FID);
tline = fgetl(FID);
LaurenHeaders = textscan(tline, '%s', 'Delimiter', ',');
LaurenHeaders = LaurenHeaders{1};

A = textscan(FID, ['%d%s%s%d%d' repmat('%f', 1, length(LaurenHeaders) - 5)], 'Delimiter', ',');

fclose(FID);

for z = 1:length(LaurenHeaders)
	LaurenData.(LaurenHeaders{z}) = A{z};
end

clear A z tline;
%keyboard;

[~, LOC] = ismember(LaurenData.IMAGEHD_CODE, NeuropsychData.ID);

for z = 1:length(LaurenHeaders)
	LaurenData.(LaurenHeaders{z}) = LaurenData.(LaurenHeaders{z})(LOC > 0);
end

NumPerms = 100000;
W = 10;
IDX = W:(NumNodes - W);

% do Lauren's analysis
% Lauren's analysis only applies to TimePoint 1
LaurenData.Group(:) = 1;
%keyboard;
% LOC is the index in the main data of Lauren's subjects
[~, LOC] = ismember(LaurenData.IMAGEHD_CODE, NeuropsychData.ID);
 
CurThicknessProfiles = bsxfun(@rdivide, ThicknessProfiles, reshape(NeuropsychData.Brain_Volumes.total_Intracranial_vol, [1, NumSubjects, NumTimePoints]));
CurThicknessProfiles = CurThicknessProfiles(:, LOC, 1);

NumLaurenSubjects = length(LaurenData.IMAGEHD_CODE);

% BETA PARAMETERS will be [constant, time to mri, age]
BETAParameters = zeros(NumNodes, 3);

for CurNode = 1:size(CurThicknessProfiles, 1)
	
	M = find(~isnan(squeeze(CurThicknessProfiles(CurNode, :))));
	if(length(M) > 2 && ~isempty(M))
		X = [ones(length(M), 1), LaurenData.TimetoMRI_1(M), NeuropsychData.Scan_date_and_age.AGE(LOC(M), 1)];
		XTXINVXT = (X' * X) \ X';
		BETAParameters(CurNode, :) = XTXINVXT * CurThicknessProfiles(CurNode, M)';
		CurThicknessProfiles(CurNode, M) = BETAParameters(CurNode, 1) + LaurenData.TimetoMRI_1(M)' * BETAParameters(CurNode, 2) + NeuropsychData.Scan_date_and_age.AGE(LOC(M), 1)' * BETAParameters(CurNode, 3) - CurThicknessProfiles(CurNode, M);
		clear XTXINVXT X;
		%keyboard;
	end
	clear M;
end

DoLauren = true;

if(DoLauren)
	OutputDir = 'lauren_correlations';
	[~, ~, ~] = mkdir(OutputDir);

	LaurenFields = {'OnsetProbability_5years', ...
		'Education', ...
		'CAG', ...
		'DBS', ...
		'UHDRS_score', ...
		'TrailsB', ...
		'BDI', ...
		'Predicted_FSIQ', ...
		'HVLT_total', ...
		'HVLT_delayed', ...
		'HVLT_rentention', ...
		'HVLT_recognition', ...
		'Tapping', ...
		'ReactionTime', ...
		'FZ_earlyLat', ...
		'FZ_earlyAmp', ...
		'CZ_earlyLat', ...
		'CZ_earlyAmp', ...
		'PZ_earlyLat', ...
		'PZ_earlyAmp', ...
		'FZ_lateLat', ...
		'FZ_lateAmp', ...
		'CZ_lateLat', ...
		'CZ_lateAmp', ...
		'PZ_lateLat', ...
		'PZ_lateAmp', ...
		'FZ_Diff2', ...
		'CZ_Diff2', ...
		'PZ_Diff2', ...
		'FZ_Slope', ...
		'FZ_Intercept', ...
		'CZ_Slope', ...
		'CZ_Intercept', ...
		'PZ_Slope', ...
		'PZ_Intercept', ...
		'Area_ECNV_FZ', ...
		'Area_ECNV_CZ', ...
		'Area_ECNV_PZ', ...
		'Area_LCNV_FZ', ...
		'Area_LCNV_CZ', ...
		'Area_LCNV_PZ'};

	NumCorrelations = length(LaurenFields);

	LaurenGroups = unique(LaurenData.Group);
	LaurenGroups = LaurenGroups(:)';

	NumLaurenGroups = length(LaurenGroups);

	CorrelationTitles = LaurenFields;
	%keyboard;
	CorrelationRValues = cell(NumLaurenGroups, NumCorrelations);
	CorrelationPValues = cell(NumLaurenGroups, NumCorrelations);
	CorrelationPFDRValues = cell(NumLaurenGroups, NumCorrelations);

	%keyboard;
	%for TimePoint = 1:3

	for CurGroupIDX = 1:length(LaurenGroups)
		CurGroup = LaurenGroups(CurGroupIDX);
		CurGroupThicknessProfiles = CurThicknessProfiles(IDX, LaurenData.Group == CurGroup);
		%keyboard;
		for CurCorrelationIDX = 1:NumCorrelations
			CurCorrelationData = LaurenData.(LaurenFields{CurCorrelationIDX})(LaurenData.Group == CurGroup);

			CorrelationRValues{CurGroup, CurCorrelationIDX} = zeros(NumNodes, 1);
			CorrelationPValues{CurGroup, CurCorrelationIDX} = ones(NumNodes, 1);
			CorrelationPFDRValues{CurGroup, CurCorrelationIDX} = ones(NumNodes, 1);
			%keyboard;

			if(any(~isnan(CurCorrelationData(:))))
				%disp([num2str(TimePoint) ' ' num2str(CurGroup) ' ' num2str(CurCorrelationIDX) ' ' num2str(sum(isnan(CurData)))]);
			%else
				for CurNode = 1:length(IDX)
					T = CurGroupThicknessProfiles(CurNode, :)';
					%keyboard;
					M = ~(isnan(CurCorrelationData) | isnan(T));
					if(any(M))
						[R, P] = corrcoef(T(M), CurCorrelationData(M));
						CorrelationRValues{CurGroup, CurCorrelationIDX}(IDX(CurNode)) = R(1, 2);
						CorrelationPValues{CurGroup, CurCorrelationIDX}(IDX(CurNode)) = P(1, 2);
						clear R P;
					end
					clear T M;
				end
				clear CurCorrelationData;
				[~, ~, CorrelationPFDRValues{CurGroup, CurCorrelationIDX}(IDX)] = fdr(CorrelationPValues{CurGroup, CurCorrelationIDX}(IDX));

				% dump results out to a file

				%keyboard;
				DisplayPics = true;
				if(all(~isnan(CorrelationRValues{CurGroup, CurCorrelationIDX})) && DisplayPics)

					cc_seg_tube_p_display(CorrelationPValues{CurGroup, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurCorrelationIDX} == 0), []);
					OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_p.png']);

					FigPos = get(gcf, 'Position');
					set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
					exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

					cc_seg_tube_p_display(CorrelationPFDRValues{CurGroup, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurCorrelationIDX} == 0), []);
					OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_pfdr.png']);

					FigPos = get(gcf, 'Position');
					set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
					exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

					cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX}, [], '\itr');
					OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_r.png']);

					FigPos = get(gcf, 'Position');
					set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
					exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

					cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX} .* (abs(CorrelationPValues{CurGroup, CurCorrelationIDX}) < 0.05), [], '\itr');
					OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpmask.png']);

					FigPos = get(gcf, 'Position');
					set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
					exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
					cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX} .* (abs(CorrelationPFDRValues{CurGroup, CurCorrelationIDX}) < 0.05), [], '\itr');
					OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpfdrmask.png']);

					FigPos = get(gcf, 'Position');
					set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
					exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
				end
			end
		end
		clear CurGroupThicknessProfiles;
	end

	CorrelationTitlesVars = genvarname(CorrelationTitles);
	for CurCorrelationIDX = 1:NumCorrelations
		FID = fopen(fullfile(OutputDir, ['stats_' CorrelationTitlesVars{CurCorrelationIDX} '.csv']), 'w');
		fprintf(FID, 'Node');
		for CurGroupIDX = 1:length(LaurenGroups)
			CurGroup = LaurenGroups(CurGroupIDX);

			T = repmat({GroupNames{CurGroup}}, 1, 3);
			fprintf(FID, ',%s r,%s p,%s p (FDR)', T{:});
			clear T;
		end
		fprintf(FID, '\n');
		T = cat(2, CorrelationRValues(:, CurCorrelationIDX), CorrelationPValues(:, CurCorrelationIDX), CorrelationPFDRValues(:, CurCorrelationIDX))';
		T = cell2mat(T(:)');

		for CurNode = 1:NumNodes
			fprintf(FID, '%d', CurNode);
			fprintf(FID, ',%f', T(CurNode, :));
			fprintf(FID, '\n');
		end
		fclose(FID);
		%keyboard;
	end
end

%keyboard;
Normalisations = cat(3, ones(NumSubjects, NumTimePoints), NeuropsychData.Brain_Volumes.total_Intracranial_vol);
Regressors = cat(3, zeros(NumSubjects, NumTimePoints), NeuropsychData.Scan_date_and_age.AGE);

NormalisationTitles = {'unnorm', 'age_and_icv'};
%keyboard;

DoOrigAnalysis = true;

if(DoOrigAnalysis)
	for z = 2
	%for z = 1:size(Normalisations, 3)

		% original analysis, all vectors in one
		
		S = Normalisations(:, :, z);
		T = Regressors(:, :, z);
		if(all(S(:) == 1) && all(T(:) == 0))
			CurThicknessProfiles = ThicknessProfiles;
			clear S T;
		else
			%keyboard;
			CurThicknessProfiles = bsxfun(@rdivide, ThicknessProfiles, reshape(Normalisations(:, :, z), [1, NumSubjects, NumTimePoints]));
			CurThicknessProfilesTimePoint = bsxfun(@rdivide, ThicknessProfiles, reshape(Normalisations(:, :, z), [1, NumSubjects, NumTimePoints]));
			clear S T;
			%keyboard;
			BETAParametersTimePoint = zeros(NumNodes, 2);
			for CurNode = 1:size(ThicknessProfiles, 1)
				for CurTimePoint = 1:3
					M = ~isnan(CurThicknessProfilesTimePoint(CurNode, :, CurTimePoint));
					if(any(M))
						X = [Regressors(M, CurTimePoint, z), ones(sum(M), 1)];
						%keyboard;
						XTXINVXT = (X' * X) \ X';
						BETAParametersTimePoint(CurNode, :) = XTXINVXT * CurThicknessProfiles(CurNode, M, CurTimePoint)';
						CurThicknessProfilesTimePoint(CurNode, M, CurTimePoint) = CurThicknessProfiles(CurNode, M, CurTimePoint) - ((BETAParametersTimePoint(CurNode, 1) * X(:, 1))' + BETAParametersTimePoint(CurNode, 2));
					end
				end
			end
			
			BETAParameters = zeros(NumNodes, 2);
			for CurNode = 1:size(ThicknessProfiles, 1)
				M = ~isnan(squeeze(CurThicknessProfiles(CurNode, :, :))) & ~isnan(Regressors(:, :, z));
				%keyboard;
				if(sum(M(:)) > 2)
					%keyboard;
					R = Regressors(:, :, z);
					%keyboard;
					X = [R(M), ones(sum(M(:)), 1)];
					clear R;
					XTXINVXT = (X' * X) \ X';
					T = squeeze(CurThicknessProfiles(CurNode, :, :));
					%keyboard;
					BETAParameters(CurNode, :) = XTXINVXT * T(M);%CurThicknessProfiles(CurNode, M, CurTimePoint)';
					%keyboard;
					T = T(M)' - ((BETAParameters(CurNode, 1) * X(:, 1))' + BETAParameters(CurNode, 2));
					S = NaN(NumSubjects, NumTimePoints);
					S(M) = T;
					clear T;
					S = reshape(S, [1, NumSubjects, NumTimePoints]);
					CurThicknessProfiles(CurNode, :, :) = S;
					clear S;

					%keyboard;
				end
			end
		end
		%keyboard;
		
		% do correlations
		if(z == 2)
			OutputDir = [NormalisationTitles{z} '_correlations'];
			[~, ~, ~] = mkdir(OutputDir);
		%	variables of interest

	% 		Disease burden score ???

			% CAG repeats		
			% NeuropsychData.Diagnostic_info.CAG_repeats
			% Stroop
			% NeuropsychData.Neuropsych_data.Stroop
			% SDMT
			% NeuropsychData.Neuropsych_data.SDMT
			% ITIPTAP-slow and –fast, UPSIT
			% NeuropsychData.Motor_function all 

			DataForCorrelation = cat(3, ...
				NeuropsychData.Diagnostic_info.UHDRS, ...
				NeuropsychData.Diagnostic_info.CAG_repeats, ...
				NeuropsychData.Neuropsych_data.Stroop, ...
				NeuropsychData.Neuropsych_data.SDMT, ...
				NeuropsychData.Motor_function.UPSIT_total, ...
				NeuropsychData.Motor_function.ITIPTAP_slow_average, ...
				NeuropsychData.Motor_function.ITIPTAP_fast_average, ...
				NeuropsychData.Motor_function.ITISTAP_Average);
			NumCorrelations = size(DataForCorrelation, 3);

			CorrelationTitles = {'UHDRS', 'CAG repeats', 'Stroop', 'SDMT', 'UPSIT', 'ITIPTAP Slow Average', 'ITIPTAP Fast Average', 'ITIPTAP Average'};
			%keyboard;
			CorrelationRValues = cell(NumGroups, NumCorrelations);
			CorrelationPValues = cell(NumGroups, NumCorrelations);
			CorrelationPFDRValues = cell(NumGroups, NumCorrelations);

			%keyboard;
			%for TimePoint = 1:3
			for CurGroup = 1:NumGroups
				CurGroupThicknessProfiles = CurThicknessProfiles(IDX, NeuropsychData.Group' == CurGroup, :);
				%keyboard;
				for CurCorrelationIDX = 1:NumCorrelations
					CurCorrelationData = DataForCorrelation(NeuropsychData.Group == CurGroup, :, CurCorrelationIDX);

					CorrelationRValues{CurGroup, CurCorrelationIDX} = zeros(NumNodes, 1);
					CorrelationPValues{CurGroup, CurCorrelationIDX} = ones(NumNodes, 1);
					CorrelationPFDRValues{CurGroup, CurCorrelationIDX} = ones(NumNodes, 1);
					%keyboard;

					if(any(~isnan(CurCorrelationData(:))))
						%disp([num2str(TimePoint) ' ' num2str(CurGroup) ' ' num2str(CurCorrelationIDX) ' ' num2str(sum(isnan(CurData)))]);
					%else
						for CurNode = 1:length(IDX)
							T = squeeze(CurGroupThicknessProfiles(CurNode, :, :));
							M = ~(isnan(CurCorrelationData) | isnan(T));
							if(any(M))
								[R, P] = corrcoef(T(M), CurCorrelationData(M));
								CorrelationRValues{CurGroup, CurCorrelationIDX}(IDX(CurNode)) = R(1, 2);
								CorrelationPValues{CurGroup, CurCorrelationIDX}(IDX(CurNode)) = P(1, 2);
								clear R P;
							end
							clear T M;
						end
						clear CurCorrelationData;
						[~, ~, CorrelationPFDRValues{CurGroup, CurCorrelationIDX}(IDX)] = fdr(CorrelationPValues{CurGroup, CurCorrelationIDX}(IDX));

						% dump results out to a file

						%keyboard;
						DisplayPics = false;
						if(all(~isnan(CorrelationRValues{CurGroup, CurCorrelationIDX})) && DisplayPics)

							cc_seg_tube_p_display(CorrelationPValues{CurGroup, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurCorrelationIDX} == 0), []);
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_p.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_p_display(CorrelationPFDRValues{CurGroup, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurCorrelationIDX} == 0), []);
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_pfdr.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX}, [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_r.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX} .* (abs(CorrelationPValues{CurGroup, CurCorrelationIDX}) < 0.05), [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpmask.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurCorrelationIDX} .* (abs(CorrelationPFDRValues{CurGroup, CurCorrelationIDX}) < 0.05), [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpfdrmask.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
						end
					end
				end
				clear CurGroupThicknessProfiles;
			end

			CorrelationTitlesVars = genvarname(CorrelationTitles);
			for CurCorrelationIDX = 1:NumCorrelations
				FID = fopen(fullfile(OutputDir, ['stats_' CorrelationTitlesVars{CurCorrelationIDX} '.csv']), 'w');
				fprintf(FID, 'Node');
				for CurGroup = 1:NumGroups
					T = repmat({GroupNames{CurGroup}}, 1, 3);
					fprintf(FID, ',%s r,%s p,%s p (FDR)', T{:});
					clear T;
				end
				fprintf(FID, '\n');
				T = cat(2, CorrelationRValues(:, CurCorrelationIDX), CorrelationPValues(:, CurCorrelationIDX), CorrelationPFDRValues(:, CurCorrelationIDX))';
				T = cell2mat(T(:)');

				for CurNode = 1:NumNodes
					fprintf(FID, '%d', CurNode);
					fprintf(FID, ',%f', T(CurNode, :));
					fprintf(FID, '\n');
				end
				fclose(FID);
				%keyboard;
			end
			%end
		end
		
		OutputDir = [NormalisationTitles{z} '_correlations_time_specific'];
		[~, ~, ~] = mkdir(OutputDir);
		%	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		CorrelationRValues = cell(NumGroups, NumTimePoints, NumCorrelations);
		CorrelationPValues = cell(NumGroups, NumTimePoints, NumCorrelations);
		CorrelationPFDRValues = cell(NumGroups, NumTimePoints, NumCorrelations);

		%keyboard;
		for CurTimePoint = 1:NumTimePoints
			for CurGroup = 1:NumGroups
				CurGroupThicknessProfiles = CurThicknessProfilesTimePoint(IDX, NeuropsychData.Group' == CurGroup, CurTimePoint);
				%keyboard;
				for CurCorrelationIDX = 1:NumCorrelations
					CurCorrelationData = DataForCorrelation(NeuropsychData.Group == CurGroup, CurTimePoint, CurCorrelationIDX);

					CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX} = zeros(NumNodes, 1);
					CorrelationPValues{CurGroup, CurTimePoint, CurCorrelationIDX} = ones(NumNodes, 1);
					CorrelationPFDRValues{CurGroup, CurTimePoint, CurCorrelationIDX} = ones(NumNodes, 1);
					%keyboard;

					if(any(~isnan(CurCorrelationData(:))))
						%disp([num2str(TimePoint) ' ' num2str(CurGroup) ' ' num2str(CurCorrelationIDX) ' ' num2str(sum(isnan(CurData)))]);
					%else
						for CurNode = 1:length(IDX)
							T = squeeze(CurGroupThicknessProfiles(CurNode, :))';
							%keyboard;
							M = ~(isnan(CurCorrelationData) | isnan(T));
							if(any(M))
								[R, P] = corrcoef(T(M), CurCorrelationData(M));
								CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX}(IDX(CurNode)) = R(1, 2);
								CorrelationPValues{CurGroup, CurTimePoint, CurCorrelationIDX}(IDX(CurNode)) = P(1, 2);
								clear R P;
							end
							clear T M;
						end
						clear CurCorrelationData;
						[~, ~, CorrelationPFDRValues{CurGroup, CurTimePoint, CurCorrelationIDX}(IDX)] = fdr(CorrelationPValues{CurGroup, CurTimePoint, CurCorrelationIDX}(IDX));

						% dump results out to a file

						%keyboard;
						DisplayPics = true;
						if(all(~isnan(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX})) && DisplayPics)

							cc_seg_tube_p_display(CorrelationPValues{CurGroup, CurTimePoint, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX} == 0), []);
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_time_' num2str(CurTimePoint) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_p.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_p_display(CorrelationPFDRValues{CurGroup, CurTimePoint, CurCorrelationIDX}, sign(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX}) + (CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX} == 0), []);
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_time_' num2str(CurTimePoint) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_pfdr.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX}, [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_time_' num2str(CurTimePoint) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_r.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX} .* (abs(CorrelationPValues{CurGroup, CurTimePoint, CurCorrelationIDX}) < 0.05), [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_time_' num2str(CurTimePoint) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpmask.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
							cc_seg_tube_scalar_display(CorrelationRValues{CurGroup, CurTimePoint, CurCorrelationIDX} .* (abs(CorrelationPFDRValues{CurGroup, CurTimePoint, CurCorrelationIDX}) < 0.05), [], '\itr');
							OutputPNG = fullfile(OutputDir, ['group_' num2str(CurGroup) '_time_' num2str(CurTimePoint) '_' strrep(CorrelationTitles{CurCorrelationIDX}, ' ', '') '_rpfdrmask.png']);

							FigPos = get(gcf, 'Position');
							set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
							exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
						end
					end
				end
				clear CurGroupThicknessProfiles;
			end
		end

		CorrelationTitlesVars = genvarname(CorrelationTitles);
		for CurTimePoint = 1:NumTimePoints
			for CurCorrelationIDX = 1:NumCorrelations
				FID = fopen(fullfile(OutputDir, ['stats_time_' num2str(CurTimePoint) '_' CorrelationTitlesVars{CurCorrelationIDX} '.csv']), 'w');
				fprintf(FID, 'Node');
				for CurGroup = 1:NumGroups
					T = repmat({GroupNames{CurGroup}}, 1, 3);
					fprintf(FID, ',%s r,%s p,%s p (FDR)', T{:});
					clear T;
				end
				fprintf(FID, '\n');
				T = cat(2, CorrelationRValues(:, CurTimePoint, CurCorrelationIDX), CorrelationPValues(:, CurTimePoint, CurCorrelationIDX), CorrelationPFDRValues(:, CurTimePoint, CurCorrelationIDX))';
				T = cell2mat(T(:)');

				for CurNode = 1:NumNodes
					fprintf(FID, '%d', CurNode);
					fprintf(FID, ',%f', T(CurNode, :));
					fprintf(FID, '\n');
				end
				fclose(FID);
				%keyboard;
			end
		end
		
		%keyboard;
		%return;
		OutputDir = [NormalisationTitles{z} '_stattest_twosample'];
		[~, ~, ~] = mkdir(OutputDir);

		PermPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedTValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPFDRValues = cell(NumGroups, NumGroups, NumTimePoints);

		clf;
		for TimePoint = 1:NumTimePoints
			for GroupOne = 1:NumGroups - 1
				for GroupTwo = GroupOne + 1:NumGroups
					disp(['Permutation testing, timepoint ' num2str(TimePoint) ', ' num2str(GroupOne) ' vs. ' num2str(GroupTwo)]);
					ThicknessA = CurThicknessProfiles(IDX, NeuropsychData.Group == GroupOne, TimePoint);
					ThicknessB = CurThicknessProfiles(IDX, NeuropsychData.Group == GroupTwo, TimePoint);
					ThicknessA = ThicknessA(:, all(~isnan(ThicknessA), 1));
					ThicknessB = ThicknessB(:, all(~isnan(ThicknessB), 1));
					[PermPValues{GroupOne, GroupTwo, TimePoint}, ...
						ObservedTValues{GroupOne, GroupTwo, TimePoint}, ...
						~, ...
						ObservedPValues{GroupOne, GroupTwo, TimePoint}] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms);
					[~, ~, ObservedPFDRValues{GroupOne, GroupTwo, TimePoint}] = fdr(ObservedPValues{GroupOne, GroupTwo, TimePoint}, 0.05);

					clear ThicknessA ThicknessB;
				end
			end
		end
		FID = fopen(fullfile(OutputDir, 'stats_group.csv'), 'w');
		fprintf(FID, 'Groups,Timepoint,Node,T,P,P (FDR),P (Permutation)\n');
		%%
		for TimePoint = 1:NumTimePoints
			for GroupOne = 1:NumGroups - 1
				for GroupTwo = GroupOne + 1:NumGroups
					GroupLabels = {GroupNames{GroupOne}, GroupNames{GroupTwo}};%['G' num2str(GroupOne)], ['G' num2str(GroupTwo)]};

					for CurNode = 1:length(IDX)
						fprintf(FID, '%svs%s,%d,%d', GroupLabels{1}, GroupLabels{2}, TimePoint, IDX(CurNode));
						fprintf(FID, ',%f,%f,%f,%f\n', ...
							ObservedTValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPFDRValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							PermPValues{GroupOne, GroupTwo, TimePoint}(CurNode));
					end
					OutputPNG = fullfile(OutputDir, ['perm_timepoint_' num2str(TimePoint) '_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						PermPValues{GroupOne, GroupTwo, TimePoint}, ...
						ObservedTValues{GroupOne, GroupTwo, TimePoint}, ...
						IDX, GroupLabels, OutputPNG);

					OutputPNG = fullfile(OutputDir, ['fdr_timepoint_' num2str(TimePoint) '_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPFDRValues{GroupOne, GroupTwo, TimePoint}, ...
						ObservedTValues{GroupOne, GroupTwo, TimePoint}, ...
						IDX, GroupLabels, OutputPNG);
					%keyboard;
					OutputPNG = fullfile(OutputDir, ['observed_timepoint_' num2str(TimePoint) '_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPValues{GroupOne, GroupTwo, TimePoint}, ...
						ObservedTValues{GroupOne, GroupTwo, TimePoint}, ...
						IDX, GroupLabels, OutputPNG);
				end
			end
		end
		fclose(FID);
		%keyboard;
		%
		PermPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedTValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPFDRValues = cell(NumGroups, NumGroups, NumTimePoints);


		clf;
		for Group = 1:NumGroups
			for TimePointOne = 1:NumTimePoints - 1
				for TimePointTwo = TimePointOne + 1:NumTimePoints
					disp(['Permutation testing, group ' num2str(Group) ', ' num2str(TimePointOne) ' vs. ' num2str(TimePointTwo)]);
					ThicknessA = CurThicknessProfiles(IDX, NeuropsychData.Group == Group, TimePointOne);
					ThicknessB = CurThicknessProfiles(IDX, NeuropsychData.Group == Group, TimePointTwo);
					ThicknessA = ThicknessA(:, all(~isnan(ThicknessA), 1));
					ThicknessB = ThicknessB(:, all(~isnan(ThicknessB), 1));
					[PermPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						~, ...
						ObservedPValues{TimePointOne, TimePointTwo, Group}] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms);
					[~, ~, ObservedPFDRValues{TimePointOne, TimePointTwo, Group}] = fdr(ObservedPValues{TimePointOne, TimePointTwo, Group}, 0.05);
					clear ThicknessA ThicknessB;
				end
			end
		end

		FID = fopen(fullfile(OutputDir, 'stats_time.csv'), 'w');
		fprintf(FID, 'Group,Timepoints,Node,T,P,P (FDR),P (Permutation)\n');

		%%
		for Group = 1:NumGroups
			for TimePointOne = 1:NumTimePoints - 1
				for TimePointTwo = TimePointOne + 1:NumTimePoints

					GroupLabels = {['t_' num2str(TimePointOne)], ['t_' num2str(TimePointTwo)]};

					for CurNode = 1:length(IDX)
						fprintf(FID, '%s,%svs%s,%d', GroupNames{Group}, GroupLabels{1}, GroupLabels{2}, IDX(CurNode));
						fprintf(FID, ',%f,%f,%f,%f\n', ...
							ObservedTValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPFDRValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							PermPValues{GroupOne, GroupTwo, TimePoint}(CurNode));
					end
					OutputPNG = fullfile(OutputDir, ['perm_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						PermPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
					OutputPNG = fullfile(OutputDir, ['fdr_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPFDRValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
					OutputPNG = fullfile(OutputDir, ['observed_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
				end
			end
		end
		fclose(FID);
		% only the subjects that had all three timepoints, paired testing
		PermPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedTValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPValues = cell(NumGroups, NumGroups, NumTimePoints);
		ObservedPFDRValues = cell(NumGroups, NumGroups, NumTimePoints);

		AllThreeTimePoints = all(~isnan(CurThicknessProfiles(1, :, :)), 3);

		clf;
		for Group = 1:NumGroups
			for TimePointOne = 1:NumTimePoints - 1
				for TimePointTwo = TimePointOne + 1:NumTimePoints
					disp(['Permutation testing, group ' num2str(Group) ', ' num2str(TimePointOne) ' vs. ' num2str(TimePointTwo)]);
					ThicknessA = CurThicknessProfiles(IDX, NeuropsychData.Group == Group & AllThreeTimePoints', TimePointOne);
					ThicknessB = CurThicknessProfiles(IDX, NeuropsychData.Group == Group & AllThreeTimePoints', TimePointTwo);
					ThicknessA = ThicknessA(:, all(~isnan(ThicknessA), 1));
					ThicknessB = ThicknessB(:, all(~isnan(ThicknessB), 1));
					if(~isequal(size(ThicknessA), size(ThicknessB)))
						keyboard;
					end
					[PermPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						~, ...
						ObservedPValues{TimePointOne, TimePointTwo, Group}] = cc_seg_pairedT_perm(ThicknessA, ThicknessB, NumPerms);
					[~, ~, ObservedPFDRValues{TimePointOne, TimePointTwo, Group}] = fdr(ObservedPValues{TimePointOne, TimePointTwo, Group}, 0.05);

					clear ThicknessA ThicknessB;

				end
			end
		end

		OutputDir = [NormalisationTitles{z} '_stattest_paired'];
		[~, ~, ~] = mkdir(OutputDir);

		FID = fopen(fullfile(OutputDir, 'stats_time.csv'), 'w');
		fprintf(FID, 'Group,Timepoints,Node,T,P,P (FDR),P (Permutation)\n');

		%%
		for Group = 1:NumGroups
			for TimePointOne = 1:NumTimePoints
				for TimePointTwo = TimePointOne + 1:NumTimePoints
					GroupLabels = {['t_' num2str(TimePointOne)], ['t_' num2str(TimePointTwo)]};

					for CurNode = 1:length(IDX)
						fprintf(FID, '%s,%svs%s,%d', GroupNames{Group}, GroupLabels{1}, GroupLabels{2}, IDX(CurNode));
						fprintf(FID, ',%f,%f,%f,%f\n', ...
							ObservedTValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							ObservedPFDRValues{GroupOne, GroupTwo, TimePoint}(CurNode), ...
							PermPValues{GroupOne, GroupTwo, TimePoint}(CurNode));
					end

					OutputPNG = fullfile(OutputDir, ['perm_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						PermPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
					OutputPNG = fullfile(OutputDir, ['fdr_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPFDRValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
					OutputPNG = fullfile(OutputDir, ['observed_group_' num2str(Group) '_TimePoint_' num2str(TimePointOne) 'vs' num2str(TimePointTwo) '.png']);
					cc_seg_image_hd_make_image(NumNodes, ...
						ObservedPValues{TimePointOne, TimePointTwo, Group}, ...
						ObservedTValues{TimePointOne, TimePointTwo, Group}, ...
						IDX, GroupLabels, OutputPNG);
				end
			end
		end
		fclose(FID);
		%keyboard;
		% perform regressions on each subject over time points
		% PermPValues = cell(3, 3);
		% ObservedTValues = cell(3, 3);
		% ObservedPValues = cell(3, 3);
		% ObservedPFDRValues = cell(3, 3);

		X = [(1:3)', ones(3, 1)];
		XTXINVXT = (X' * X) \ X';

		T = permute(CurThicknessProfiles, [3, 1, 2]);
		T = reshape(T, [3, size(T, 2) * size(T, 3)]);
		B = XTXINVXT * T;

		BETAProfiles = reshape(B(1, :), size(CurThicknessProfiles, 1), size(CurThicknessProfiles, 2));

		% for I = 1:size(CurThicknessProfiles, 1)
		% 	for I = 1:size(CurThicknessProfiles, 1)
		% 	end
		% end
		%%
		OutputDir = [NormalisationTitles{z} '_stattest_regression'];
		[~, ~, ~] = mkdir(OutputDir);
		NumNodes = size(CurThicknessProfiles, 1);
		GroupBETAMeans = zeros(NumNodes, 3);
		for Group = 1:NumGroups
			GroupBETAMeans(:, Group) = mean(BETAProfiles(:, NeuropsychData.Group == Group & AllThreeTimePoints(:)), 2);
			OutputPNG = fullfile(OutputDir, ['beta_group' num2str(Group)]);
			cc_seg_tube_scalar_display(GroupBETAMeans(:, Group), [], '\beta');
			FigPos = get(gcf, 'Position');
			set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'Color', 'k', 'InvertHardCopy', 'off');
			exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
			delete(gcf);
		end

		PermPValues = cell(3, 3);
		ObservedTValues = cell(3, 3);
		ObservedPValues = cell(3, 3);
		ObservedPFDRValues = cell(3, 3);

		clf;
		for GroupOne = 1:NumGroups - 1
			for GroupTwo = GroupOne + 1:NumGroups
				disp(['Permutation testing for beta, group ' num2str(GroupOne) ' vs. ' num2str(GroupTwo)]);
				ThicknessA = abs(BETAProfiles(IDX, NeuropsychData.Group == GroupOne & AllThreeTimePoints'));
				ThicknessB = abs(BETAProfiles(IDX, NeuropsychData.Group == GroupTwo & AllThreeTimePoints'));
				ThicknessA = ThicknessA(:, all(~isnan(ThicknessA), 1));
				ThicknessB = ThicknessB(:, all(~isnan(ThicknessB), 1));
				[PermPValues{GroupOne, GroupTwo}, ...
					ObservedTValues{GroupOne, GroupTwo}, ...
					~, ...
					ObservedPValues{GroupOne, GroupTwo}] = cc_seg_twosampleT_perm(ThicknessA, ThicknessB, NumPerms);
				[~, ~, ObservedPFDRValues{GroupOne, GroupTwo}] = fdr(ObservedPValues{GroupOne, GroupTwo}, 0.05);

		% 			subplot(SR, SC, PlotIDX(GroupOne, GroupTwo));
		% 			plot_pos_neg_p_values(PermPValues{TimePointOne, TimePointTwo, Group}, ObservedTValues{GroupOne, GroupTwo, TimePoint}, IDX, {['Group ' num2str(GroupOne)], ['Group ' num2str(GroupTwo)]});
		% 			title([num2str(size(ThicknessA, 2)) ' vs. ' num2str(size(ThicknessB, 2))]);
				clear ThicknessA ThicknessB;

			end
		end

		for GroupOne = 1:NumGroups - 1
			for GroupTwo = GroupOne + 1:NumGroups
				GroupLabels = {GroupNames{GroupOne}, GroupNames{GroupTwo}};%['G' num2str(GroupOne)], ['G' num2str(GroupTwo)]};

				OutputPNG = fullfile(OutputDir, ['perm_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
				cc_seg_image_hd_make_image(NumNodes, ...
					PermPValues{GroupOne, GroupTwo}, ...
					ObservedTValues{GroupOne, GroupTwo}, ...
					IDX, GroupLabels, OutputPNG);

				OutputPNG = fullfile(OutputDir, ['fdr_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
				cc_seg_image_hd_make_image(NumNodes, ...
					ObservedPFDRValues{GroupOne, GroupTwo}, ...
					ObservedTValues{GroupOne, GroupTwo}, ...
					IDX, GroupLabels, OutputPNG);
				OutputPNG = fullfile(OutputDir, ['observed_group_' num2str(GroupOne) 'vs' num2str(GroupTwo) '.png']);
				cc_seg_image_hd_make_image(NumNodes, ...
					ObservedPValues{GroupOne, GroupTwo}, ...
					ObservedTValues{GroupOne, GroupTwo}, ...
					IDX, GroupLabels, OutputPNG);
			end
		end

		%keyboard;
		%%
		%cc_seg_tube_p_display(P, PSign, GroupLabels)
		delete(gcf);

		%%
	end
end

return;
% % 
% % MeanAllThree = squeeze(mean(CurThicknessProfiles(:, AllThreeTimePoints, :), 2));
% % 
% % surf(MeanAllThree');
% % axis equal;
% % axis ij;
% % view(-71.5, 12);
% % xlabel('Node');
% % ylabel('Timepoint');
% % FigPos = fullscreen_fig_pos;
% % OutputPNG = 'mean_thicknesses_1_to_3.png';
% % set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% % exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% % 
% % clf;
% % plot(MeanAllThree);
% % xlabel('Node');
% % ylabel('Thickness');
% % FigPos = fullscreen_fig_pos;
% % OutputPNG = 'mean_thicknesses_1_to_3_plot.png';
% % set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% % exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% % keyboard;
% NumNodes = size(CurThicknessProfiles, 1);
% 
% 
% OneMinusTwoPValues = zeros(1, NumNodes);
% OneMinusTwoTValues = zeros(1, NumNodes);
% OneMinusThreePValues = zeros(1, NumNodes);
% OneMinusThreeTValues = zeros(1, NumNodes);
% TwoMinusThreePValues = zeros(1, NumNodes);
% TwoMinusThreeTValues = zeros(1, NumNodes);
% 
% for z = 1:NumNodes
% 	[~,P,~,STATS] = ttest(CurThicknessProfiles(z, AllThreeTimePoints, 1), CurThicknessProfiles(z, AllThreeTimePoints, 2), 0.05, 'right');
% 	OneMinusTwoPValues(z) = P;
% 	OneMinusTwoTValues(z) = STATS.tstat;
% 	
% 	[~,P,~,STATS] = ttest(CurThicknessProfiles(z, AllThreeTimePoints, 1), CurThicknessProfiles(z, AllThreeTimePoints, 3), 0.05, 'right');
% 	OneMinusThreePValues(z) = P;
% 	OneMinusThreeTValues(z) = STATS.tstat;
% 	[~,P,~,STATS] = ttest(CurThicknessProfiles(z, AllThreeTimePoints, 2), CurThicknessProfiles(z, AllThreeTimePoints, 3), 0.05, 'right');
% 	TwoMinusThreePValues(z) = P;
% 	TwoMinusThreeTValues(z) = STATS.tstat;
% end
% 
% IDX = 10:90;
% 
% OneMinusTwoPValues = OneMinusTwoPValues(IDX);
% OneMinusTwoTValues = OneMinusTwoTValues(IDX);
% OneMinusThreePValues = OneMinusThreePValues(IDX);
% OneMinusThreeTValues = OneMinusThreeTValues(IDX);
% TwoMinusThreePValues = TwoMinusThreePValues(IDX);
% TwoMinusThreeTValues = TwoMinusThreeTValues(IDX);
% 
% OneMinusTwoPValuesFDRThresh = fdr(OneMinusTwoPValues, 0.05);
% OneMinusThreePValuesFDRThresh = fdr(OneMinusThreePValues, 0.05);
% TwoMinusThreePValuesFDRThresh = fdr(TwoMinusThreePValues, 0.05);
% 
% if(~isempty(OneMinusTwoPValuesFDRThresh))
% 	OneMinusTwoPValuesFDR = min(1, (0.05 ./ OneMinusTwoPValuesFDRThresh) .* OneMinusTwoPValues);
% else
% 	OneMinusTwoPValuesFDR = ones(size(OneMinusTwoPValues));
% end
% 
% if(~isempty(OneMinusThreePValuesFDRThresh))
% 	OneMinusThreePValuesFDR = min(1, (0.05 ./ OneMinusThreePValuesFDRThresh) .* OneMinusThreePValues);
% else
% 	OneMinusThreePValuesFDR = ones(size(OneMinusThreePValues));
% end
% 
% if(~isempty(TwoMinusThreePValuesFDRThresh))
% 	TwoMinusThreePValuesFDR = min(1, (0.05 ./ TwoMinusThreePValuesFDRThresh) .* TwoMinusThreePValues);
% else
% 	TwoMinusThreePValuesFDR = ones(size(TwoMinusThreePValues));
% end
% 
% disp('Permutation testing, 1 - 2');
% [OneMinusTwoPermPValues, OneMinusTwoPermSign, ~, OneMinusTwoPermPValuesP] = cc_seg_pairedT_perm(CurThicknessProfiles(IDX, AllThreeTimePoints, 1), CurThicknessProfiles(IDX, AllThreeTimePoints, 2));
% disp('Permutation testing, 1 - 3');
% [OneMinusThreePermPValues, OneMinusThreePermSign, ~, OneMinusThreePermPValuesP] = cc_seg_pairedT_perm(CurThicknessProfiles(IDX, AllThreeTimePoints, 1), CurThicknessProfiles(IDX, AllThreeTimePoints, 3));
% disp('Permutation testing, 2 - 3');
% [TwoMinusThreePermPValues, TwoMinusThreePermSign, ~, TwoMinusThreePermPValuesP] = cc_seg_pairedT_perm(CurThicknessProfiles(IDX, AllThreeTimePoints, 2), CurThicknessProfiles(IDX, AllThreeTimePoints, 3));
% %%
% clf;
% AX = zeros(1, 3);
% AX(1) = subplot(2, 2, 1);
% plot_pos_neg_p_values(OneMinusTwoPValues, OneMinusTwoTValues, IDX);
% %plot(1:NumNodes, OneMinusTwoPValues .* sign(OneMinusTwoTValues), '*');
% title('Time 1 - Time 2');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 2', 'Time 1 < Time 2'});
% AX(2) = subplot(2, 2, 2);
% plot_pos_neg_p_values(OneMinusThreePValues, OneMinusThreeTValues, IDX);
% %plot(1:NumNodes, OneMinusThreePValues .* sign(OneMinusThreeTValues), '*');
% title('Time 1 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 3', 'Time 1 < Time 3'});
% AX(3) = subplot(2, 2, 4);
% plot_pos_neg_p_values(TwoMinusThreePValues, TwoMinusThreeTValues, IDX);
% %plot(1:NumNodes, TwoMinusThreePValues .* sign(TwoMinusThreeTValues), '*');
% title('Time 2 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 2 > Time 3', 'Time 2 < Time 3'});
% 
% set(AX, 'YTick', [0.05, 0.5, 1]);
% OutputPNG = 'paired_t_1_to_3.png';
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 
% clf;
% AX = zeros(1, 3);
% AX(1) = subplot(2, 2, 1);
% plot_pos_neg_p_values(OneMinusTwoPValuesFDR, OneMinusTwoTValues, IDX);
% %plot(1:NumNodes, OneMinusTwoPValues .* sign(OneMinusTwoTValues), '*');
% title('Time 1 - Time 2');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 2', 'Time 1 < Time 2'});
% AX(2) = subplot(2, 2, 2);
% plot_pos_neg_p_values(OneMinusThreePValuesFDR, OneMinusThreeTValues, IDX);
% %plot(1:NumNodes, OneMinusThreePValues .* sign(OneMinusThreeTValues), '*');
% title('Time 1 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 3', 'Time 1 < Time 3'});
% AX(3) = subplot(2, 2, 4);
% plot_pos_neg_p_values(TwoMinusThreePValuesFDR, TwoMinusThreeTValues, IDX);
% %plot(1:NumNodes, TwoMinusThreePValues .* sign(TwoMinusThreeTValues), '*');
% title('Time 2 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 2 > Time 3', 'Time 2 < Time 3'});
% 
% set(AX, 'YTick', [0.05, 0.5, 1]);
% OutputPNG = 'paired_t_1_to_3_fdr.png';
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% clf;
% AX = zeros(1, 3);
% AX(1) = subplot(2, 2, 1);
% plot_pos_neg_p_values(OneMinusTwoPermPValuesP, OneMinusTwoPermSign, IDX);
% %plot(1:NumNodes, OneMinusTwoPValues .* sign(OneMinusTwoTValues), '*');
% title('Time 1 - Time 2');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 2', 'Time 1 < Time 2'});
% AX(2) = subplot(2, 2, 2);
% plot_pos_neg_p_values(OneMinusThreePermPValuesP, OneMinusThreePermSign, IDX);
% %plot(1:NumNodes, OneMinusThreePValues .* sign(OneMinusThreeTValues), '*');
% title('Time 1 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 3', 'Time 1 < Time 3'});
% AX(3) = subplot(2, 2, 4);
% plot_pos_neg_p_values(TwoMinusThreePermPValuesP, TwoMinusThreePermSign, IDX);
% %plot(1:NumNodes, TwoMinusThreePValues .* sign(TwoMinusThreeTValues), '*');
% title('Time 2 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 2 > Time 3', 'Time 2 < Time 3'});
% 
% set(AX, 'YTick', [0.05, 0.5, 1]);
% OutputPNG = 'paired_t_1_to_3_observed.png';
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% %keyboard;
% %%%%% permutation test
% 
% clf;
% AX = zeros(1, 3);
% AX(1) = subplot(2, 2, 1);
% plot_pos_neg_p_values(OneMinusTwoPermPValues, OneMinusTwoPermSign, IDX);
% 
% 
% %plot(1:NumNodes, OneMinusTwoPValues .* sign(OneMinusTwoTValues), '*');
% title('Time 1 - Time 2');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 2', 'Time 1 < Time 2'});
% AX(2) = subplot(2, 2, 2);
% plot_pos_neg_p_values(OneMinusThreePermPValues, OneMinusThreePermSign, IDX);
% %plot(1:NumNodes, OneMinusThreePValues .* sign(OneMinusThreeTValues), '*');
% title('Time 1 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 1 > Time 3', 'Time 1 < Time 3'});
% AX(3) = subplot(2, 2, 4);
% plot_pos_neg_p_values(TwoMinusThreePermPValues, TwoMinusThreePermSign, IDX);
% %plot(1:NumNodes, TwoMinusThreePValues .* sign(TwoMinusThreeTValues), '*');
% title('Time 2 - Time 3');
% line([0, NumNodes], [0.05, 0.05], 'Color', 'k');
% xlabel('Node');
% ylabel('\itp');
% legend({'Time 2 > Time 3', 'Time 2 < Time 3'});
% 
% set(AX, 'YTick', [0.05, 0.5, 1]);
% OutputPNG = 'paired_t_1_to_3_perm.png';
% set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% delete(gcf);
% %%
% % [~, I] = CurThicknessProfiles(:, AllThreeTimePoints, :),min([Output.FinalDice]);
% % hold off;
% % imshow(Output(I).IMG, []);
% % hold on;
% % [~, CC] = contour(Output(I).GroundSeg, [0.5, 0.5]);
% % set(CC, 'Color', 'b');
% % [~, CC] = contour(Output(I).FinalSeg, [0.5, 0.5]);
% % set(CC, 'Color', 'r');
% % axis equal ij;
% %delete(gcf);
% % 
% % FID = fopen('convex_hull_used.txt', 'w');
% % fprintf(FID, 'Used convex hull:\n');
% % fprintf(FID, '%d\n', find(ConvexHullUsed));
% % fprintf(FID, 'Didnt use convex hull:\n');
% % fprintf(FID, '%d\n', find(~ConvexHullUsed));
% % fclose(FID);
% 
% % for IDX = 1:length(ConvexHullSubjects)
% % 	z = ConvexHullSubjects(IDX);
% % 	disp(['Subject ' num2str(z)]);
% % 	BaseName = strrep(D(z).name, '.nii.gz', '');
% % 	CCName = strrep(D(z).name, 'msp', 'cc');
% % 	%cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), fullfile('oasis_database', 'flipped', CCName), fullfile(FigureDir, [sprintf('%04d', z) '-' BaseName]));
% % 	ConvexHullUsed(z) = cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), ...
% % 	fullfile('oasis_database', 'flipped', D(z).name), ...
% % 	OutputDir, ...
% % 	[sprintf('%04d', z) '-' BaseName]);
% % end
% % delete(gcf);
% 
% %return;
% % MeanOriginal = mean([Overlap.Original]);
% % MeanOriginalW = mean([Overlap.OriginalW]);
% % MeanSeg = nanmean([Overlap.Seg]);
% % MeanSegW = nanmean([Overlap.SegW]);
% % 
% % MinOriginal = min([Overlap.Original]);
% % MinOriginalW = min([Overlap.OriginalW]);
% % MinSeg = nanmin([Overlap.Seg]);
% % MinSegW = nanmin([Overlap.SegW]);
% % 
% % MaxOriginal = max([Overlap.Original]);
% % MaxOriginalW = max([Overlap.OriginalW]);
% % MaxSeg = nanmax([Overlap.Seg]);
% % MaxSegW = nanmax([Overlap.SegW]);
% % 
% % disp(['Original: ' num2str(MeanOriginal) ' [' num2str(MinOriginal) ', ' num2str(MaxOriginal) ']']);
% % disp(['OriginalW: ' num2str(MeanOriginalW) ' [' num2str(MinOriginalW) ', ' num2str(MaxOriginalW) ']']);
% % disp(['Seg: ' num2str(MeanSeg) ' [' num2str(MinSeg) ', ' num2str(MaxSeg) ']']);
% % disp(['SegW: ' num2str(MeanSegW) ' [' num2str(MinSegW) ', ' num2str(MaxSegW) ']']);
% % 
% % save lk_overlap Overlap Mean* Min* Max*;