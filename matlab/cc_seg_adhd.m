clear;

InputDir = 'RawData';
D = dir(fullfile(InputDir, '*.nii.gz'));
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;

ExistingOutputFiles = findallfiles(OutputDir);
ExistingOutputFilePaths = cell(1, length(ExistingOutputFiles));
ExistingOutputFileNames = cell(1, length(ExistingOutputFiles));
ExistingOutputFileExtensions = cell(1, length(ExistingOutputFiles));

for z = 1:length(ExistingOutputFiles)
	[ExistingOutputFilePaths{z}, ExistingOutputFileNames{z}, ExistingOutputFileExtensions{z}] = fileparts(ExistingOutputFiles{z});
end

tokens = regexp(ExistingOutputFileNames, '^(\d\d\d\d)-(.*)$', 'tokens');

ExistingOutputFileNumbers = zeros(1, length(ExistingOutputFiles));
ExistingOutputFileSuffixes = cell(1, length(ExistingOutputFiles));

for z = 1:length(ExistingOutputFiles)
	ExistingOutputFileNumbers(z) = str2double(tokens{z}{1}{1});
	ExistingOutputFileSuffixes{z} = tokens{z}{1}{2};
end

Mode = 'skip';

IDX = 1:length(D);

[OutputBases{1:length(D)}] = deal(D(IDX).name);
OutputBases = strrep(OutputBases, '.nii.gz', '');

% for each existing output file, check to see where it is in the current
% list
ExistingOutputFileNewNumbers = zeros(1, length(ExistingOutputFiles));

for z = 1:length(OutputBases)
	T = strncmp(OutputBases{z}, ExistingOutputFileSuffixes, length(OutputBases{z}));
	F = find(T);
	if(~isempty(F))
		ExistingOutputFileNewNumbers(F) = z;
	end
	clear T F;
end

ExistingOutputFileNumbersNotEqual = find(ExistingOutputFileNumbers ~= ExistingOutputFileNewNumbers);

if(~isempty(ExistingOutputFileNumbersNotEqual))
	for z = 1:length(ExistingOutputFileNumbersNotEqual)
		disp([ExistingOutputFiles{ExistingOutputFileNumbersNotEqual(z)} ' -> ' ...
			fullfile(ExistingOutputFilePaths{ExistingOutputFileNumbersNotEqual(z)}, ...
			[sprintf('%04d', ExistingOutputFileNewNumbers(ExistingOutputFileNumbersNotEqual(z))) '-' ExistingOutputFileSuffixes{ExistingOutputFileNumbersNotEqual(z)} ExistingOutputFileExtensions{ExistingOutputFileNumbersNotEqual(z)}])]);
		movefile(ExistingOutputFiles{ExistingOutputFileNumbersNotEqual(z)}, ...
			fullfile(ExistingOutputFilePaths{ExistingOutputFileNumbersNotEqual(z)}, ...
			[sprintf('%04d', ExistingOutputFileNewNumbers(ExistingOutputFileNumbersNotEqual(z))) '-' ExistingOutputFileSuffixes{ExistingOutputFileNumbersNotEqual(z)} ExistingOutputFileExtensions{ExistingOutputFileNumbersNotEqual(z)}]));
		
	end
end
%keyboard;
%IDX = 29:length(D);
%IDX = 9;
for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	switch(lower(Mode))
		case 'seg_and_thickness'
		% 	if(exist([OutputBase '_seg.mat'], 'file') ~= 2)
% 				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
% 				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
% 				[], ...
% 				OutputDir, ...
% 				[sprintf('%04d', z) '-' BaseName], 'flirt');
% 				if(SegReturnCode == 0)
% 					ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
% 				end
		%	end
			if(exist([OutputBase '_seg.mat'], 'file') ~= 2 || exist([OutputBase '_thickness.mat'], 'file') ~= 2 )
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName], 'flirt');
				%ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			end
		case 'manedit_thickness'
			%if(exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			%end
	end
	%Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
end

if(~strcmp(Mode, 'skip'))
	delete(gcf);
end

% load all thickness profiles
ThicknessProfiles = cell(length(D), 1);

for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	T = load([OutputBase '_thickness'], 'Thickness');
	ThicknessProfiles{z} = T.Thickness;
	clear T;
end


FID = fopen('ADHD_CC_pairings.csv', 'r');
tline = fgetl(FID);
CSVHeaders = textscan(tline, '%q', 'Delimiter', ',');
CSVHeaders = CSVHeaders{1};

CSVHeaders = strrep(CSVHeaders, ' ', '');
%"24/10/2006","5/11/1996",10,"ADHD_MPH_101",100,104,102,104,"R"
A = textscan(FID, '%q%q%f%q%f%f%f%f%q', 'Delimiter', ',');

for z = 1:length(A)
	CSVData.(CSVHeaders{z}) = A{z};
end
clear A;

fclose(FID);

% find indices for each of the subjects in the CSV file
% select the first one if there are multiples

CSVSegID = zeros(1, length(CSVData.code));

for z = 1:length(CSVData.code)
	I = find(strncmp(CSVData.code{z}, OutputBases, length(CSVData.code{z})), 1, 'first');
	if(~isempty(I))
		CSVSegID(z) = I;
	end
	%keyboard;
end

F = fieldnames(CSVData);

for z = 1:length(F)
	CSVData.(F{z}) = CSVData.(F{z})(CSVSegID > 0);
end
CSVSegID = CSVSegID(CSVSegID > 0);

Age = NaN(1, length(D));
Age(CSVSegID) = CSVData.age;

%keyboard;
Exclude = [44 114];

CTLIDX = find(strncmp('CTL', OutputBases, 3));
ADHDIDX = find(strncmp('ADHD', OutputBases, 4));

ADHDIDX = [ADHDIDX 118:125];

BRIIDX = cellfun(@(x) (strfind(x, '_BRI_')), OutputBases, 'UniformOutput', false);
BRIIDX = find(~cellfun(@isempty, BRIIDX));

%ADHDIDX = [1:24, 58:65];
%CTLIDX = 25:57;

ADHDIDX = setdiff(ADHDIDX, Exclude);
CTLIDX = setdiff(CTLIDX, Exclude);

ADHDIDX = setdiff(ADHDIDX, BRIIDX);
CTLIDX = setdiff(CTLIDX, BRIIDX);

% select the IDs from the CTL and ADHD that are in the CSV file
ADHDIDX = intersect(ADHDIDX, CSVSegID);
CTLIDX = intersect(CTLIDX, CSVSegID);

%keyboard;
ThicknessCTL = cat(2, ThicknessProfiles{CTLIDX});
ThicknessADHD = cat(2, ThicknessProfiles{ADHDIDX});

%%
XC = bsxfun(@minus, ThicknessCTL, mean(ThicknessCTL));

%XC = bsxfun(@rdivide, ThicknessCTL, sum(ThicknessCTL));
PolyFitRaw = zeros(size(ThicknessCTL, 1), 2);
PolyFitXC = zeros(size(ThicknessCTL, 1), 2);

RegressP = zeros(size(ThicknessCTL, 1), 1);
RegressB = zeros(size(ThicknessCTL, 1), 2);
for z = 1:size(ThicknessCTL, 1)
	PolyFitRaw(z, :) = polyfit(Age(CTLIDX), ThicknessCTL(z, :), 1);
	[B,~,~,~,STATS] = regress(ThicknessCTL(z, :)', [Age(CTLIDX)', ones(length(CTLIDX), 1)]);
	RegressB(z, :) = B';
	RegressP(z) = STATS(3);
	
	PolyFitXC(z, :) = polyfit(Age(CTLIDX), XC(z, :), 1);
end
% 
% subplot 121;
% hold off;
% plot(PolyFitRaw(:, 1));
% 
% subplot 122;
% hold off;
% plot(PolyFitXC(:, 1));
% hold off;
% plot(Age(CTLIDX), XC(C, :), '*');
% hold on;
% B = polyfit(Age(CTLIDX), XC(C, :), 1);

cc_seg_tube_scalar_display(PolyFitRaw(:, 1), '\beta');
FigPos = fullscreen_fig_pos;
set(gcf, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off');
OutputFile = 'thickness_slope_age';
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

IMG = imread([OutputFile, '.png']);
IMG = imcomplement(IMG);
IMG = imautocropwhite(IMG, 20);
IMG = imcomplement(IMG);
imwrite(IMG, [OutputFile, '.png']);
delete(gcf);
% return;
% %%
% C = 13;
% 
% 
% 
% %B = regress(Age(CTLIDX), cat(2, ThicknessCTL(C, :)', ones(size(ThicknessCTL, 2), 1)));
% %BXC = regress(Age(CTLIDX), cat(2, XC(C, :)', ones(size(ThicknessCTL, 2), 1)));
% 
% CTLAgeRange = [min(Age(CTLIDX)) max(Age(CTLIDX))];
% subplot 121;
% hold off;
% plot(Age(CTLIDX), ThicknessCTL(C, :), '*');
% hold on;
% B = polyfit(Age(CTLIDX), ThicknessCTL(C, :), 1);
% plot(CTLAgeRange, polyval(B, CTLAgeRange));
% 
% subplot 122;
% hold off;
% plot(Age(CTLIDX), XC(C, :), '*');
% hold on;
% B = polyfit(Age(CTLIDX), XC(C, :), 1);
% plot(CTLAgeRange, polyval(B, CTLAgeRange));
% 
% % set(H1, 'LineStyle', 'none', 'Marker', '*');
% % set(H2, 'LineStyle', 'none', 'Marker', '*');
%%
NumPerms = 100000;

W = 10;
MIDX = W:size(ThicknessCTL, 1) - W;
disp(['Permutation testing: ' num2str(size(ThicknessCTL, 2)) ' vs. ' num2str(size(ThicknessADHD, 2))]);
[PermPValues, PermTSignValues, ~, ObservedPValues] = cc_seg_twosampleT_perm(ThicknessCTL(MIDX, :), ThicknessADHD(MIDX, :), NumPerms);
[~, ~, ObservedPFDRValues] = fdr(ObservedPValues, 0.05);
%%
subplot 221;
plot([nanmean(ThicknessCTL, 2), nanmean(ThicknessADHD, 2)]);
legend({'CTL', 'ADHD'}, 'Location', 'Best');
subplot 222;
plot_pos_neg_p_values(ObservedPValues, PermTSignValues, MIDX, {'CTL', 'ADHD'});
title('Observed p-values');

subplot 223;
plot_pos_neg_p_values(PermPValues, PermTSignValues, MIDX, {'CTL', 'ADHD'});
title('Permutation p-values');

subplot 224;
plot_pos_neg_p_values(ObservedPFDRValues, PermTSignValues, MIDX, {'CTL', 'ADHD'});
title('FDR p-values');
%title([num2str(size(ThicknessA, 2)) ' vs. ' num2str(size(ThicknessB, 2))]);
FigPos = fullscreen_fig_pos;
set(gcf, 'PaperPosition', FigPos, 'PaperUnits', 'points', 'InvertHardCopy', 'off');
OutputFile = 'perm_adhd_vs_ctl';
exportfig(gcf, OutputFile, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

IMG = imread([OutputFile, '.png']);
%IMG = imcomplement(IMG);
IMG = imautocropwhite(IMG, 20);
%IMG = imcomplement(IMG);
imwrite(IMG, [OutputFile, '.png']);
delete(gcf);