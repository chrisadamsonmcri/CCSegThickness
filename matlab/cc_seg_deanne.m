clear;

InputDir = 'RawData';
D = dir(fullfile(InputDir, '*.nii.gz'));
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

Mode = 'stats';

IDX = 1:length(D);
Thicknesses = cell(1, length(IDX));
AllBaseNames = cell(1, length(IDX));
for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	AllBaseNames{z} = BaseName;
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	switch(lower(Mode))
		case 'seg_and_thickness'
		% 	if(exist([OutputBase '_seg.mat'], 'file') ~= 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName], 'acpcdetect');
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
		case 'binmask_images_thickness'
			disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
			% this is where the input images are already masks, so we just have to make a MAT file with the relevant variables
			SegReturnCode = cc_seg_one_subject_pipe_binmask(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName]);
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
	end
	T = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_thickness']));
	Thicknesses{z} = T.Thickness;
	%keyboard;
end
%delete(gcf);

Thicknesses = cat(2, Thicknesses{:});
NumNodes = size(Thicknesses, 1);
[tokens] = regexp(AllBaseNames, '^(NMR\d\d\d_)?(V\d\d\d)_cc_(\d)yr$', 'tokens');

NMRCodes = cell(1, length(AllBaseNames));
Years = zeros(1, length(AllBaseNames));

for z = 1:length(AllBaseNames)
	NMRCodes{z} = tokens{z}{1}{2};
	Years(z) = str2double(tokens{z}{1}{3});
end

[UniqueNMRCodes, ~, J] = unique(NMRCodes);

ZeroSevenYrThicknesses = NaN(size(Thicknesses, 1), length(UniqueNMRCodes), 2);

for z = 1:length(J)
	if(Years(z) == 0)
		ZeroSevenYrThicknesses(:, J(z), 1) = Thicknesses(:, z);
	else
		ZeroSevenYrThicknesses(:, J(z), 2) = Thicknesses(:, z);
	end
end
%keyboard;
%BothTimePoints = find(histc(J, 1:length(UniqueNMRCodes)) > 1);

FID = fopen('CClongitudinal4Chris.csv', 'r');

% get rid of the header line
tline = fgetl(FID);

CSVData.NMRCodes = cell(1, length(NMRCodes));
CSVData.ZeroYrExclude = false(1, length(NMRCodes));
CSVData.IsPreTerm = false(1, length(NMRCodes));
CSVData.ZeroYrICV = NaN(1, length(NMRCodes));
CSVData.SevenYrICV = NaN(1, length(NMRCodes));

IDX = 1;

while 1
    tline = fgetl(FID);
	if(~ischar(tline))
		break;
	end
	[tokens] = regexp(tline, '^(\d+),([01])?,([01])?,([01])?,([01])?,(\d+(\.\d+)?)?,(\d+(\.\d+)?)?$', 'tokens');
	CSVData.NMRCodes{IDX} = ['V', sprintf('%03d', str2double(tokens{1}{1}))];
	CSVData.ZeroYrExclude(IDX) = ~strcmp(tokens{1}{4}, '1');
	CSVData.IsPreTerm(IDX) = strcmp(tokens{1}{5}, '1');
	if(~isempty(tokens{1}{6}))
		CSVData.ZeroYrICV(IDX) = str2double(tokens{1}{6});
	end
	if(CSVData.ZeroYrExclude(IDX))
		CSVData.ZeroYrICV(IDX) = NaN;
	end
	
	if(~isempty(tokens{1}{7}))
		CSVData.SevenYrICV(IDX) = str2double(tokens{1}{7}) ./ 1000;
	end
	IDX = IDX + 1;
end
CSVData.NMRCodes = CSVData.NMRCodes(1:IDX - 1);
CSVData.ZeroYrExclude = CSVData.ZeroYrExclude(1:IDX - 1);
CSVData.IsPreTerm = CSVData.IsPreTerm(1:IDX - 1);
CSVData.ZeroYrICV = CSVData.ZeroYrICV(1:IDX - 1);
CSVData.SevenYrICV = CSVData.SevenYrICV(1:IDX - 1);
fclose(FID);

I = ~(isnan(CSVData.ZeroYrICV) & isnan(CSVData.SevenYrICV));

CSVData.NMRCodes = CSVData.NMRCodes(I);
CSVData.ZeroYrExclude = CSVData.ZeroYrExclude(I);
CSVData.IsPreTerm = CSVData.IsPreTerm(I);
CSVData.ZeroYrICV = CSVData.ZeroYrICV(I);
CSVData.SevenYrICV = CSVData.SevenYrICV(I);

clear I;

[I, LOC] = ismember(CSVData.NMRCodes, UniqueNMRCodes);

CSVData.NMRCodes = CSVData.NMRCodes(I);
CSVData.ZeroYrExclude = CSVData.ZeroYrExclude(I);
CSVData.IsPreTerm = CSVData.IsPreTerm(I);
CSVData.ZeroYrICV = CSVData.ZeroYrICV(I);
CSVData.SevenYrICV = CSVData.SevenYrICV(I);


UniqueNMRCodes = UniqueNMRCodes(LOC(I));
%Thicknesses = Thicknesses(:, LOC);
ZeroSevenYrThicknesses = ZeroSevenYrThicknesses(:, LOC(I), :);

CSVData.ZeroYrICV(isnan(ZeroSevenYrThicknesses(1, :, 1))) = NaN;
CSVData.SevenYrICV(isnan(ZeroSevenYrThicknesses(1, :, 2))) = NaN;

TwoTimePoints = ~(isnan(CSVData.ZeroYrICV) | isnan(CSVData.SevenYrICV));

TwoSampleThicknessA = cell(1, 3);
TwoSampleThicknessB = cell(1, 3);
TwoSampleTitles = {'0yr PT vs. 0yr FT', '7yr PT vs. 7yr FT', '7yr FT - 0yr FT vs. 7yr PT - 0yr PT'};
TwoSampleGroupNames = regexp(TwoSampleTitles, '^(.+) vs\. (.+)$', 'tokens');

% Comparisons two-sample Welch's T-test
%	0yr PT vs. 0yr FT
% ThicknessA is 0yr PT
TwoSampleThicknessA{1} = ZeroSevenYrThicknesses(:,  CSVData.IsPreTerm & ~CSVData.ZeroYrExclude & ~isnan(CSVData.ZeroYrICV), 1); % PT
TwoSampleThicknessB{1} = ZeroSevenYrThicknesses(:, ~CSVData.IsPreTerm & ~CSVData.ZeroYrExclude & ~isnan(CSVData.ZeroYrICV), 1); % FT

%	7yr PT vs. 7yr FT
TwoSampleThicknessA{2} = ZeroSevenYrThicknesses(:,  CSVData.IsPreTerm & ~isnan(CSVData.SevenYrICV), 2); % PT
TwoSampleThicknessB{2} = ZeroSevenYrThicknesses(:, ~CSVData.IsPreTerm & ~isnan(CSVData.SevenYrICV), 2); % FT

%	7yr FT - 0yr FT vs 7yr PT - 0yr PT
I = TwoTimePoints & ~CSVData.IsPreTerm;
TwoSampleThicknessA{3} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 2), CSVData.SevenYrICV(I)) - bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 1), CSVData.ZeroYrICV(I));
I = TwoTimePoints & CSVData.IsPreTerm;
TwoSampleThicknessB{3} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 2), CSVData.SevenYrICV(I)) - bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 1), CSVData.ZeroYrICV(I));

OneSampleThicknessA = cell(1, 3);
OneSampleThicknessB = cell(1, 3);
OneSampleTitles = {'0yr vs. 7yr', '0yr FT vs. 7yr FT', '0yr PT vs. 7yr PT'};
OneSampleGroupNames = regexp(OneSampleTitles, '^(.+) vs\. (.+)$', 'tokens');
% One sample Student's T-test
%	0yr vs. 7yr
OneSampleThicknessA{1} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, TwoTimePoints, 1), CSVData.ZeroYrICV(TwoTimePoints));
OneSampleThicknessB{1} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, TwoTimePoints, 2), CSVData.SevenYrICV(TwoTimePoints));
%	0yr FT vs. 7yr FT
I = TwoTimePoints & ~CSVData.IsPreTerm;
OneSampleThicknessA{2} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 1), CSVData.ZeroYrICV(I));
OneSampleThicknessB{2} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 2), CSVData.SevenYrICV(I));
%	0yr PT vs. 7yr PT
I = TwoTimePoints & CSVData.IsPreTerm;
OneSampleThicknessA{3} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 1), CSVData.ZeroYrICV(I));
OneSampleThicknessB{3} = bsxfun(@rdivide, ZeroSevenYrThicknesses(:, I, 2), CSVData.SevenYrICV(I));

W = 10;
IDX = W:(NumNodes - W);

[PValueTwoSample, ObservedTSignTwoSample, omniPTwoSample, ObservedPTwoSample, ObservedPFDRTwoSample] = deal(cell(1, length(TwoSampleThicknessA)));

for z = 1:length(TwoSampleThicknessA)
	disp(['Two sample test ' num2str(z), ': ' TwoSampleTitles{z}]);
	disp(['Group A subjects: ' num2str(size(TwoSampleThicknessA{z}, 2)) ', Group B subjects: ' num2str(size(TwoSampleThicknessB{z}, 2))]);
	[PValueTwoSample{z}, ...
		ObservedTSignTwoSample{z}, ...
		omniPTwoSample{z}, ...
		ObservedPTwoSample{z}] = cc_seg_twosampleT_perm(TwoSampleThicknessA{z}(IDX, :), TwoSampleThicknessB{z}(IDX, :), 100000);
	[~, ~, ObservedPFDRTwoSample{z}] = fdr(ObservedPTwoSample{z}, 0.05);
% 	subplot(SR, SC, z);
% 
% 	%plot([mean(TwoSampleThicknessA{z}(IDX, :), 2), mean(TwoSampleThicknessB{z}(IDX, :), 2)]);
% 	title(TwoSampleTitles{z});
% 	legend({'A < B', 'B < A'});
end

[PValueOneSample, ObservedTSignOneSample, omniPOneSample, ObservedPOneSample, ObservedPFDROneSample] = deal(cell(1, length(OneSampleThicknessA)));

for z = 1:length(OneSampleThicknessA)
	disp(['One sample test ' num2str(z), ': ' OneSampleTitles{z}]);
	disp(['Group A subjects: ' num2str(size(OneSampleThicknessA{z}, 2)) ', Group B subjects: ' num2str(size(OneSampleThicknessB{z}, 2))]);
	[PValueOneSample{z}, ...
		ObservedTSignOneSample{z}, ...
		omniPOneSample{z}, ...
		ObservedPOneSample{z}] = cc_seg_pairedT_perm(OneSampleThicknessA{z}(IDX, :), OneSampleThicknessB{z}(IDX, :), 100000);
	
% 	subplot(SR, SC, z + 3);
% 	plot_pos_neg_p_values(PValueOneSample{z}, ObservedTSignOneSample{z}, IDX);
% 	%plot([mean(OneSampleThicknessA{z}(IDX, :), 2), mean(OneSampleThicknessB{z}(IDX, :), 2)]);
% 	title(OneSampleTitles{z});
% 	legend({'A < B', 'B < A'});
end
%%
for z = 1:length(TwoSampleThicknessA)
	[~, ~, ObservedPFDRTwoSample{z}] = fdr(ObservedPTwoSample{z}, 0.05);
end

for z = 1:length(OneSampleThicknessA)
	[~, ~, ObservedPFDROneSample{z}] = fdr(ObservedPOneSample{z}, 0.05);
end


	
SR = 2; SC = 3;


figure;
for z = 1:length(TwoSampleThicknessA)
	subplot(SR, SC, z);
	plot_pos_neg_p_values(ObservedPTwoSample{z}, ObservedTSignTwoSample{z}, IDX, TwoSampleGroupNames{z}{1});
	title(TwoSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

for z = 1:length(OneSampleThicknessA)
	subplot(SR, SC, z + 3);
	plot_pos_neg_p_values(ObservedPOneSample{z}, ObservedTSignOneSample{z}, IDX, OneSampleGroupNames{z}{1});
	%plot([mean(OneSampleThicknessA{z}(IDX, :), 2), mean(OneSampleThicknessB{z}(IDX, :), 2)]);
	title(OneSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

figure;
for z = 1:length(TwoSampleThicknessA)
	subplot(SR, SC, z);
	plot_pos_neg_p_values(PValueTwoSample{z}, ObservedTSignTwoSample{z}, IDX, TwoSampleGroupNames{z}{1});
	title(TwoSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

for z = 1:length(OneSampleThicknessA)
	subplot(SR, SC, z + 3);
	plot_pos_neg_p_values(PValueOneSample{z}, ObservedTSignOneSample{z}, IDX, OneSampleGroupNames{z}{1});
	%plot([mean(OneSampleThicknessA{z}(IDX, :), 2), mean(OneSampleThicknessB{z}(IDX, :), 2)]);
	title(OneSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

figure;
for z = 1:length(TwoSampleThicknessA)
	subplot(SR, SC, z);
	plot_pos_neg_p_values(ObservedPFDRTwoSample{z}, ObservedTSignTwoSample{z}, IDX, TwoSampleGroupNames{z}{1});
	title(TwoSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

for z = 1:length(OneSampleThicknessA)
	subplot(SR, SC, z + 3);
	plot_pos_neg_p_values(ObservedPFDROneSample{z}, ObservedTSignOneSample{z}, IDX, OneSampleGroupNames{z}{1});
	%plot([mean(OneSampleThicknessA{z}(IDX, :), 2), mean(OneSampleThicknessB{z}(IDX, :), 2)]);
	title(OneSampleTitles{z});
	%legend({'A > B', 'B > A'});
end

%%
figure;
for z = 1:length(TwoSampleThicknessA)
	subplot(SR, SC, z);
	plot([mean(TwoSampleThicknessA{z}(IDX, :), 2), mean(TwoSampleThicknessB{z}(IDX, :), 2)]);
	title(TwoSampleTitles{z});
	legend(TwoSampleGroupNames{z}{1});
end

for z = 1:length(OneSampleThicknessA)
	subplot(SR, SC, z + 3);
	plot([mean(OneSampleThicknessA{z}(IDX, :), 2), mean(OneSampleThicknessB{z}(IDX, :), 2)]);
	title(OneSampleTitles{z});
	legend(OneSampleGroupNames{z}{1});
end
%%
% 
% ThicknessA = cell(1, length(BothTimePoints));
% ThicknessB = cell(1, length(BothTimePoints));
% 
% for z = 1:length(BothTimePoints)
% 	L = find(J == BothTimePoints(z));
% 	if(Years(L(1)) == 0)
% 		ThicknessA{z} = Thicknesses(:, L(1));
% 		ThicknessB{z} = Thicknesses(:, L(2));
% 	else
% 		ThicknessA{z} = Thicknesses(:, L(2));
% 		ThicknessB{z} = Thicknesses(:, L(1));
% 	end
% 	clear L;
% end
% 
% ThicknessA = cat(2, ThicknessA{:});
% ThicknessB = cat(2, ThicknessB{:});
%plot([mean(ThicknessA, 2), mean(ThicknessB, 2)]);
%[PValue, ObservedTSign, omniP] = cc_seg_pairedT_perm(ThicknessA(IDX, :), ThicknessB(IDX, :), 100000);

%[PValueTwoSample, ObservedTSignTwoSample, omniPTwoSample, ObservedP] = cc_seg_twosampleT_perm(Thicknesses(IDX, Years == 0), Thicknesses(IDX, Years == 7), 100000);


%keyboard;

% FID = fopen('/data/addo/cc_seg/deanne_cc_longitudinal/thicknesses.csv', 'w');
% fprintf(FID, [repmat('%s,', 1, length(AllBaseNames) - 1) '%s\n'], AllBaseNames{:});
% fprintf(FID, [repmat('%f,', 1, length(AllBaseNames) - 1) '%f\n'], Thicknesses');
% fclose(FID);