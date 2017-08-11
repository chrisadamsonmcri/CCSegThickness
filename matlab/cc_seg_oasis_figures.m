clear;

D = dir(fullfile('oasis_database', 'flipped', 'OAS*_msp.nii.gz'));
D = dir(fullfile('oasis_msp', 'OAS*_msp.nii.gz'));
FigureDir = 'oasis_lkccseg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

ConvexHullSubjects = [2,14,27,40,61,67,69,75,77,79,95,106,115,116,120,144,152,183,199,209,234,239];

ConvexHullUsed = false(1, length(D));

%Output
IDX = 1:length(D);
%IDX = 132;%:length(D);
%IDX = 144;
%IDX = 1; % this has an artefact
%IDX = 57;
% for vein removal
% VeinRemovalFiles = dir(fullfile(FigureDir, 'paper_figures', '*VeinRemoval*.png'));
% [VeinRemovalFileNames{1:length(VeinRemovalFiles)}] = deal(VeinRemovalFiles.name);
% 
% tokens = regexp(VeinRemovalFileNames, '^(\d\d\d\d)-.*$', 'tokens');
% VeinRemovalNumbers = zeros(1, length(VeinRemovalFiles));
% 
% for z = 1:length(VeinRemovalFiles)
% 	VeinRemovalNumbers(z) = str2double(tokens{z}{1});
% end
VeinRemovalNumbers = [6 9 12 24 29 33 36 37 38 48 54 56 57 61 62 69 71 74 76 77 80 85 88 89 94 100 101 107 117 146 150 153 155 159 163 183 192 198 201 211 213 221 223 226 231 242 258 274 275 276 280 282 284 286 290 293 303 305 307 311];
%keyboard;

%FigureIDX = [1, 2, 12, 57, 144, 311];
%IDX = FigureIDX;
%for z = 57%VeinRemovalNumbers
%FID = fopen('cc_seg_oasis.parallel', 'w');

for z = IDX
	
	BaseName = strrep(D(z).name, '.nii.gz', '');
	CCName = strrep(D(z).name, 'msp', 'cc');
	%cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), fullfile('oasis_database', 'flipped', CCName), fullfile(FigureDir, [sprintf('%04d', z) '-' BaseName]));
	
	%if(exist(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '.mat']), 'file') ~= 2)
		disp(['Subject ' num2str(z)]);
%   		fprintf(FID, 'matlab -nosplash -nodesktop -r "cc_seg_one_subject_pipe_seg_figures(''%s'', ''%s'', ''%s'', ''%s'', ''%s'', %d); exit;"\n', ...
%  			fullfile('oasis_database', 'flipped', D(z).name), ...
%  			fullfile('oasis_database', 'flipped', CCName), ...
%  			OutputDir, ...
%  			[sprintf('%04d', z) '-' BaseName], ...
%  			'acpcdetect', ...
%  			z == 1); ... this is the angle figure, we only do it for the first subject
% 		
		cc_seg_one_subject_pipe_seg_figures(fullfile('oasis_database', 'flipped', D(z).name), ...
 		fullfile('oasis_database', 'flipped', CCName), ...
 		OutputDir, ...
 		[sprintf('%04d', z) '-' BaseName], 'acpcdetect', z == 1); ... this is the angle figure, we only do it for the first subject
		%fprintf(FID, 'matlab -nosplash -nodesktop -r "cc_seg_one_subject_pipe_thickness_figures(''%s''); exit;"\n', ...
		%fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]));
		%cc_seg_one_subject_pipe_thickness_figures(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]));
	%end
	%Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
end
%fclose(FID);
delete(gcf);