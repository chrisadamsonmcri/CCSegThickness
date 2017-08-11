clear;

InputDir = 'RawT1ReorientCropped';
D = dir(fullfile(InputDir, '*.nii.gz'));

[BaseNames{1:length(D)}] = deal(D.name);
BaseNames = strrep(BaseNames, '.nii.gz', '');
SubjectIDs = strrep(BaseNames, '_sc', '');

%FigureDir = 'cc_seg_flirt';
FigureDir = 'cc_seg_python_mat';
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

% for z = IDX
% 	BaseName = strrep(D(z).name, '.nii.gz', '');
% 	OutputBase = fullfile(OutputDir, BaseName);
% 	
% 	InputHDF5File = [OutputBase '_thickness.hdf5'];
% 	
% 	T.Thickness = h5read(InputHDF5File, '/thicknessProfile');
% 	T.ValidStreamlines = h5read(InputHDF5File, '/validStreamlines');
% 	T.ValidStreamlines = (T.ValidStreamlines > 0);
% 	
% 	T.Streamlines = cell(numel(T.ValidStreamlines), 1);
% 	A = zeros(numel(T.ValidStreamlines), 1);
% 	for CurStreamline = 1:numel(T.ValidStreamlines)
% 		T.Streamlines{CurStreamline} = h5read(InputHDF5File, ['/streamlines/' num2str(CurStreamline - 1)]);
% 		A(CurStreamline) = arc_length(T.Streamlines{CurStreamline});
% 	end
% 	
% 	%keyboard;
% 	
% 	%T = load([OutputBase '_thickness.mat']);
% 	
% 	S = T.Thickness;
% 	S(~T.ValidStreamlines) = NaN;
% 	ThicknessProfiles(:, z) = S;
% 	clear S T;
% end

% FID = fopen('thickness_profiles.csv', 'w');
% fprintf(FID, '%s', BaseNames{1});
% fprintf(FID, ',%s', BaseNames{2:end});
% fprintf(FID, '\n');
% 
% fprintf(FID, ['%f' repmat(',%f', 1, size(ThicknessProfiles, 2) - 1) '\n'], ThicknessProfiles');
% fclose(FID);

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

[GroupIDX, GroupLabels] = grp2idx(ADNIData.dx);

ADNIData.PTID = cellfun(@(x) ([x '_sc']), ADNIData.PTID, 'UniformOutput', false);

for GroupI = 1:length(GroupLabels) - 1
	for GroupJ = GroupI + 1:length(GroupLabels)
		FID = fopen(['paired_groups_' GroupLabels{GroupI} '_' GroupLabels{GroupJ} '.txt'], 'w');
		T = cat(1, repmat(GroupLabels(GroupI), 1, sum(GroupIDX == GroupI)), ADNIData.PTID(GroupIDX == GroupI)');
		fprintf(FID, '%s\t%s\n', T{:});
		T = cat(1, repmat(GroupLabels(GroupJ), 1, sum(GroupIDX == GroupJ)), ADNIData.PTID(GroupIDX == GroupJ)');
		fprintf(FID, '%s\t%s\n', T{:});
		clear T;
		fclose(FID);
		%keyboard;
	end
end