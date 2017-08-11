clear;

InputDir = 'InputDir';
%D = dir(fullfile(InputDir, '*.nii.gz'));

%[BaseNames{1:length(D)}] = deal(D.name);
%BaseNames = strrep(BaseNames, '.nii.gz', '');

FID = fopen('dwi_good_ids.txt', 'r');
BaseNames = textscan(FID, '%s');
BaseNames = BaseNames{1};
fclose(FID);
%return;

%FigureDir = 'cc_seg_flirt';
OutputDir = 'OutputDir';
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

%Mode = 'seg_and_thickness';
Mode = 'skip';
%Mode = 'seg';
%Mode = 'manedit_thickness';

IDX = 1:length(BaseNames);

%[~, IDX] = ismember({'0051036_T1_mprage', '0051322_T1_mprage'}, BaseNames);
%keyboard;
%IDX = setdiff(IDX, [275, 353]);
% IDX = 762:length(D);
% IDX = 761;
% % flirt
% IDX = [226, 254, 340, 519, 530, 582, 592, 665, 674, 711, 756];
% IDX2 = [226, 519, 756];
% IDX = setdiff(IDX, IDX2);
% IDX = 592;
%IDX = 672;
% IDX = 1;

% IDs = cell(1, length(BaseNames));
% TimePoints = zeros(1, length());
CurSubject = 1;
%BaseName = strrep(D(z).name, '.nii.gz', '');

FID = fopen('midsag_matlab_redo.txt', 'w');
for z = IDX

	BaseName = BaseNames{z};
	%OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	OutputBase = fullfile(OutputDir, BaseName);
	switch(lower(Mode))
		case 'midsag'
			cc_seg_one_subject_pipe_midsag(fullfile(InputDir, BaseNames{z}), ...
				OutputDir, ...
				BaseName, 'flirt');
			
		case 'seg_and_thickness'
		% 	if(exist([OutputBase '_seg.mat'], 'file') ~= 2)
				disp(['Subject ' num2str(CurSubject) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, BaseNames{z}), ...
				[], ...true
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
			%if(exist([OutputBase '_seg.mat'], 'file') ~= 2 || exist([OutputBase '_thickness.mat'], 'file') ~= 2 )
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(OutputBase, [], true);
				ThicknessReturnCode = 0;
				if(SegReturnCode == 0)
					ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
				end
				if(SegReturnCode == 0 || ThicknessReturnCode == 0)
					fprintf(FID, '%s\n', BaseName);
				end
			%end
		case 'thickness'
			disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
			ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
		case 'manedit_thickness'
			if(exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
				ManeditAttribs = dir([OutputBase '_seg_manedit.mat']);
				if(exist([OutputBase '_thickness.mat'], 'file') == 2)
					ThicknessAttribs = dir([OutputBase '_thickness.mat']);
					%keyboard;
				else
					ThicknessAttribs.datenum = -Inf;
				end
				if(ManeditAttribs.datenum > ThicknessAttribs.datenum)
					disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
					ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
				end
				
				%keyboard;
				
				
			end
	end
	CurSubject = CurSubject + 1;
	%
end
delete(gcf);
fclose(FID);

NumNodes = 100;

ThicknessProfiles = NaN(NumNodes, length(BaseNames));

for z = 1:length(BaseNames)
	BaseName = BaseNames{z};
	OutputBase = fullfile(OutputDir, BaseName);
	if(exist([OutputBase '_thickness.mat'], 'file') == 2)
		T = load([OutputBase '_thickness.mat'], 'Thickness', 'ValidStreamlines');
		ThicknessProfiles(:, z) = T.Thickness(:);
		ThicknessProfiles(~T.ValidStreamlines, z) = NaN;
		%keyboard;
	end
end

FID = fopen('thickness_profiles.csv', 'w');
fprintf(FID, '%s', BaseNames{1});
fprintf(FID, ',%s', BaseNames{2:end});
fprintf(FID, '\n');

fprintf(FID, ['%f' repmat(',%f', 1, NumNodes - 1) '\n'], ThicknessProfiles');
fclose(FID);