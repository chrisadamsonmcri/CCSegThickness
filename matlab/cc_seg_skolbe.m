clear;

InputDir = 'RawData';
D = dir(fullfile(InputDir, '*.nii.gz'));
FigureDir = 'cc_seg';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

Mode = 'seg_only';

IDX = 1:length(D);
%IDX = 173;
for z = IDX
	BaseName = strrep(D(z).name, '.nii.gz', '');
	OutputBase = fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]);
	
	switch(lower(Mode))
		case 'seg_only'
			disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				SegReturnCode = cc_seg_one_subject_pipe_seg(fullfile(InputDir, D(z).name), ...
				[], ...
				OutputDir, ...
				[sprintf('%04d', z) '-' BaseName], 'flirt');
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
		case 'manedit_thickness'
			if(exist([OutputBase '_seg_manedit.mat'], 'file') == 2)
				disp(['Subject ' num2str(z) ' of ' num2str(length(IDX)) ', ' BaseName]);
				ThicknessReturnCode = cc_seg_one_subject_pipe_thickness(OutputBase);
			end
	end
	%Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '_seg']));
end
delete(gcf);
% [~, I] = min([Output.FinalDice]);
% hold off;
% imshow(Output(I).IMG, []);
% hold on;
% [~, CC] = contour(Output(I).GroundSeg, [0.5, 0.5]);
% set(CC, 'Color', 'b');
% [~, CC] = contour(Output(I).FinalSeg, [0.5, 0.5]);
% set(CC, 'Color', 'r');
% axis equal ij;
%delete(gcf);
% 
% FID = fopen('convex_hull_used.txt', 'w');
% fprintf(FID, 'Used convex hull:\n');
% fprintf(FID, '%d\n', find(ConvexHullUsed));
% fprintf(FID, 'Didnt use convex hull:\n');
% fprintf(FID, '%d\n', find(~ConvexHullUsed));
% fclose(FID);

% for IDX = 1:length(ConvexHullSubjects)
% 	z = ConvexHullSubjects(IDX);
% 	disp(['Subject ' num2str(z)]);
% 	BaseName = strrep(D(z).name, '.nii.gz', '');
% 	CCName = strrep(D(z).name, 'msp', 'cc');
% 	%cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), fullfile('oasis_database', 'flipped', CCName), fullfile(FigureDir, [sprintf('%04d', z) '-' BaseName]));
% 	ConvexHullUsed(z) = cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), ...
% 	fullfile('oasis_database', 'flipped', D(z).name), ...
% 	OutputDir, ...
% 	[sprintf('%04d', z) '-' BaseName]);
% end
% delete(gcf);

%return;
% MeanOriginal = mean([Overlap.Original]);
% MeanOriginalW = mean([Overlap.OriginalW]);
% MeanSeg = nanmean([Overlap.Seg]);
% MeanSegW = nanmean([Overlap.SegW]);
% 
% MinOriginal = min([Overlap.Original]);
% MinOriginalW = min([Overlap.OriginalW]);
% MinSeg = nanmin([Overlap.Seg]);
% MinSegW = nanmin([Overlap.SegW]);
% 
% MaxOriginal = max([Overlap.Original]);
% MaxOriginalW = max([Overlap.OriginalW]);
% MaxSeg = nanmax([Overlap.Seg]);
% MaxSegW = nanmax([Overlap.SegW]);
% 
% disp(['Original: ' num2str(MeanOriginal) ' [' num2str(MinOriginal) ', ' num2str(MaxOriginal) ']']);
% disp(['OriginalW: ' num2str(MeanOriginalW) ' [' num2str(MinOriginalW) ', ' num2str(MaxOriginalW) ']']);
% disp(['Seg: ' num2str(MeanSeg) ' [' num2str(MinSeg) ', ' num2str(MaxSeg) ']']);
% disp(['SegW: ' num2str(MeanSegW) ' [' num2str(MinSegW) ', ' num2str(MaxSegW) ']']);
% 
% save lk_overlap Overlap Mean* Min* Max*;