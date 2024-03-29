clear;

D = dir(fullfile('datasets', 'mindboggle101_N4', '*.nii.gz'));

% take out the colin-27 image, it is weird and not a standard T1
D = D([1, 3:end]);
FigureDir = 'mindboggle_N4_test';
[~, ~, ~] = mkdir(FigureDir);
OutputDir = FigureDir;
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');

%ConvexHullSubjects = [2,14,27,40,61,67,69,75,77,79,95,106,115,116,120,144,152,183,199,209,234,239];

ConvexHullUsed = false(1, length(D));

%Output
for z = 1:length(D)
	disp(['Subject ' num2str(z) ' of ' num2str(length(D)) ', ' D(z).name]);
	BaseName = strrep(D(z).name, '.nii.gz', '');
	CCName = strrep(D(z).name, 'msp', 'cc');
	%cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), fullfile('oasis_database', 'flipped', CCName), fullfile(FigureDir, [sprintf('%04d', z) '-' BaseName]));
	
	%if(exist(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '.mat']), 'file') ~= 2)
		cc_seg_one_subject_pipe(fullfile('datasets', 'mindboggle101_N4', D(z).name), ...
		[], ...
		OutputDir, ...
		[sprintf('%04d', z) '-' BaseName]);
	%end
	%Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName]));
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