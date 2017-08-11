clear;

D = dir(fullfile('oasis_database', 'flipped', 'OAS*_msp.nii.gz'));

OutputDir = 'oasis_output_weighted_edge';
Bad = [2    61    69   101   115   120   131   183   209];
IDX = 1;
z = Bad(IDX);

BaseName = strrep(D(z).name, '.nii.gz', '');
CCName = strrep(D(z).name, 'msp', 'cc');

cc_seg_one_subject_final(fullfile('oasis_database', 'flipped', D(z).name), ...
	fullfile('oasis_database', 'flipped', D(z).name), ...
	OutputDir, ...
	[sprintf('%04d', z) '-' BaseName]);

%delete(gcf);
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