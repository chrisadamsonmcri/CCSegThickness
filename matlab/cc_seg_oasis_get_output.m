clear;

D = dir(fullfile('oasis_database', 'flipped', 'OAS*_msp.nii.gz'));
%FigureDir = 'oasis_figures_weighted_edge_exclusion_normxcorr_reg';
%[~, ~, ~] = mkdir(FigureDir);
OutputDir = 'oasis_output_weighted_edge';
%FID = fopen('cc_seg_oasis.parallel', 'w');
%cc_seg_one_subject('0001_acpc.nii.gz', 'test.png');


for z = 1:length(D)
	disp(['Subject ' num2str(z)]);
	BaseName = strrep(D(z).name, '.nii.gz', '');
	CCName = strrep(D(z).name, 'msp', 'cc');
	
	Output(z) = load(fullfile(OutputDir, [sprintf('%04d', z) '-' BaseName '.mat']));
end
%fclose(FID);
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