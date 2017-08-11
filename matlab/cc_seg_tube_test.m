clear;

% this tests the misalignment between the top and bottom of the CC
P = repmat(linspace(0, 0.05, 50), 1, 2)';
PSign = cat(2, ones(1, 50), -ones(1, 50))';


cc_seg_tube_p_display(P, PSign, {'Group 1', 'Group 2'});