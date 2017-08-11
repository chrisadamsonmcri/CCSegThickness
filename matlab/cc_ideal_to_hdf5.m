clear;

OutputHDF5 = 'ideal_cc.hdf5';
delete(OutputHDF5);

S = load('ideal_callosum_contours.mat');
F = fieldnames(S);

for z = 1:length(F)
	h5create(OutputHDF5, ['/contours/' F{z}], size(S.(F{z})));
	h5write(OutputHDF5, ['/contours/' F{z}], S.(F{z}));
end

S = load('ideal_callosum_interp_points_curved.mat');
F = fieldnames(S);

for z = 1:length(F)
	h5create(OutputHDF5, ['/' F{z}], size(S.(F{z})));
	h5write(OutputHDF5, ['/' F{z}], S.(F{z}));
end