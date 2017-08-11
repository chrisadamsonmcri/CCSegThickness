function [ReturnCode] = cc_seg_one_subject_pipe_thickness(MATFilePrefix, NumThicknessNodes)

% cc_seg_one_subject_pipe_thickness(MATFile)
%
% DESCRIPTION
%	Performs the thickness computation part of the CC pipeline.
%
% PARAMETERS
%	MATFilePrefix (char): the base of the output file created by
%	cc_seg_one_subject_pipe_seg, usually run by 
%	cc_seg_process(..., ..., 'seg') or
%	cc_seg_process(..., ..., 'seg_and_thickness') or
%	NumThicknessNodes [1]: the number of streamlines (nodes) to generate,
%	100 by default
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
%

if(nargin < 2)
	NumThicknessNodes = 100;
end

if(~isnumeric(NumThicknessNodes) || ~isscalar(NumThicknessNodes))
	ReturnCode = 1;
	return;
end
if(NumThicknessNodes < 1 || floor(NumThicknessNodes) ~= NumThicknessNodes)
	ReturnCode = 1;
	return;
end

ReturnCode = 0;

[pathstr, name, ext] = fileparts(MATFilePrefix);
[~, ~, ~] = mkdir(fullfile(pathstr, 'endpoint'));
EndpointFile = fullfile(pathstr, 'endpoint', name);
disp(EndpointFile);

if(strcmp(ext, '.mat'))
	FileName = name;
else
	FileName = [name ext];
end

ManEdit = [];

SegMatFileNamePrefix = fullfile(pathstr, [FileName '_seg']);
SegManeditMatFileNamePrefix = fullfile(pathstr, [FileName '_seg_manedit']);
ThicknessMatFileName = fullfile(pathstr, [FileName '_thickness']);

if(exist([SegMatFileNamePrefix '.mat'], 'file') == 2)
	SegMatFileName = [SegMatFileNamePrefix '.mat'];
	SegMatFileExt = 'mat';
elseif(exist([SegMatFileNamePrefix '.hdf5'], 'file') == 2)
	SegMatFileName = [SegMatFileNamePrefix '.hdf5'];
	SegMatFileExt = 'hdf5';
else
	error('Seg file not found');
end

switch(SegMatFileExt)
	case 'mat'
		Seg = load(SegMatFileName);
	case 'hdf5'
		%keyboard;
		
		try
			I = h5info(SegMatFileName);
			clear I;
			HDFReadFunc = @h5read;
		catch e
			HDFReadFunc = @hdf5read;
		end
		
		Seg.IMG = HDFReadFunc(SegMatFileName, '/IMG'); Seg.IMG = Seg.IMG';
		Seg.FinalSeg = HDFReadFunc(SegMatFileName, '/finalSeg'); Seg.FinalSeg = logical(Seg.FinalSeg');
		Seg.FinalSegArtefactsRemoved = HDFReadFunc(SegMatFileName, '/finalSegNoArtefacts'); Seg.FinalSegArtefactsRemoved = logical(Seg.FinalSegArtefactsRemoved');
		Seg.InitialSeg = HDFReadFunc(SegMatFileName, '/initialSeg'); Seg.InitialSeg = logical(Seg.InitialSeg');
		Seg.TemplatePixdims = HDFReadFunc(SegMatFileName, '/templatePixdims'); Seg.TemplatePixdims = Seg.TemplatePixdims';
end

if(exist([SegManeditMatFileNamePrefix '.mat'], 'file') == 2)
	SegManeditMatFileName = [SegManeditMatFileNamePrefix '.mat'];
	SegManeditMatFileExt = 'mat';
elseif(exist([SegManeditMatFileNamePrefix '.hdf5'], 'file') == 2)
	SegManeditMatFileName = [SegManeditMatFileNamePrefix '.hdf5'];
	SegManeditMatFileExt = 'hdf5';
else
	SegManeditMatFileName = [];
	SegManeditMatFileExt = [];
end

if(ischar(SegManeditMatFileName))
	switch(SegManeditMatFileExt)
	case 'mat'
		ManEdit = load(SegManeditMatFileName);
	case 'hdf5'
		%keyboard;
		try
			I = h5info(SegManeditMatFileName);
			clear I;
			HDFReadFunc = @h5read;
		catch e
			HDFReadFunc = @hdf5read;
		end
		ManEdit.FinalSegManEdit = hdf5read(SegMatFileName, '/finalSegManedit'); ManEdit.FinalSegManEdit = logical(ManEdit.FinalSegManEdit');
	end
else
	ManEdit = [];
end

InitialSegVariableName = 'FinalSegArtefactsRemoved';
ManEditVariableName = 'FinalSegManEdit';

if(~(isfield(Seg, InitialSegVariableName) && isfield(Seg, 'TemplatePixdims')))
	error('The seg file does not have the required fields');
else
	FinalSeg = Seg.(InitialSegVariableName);
end

if(~isempty(ManEdit))
	if(isfield(ManEdit, ManEditVariableName))
		FinalSeg = ManEdit.(ManEditVariableName);
	else
		disp('Manedit file does not have required field');
	end
end

GroundTruthGiven = isfield(Seg, 'GroundSeg');

ThicknessCalculationError = false;
try
	[Output.Contours, ...
		Output.SolvedImage, ...
		Output.InitialPoints, ...
		Output.ArcLengthsLeft, ...
		Output.ArcLengthsRight, ...
		Output.LeftOffsetIDXStart, ...
		Output.CurLeftOffsetIDX, ...
		Output.RightOffsetIDXStart, ...
		Output.CurRightOffsetIDX] = endpoint_find_one_subject(FinalSeg, EndpointFile);
	ReturnCode = 1;
	return;
	[Output.MaskInnerBoundary, ...
		Output.MaskOuterBoundary, ...
		Output.MaskFree, ...
		Output.SolvedImage, ...
		Output.Streamlines, ...
		Output.X, ...
		Output.Y, ...
		Output.NormFX, ...
		Output.NormFY, ...
		Output.ValidStreamlines, ...
		Output.StartV] = laplace_get_points_2d_auto_mw(Output.Contours.xi, Output.Contours.yi, Output.Contours.xo(end:-1:1), Output.Contours.yo(end:-1:1), 1, 1, 1, 100);

	StreamlinesForArcLengths = cell(size(Output.Streamlines));

	for z = 1:length(Output.Streamlines)
		StreamlinesForArcLengths{z} = bsxfun(@times, Output.Streamlines{z}, Seg.TemplatePixdims(:)');
	end
	
	T = cellfun(@arc_length, StreamlinesForArcLengths(Output.ValidStreamlines));
	Output.Thickness = zeros(size(Output.Streamlines));
	Output.Thickness(Output.ValidStreamlines) = T;
	
	
	clear T;
	% production release, remove variables that are intermediate
	Output = rmfield(Output, {'StartV', 'NormFX', 'NormFY', 'X', 'Y', 'MaskFree', 'MaskInnerBoundary', 'MaskOuterBoundary', ...
		'SolvedImage', ...
		'InitialPoints', ...
		'ArcLengthsLeft', ...
		'ArcLengthsRight', ...
		'LeftOffsetIDXStart', ...
		'CurLeftOffsetIDX', ...
		'RightOffsetIDXStart', ...
		'CurRightOffsetIDX'});
		
catch ME
	disp('There was a problem with the thickness calculations');
	%ThicknessCalculationError = true;
	disp(ME.message);
	disp('Stack:');
	for z = 1:length(ME.stack)
		disp([ME.stack(z).file '>' ME.stack(z).name ': line ' num2str(ME.stack(z).line)]);
	end
	ReturnCode = 1;
	return;
end

if(GroundTruthGiven)
	try
		ThicknessCalculationErrorGround = false;
		[Output.GND.Contours, ...
			Output.GND.SolvedImage, ...
			Output.GND.InitialPoints, ...
			Output.GND.ArcLengthsLeft, ...
			Output.GND.ArcLengthsRight, ...
			Output.GND.LeftOffsetIDXStart, ...
			Output.GND.CurLeftOffsetIDX, ...
			Output.GND.RightOffsetIDXStart, ...
			Output.GND.CurRightOffsetIDX] = endpoint_find_one_subject(imdilate(Seg.GroundSeg, strel('disk', 1)));
		[Output.GND.MaskInnerBoundary, ...
			Output.GND.MaskOuterBoundary, ...
			Output.GND.MaskFree, ...
			Output.GND.SolvedImage, ...
			Output.GND.Streamlines, ...
			Output.GND.X, ...
			Output.GND.Y, ...
			Output.GND.NormFX, ...
			Output.GND.NormFY, ...
			Output.GND.ValidStreamlines, ...
			Output.GND.StartV] = laplace_get_points_2d_auto_mw(Output.GND.Contours.xi, Output.GND.Contours.yi, Output.GND.Contours.xo(end:-1:1), Output.GND.Contours.yo(end:-1:1), 1, 1, 1, NumThicknessNodes);
			StreamlinesForArcLengths = cell(size(Output.GND.Streamlines));
		%keyboard;
		for z = 1:length(Output.GND.Streamlines)
			StreamlinesForArcLengths{z} = bsxfun(@times, Output.GND.Streamlines{z}, Seg.TemplatePixdims(:)');
		end

		T = cellfun(@arc_length, StreamlinesForArcLengths(Output.GND.ValidStreamlines));
		Output.GND.Thickness = zeros(size(Output.GND.Streamlines));
		Output.GND.Thickness(Output.GND.ValidStreamlines) = T;
		clear T;
	catch ME
		disp('There was a problem with the thickness calculations, ground truth image');
		%ThicknessCalculationErrorGround = true;
		disp(ME.message);
		disp('Stack:');
		for z = 1:length(ME.stack)
			disp([ME.stack(z).file '>' ME.stack(z).name ': line ' num2str(ME.stack(z).line)]);
		end
		%Output = rmfield(Output, 'GND');
		ReturnCode = 1;
		return;
	end
end

if(~ThicknessCalculationError)
	save(ThicknessMatFileName, '-struct', 'Output');
end

PlotGraphics = true;

if(PlotGraphics)
	
	%%
	clf;
	SR = 2;
	SC = 2;
	
	if(~ThicknessCalculationError)
		subplot(SR, SC, 1);
		imshow(Seg.IMG, []);
		hold on;
		%keyboard;
		S = streamline(Output.Streamlines(Output.ValidStreamlines));
		set(S, 'Color', 'r', 'LineWidth', 2);
		Y = plot(Output.Contours.xo, Output.Contours.yo, Output.Contours.xi, Output.Contours.yi);
		set(Y, 'LineWidth', 2);
		title('Streamlines');
		set(gca, 'XTick', [], 'YTick', []);
		
		subplot(SR, SC, 3);
		plot(Output.Thickness);
		xlabel('Node');
		ylabel('Thickness (mm)');
		title('Thickness');
	end
	
	if(GroundTruthGiven)
		if(~ThicknessCalculationErrorGround)
			subplot(SR, SC, 2);
			imshow(Seg.IMG, []);
			hold on;
			S = streamline(Output.Streamlines(Output.ValidStreamlines));
			set(S, 'Color', 'r', 'LineWidth', 2);
			S = streamline(Output.GND.Streamlines(Output.GND.ValidStreamlines));
			set(S, 'Color', 'b', 'LineWidth', 2);
			title('Streamlines + Ground Truth');
			set(gca, 'XTick', [], 'YTick', []);

			subplot(SR, SC, 4);
			Y = plot([Output.Thickness, Output.GND.Thickness(:)]);
			set(Y, 'LineWidth', 2);
			xlabel('Node');
			ylabel('Thickness (mm)');
			title('Thickness (blue) + Ground Truth (green)');
		end
	end
	if(~ThicknessCalculationError)
		OutputPNG = fullfile(pathstr, 'thickness', [FileName '_thickness.png']);
		[~, ~, ~] = mkdir(fullfile(pathstr, 'thickness'));
		FigPos = fullscreen_fig_pos;
		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');

		IMG = imread(OutputPNG);
		IMG = imautocropwhite(IMG);
		imwrite(IMG, OutputPNG);
	end
	%delete(gcf);
end

%save(OutputMAT, '-struct', 'OutputS');
