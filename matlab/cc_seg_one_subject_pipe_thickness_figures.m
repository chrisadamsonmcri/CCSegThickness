function [ReturnCode] = cc_seg_one_subject_pipe_thickness_figures(MATFilePrefix)

% cc_seg_one_subject_pipe_thickness(MATFile)
%
% DESCRIPTION
%	Performs the thickness computation part of the CC pipeline.
%
% PARAMETERS
%	First argument
%	NIIFileName (char): expected to be a NIFTI file
%	Data (struct) with fields:
%		.AVW: midsagittal slice
%		.PixDims: pixel dimensions (X, Y) of midasgittal slice

%F = who('-file', MATFile);

ReturnCode = 0;

[pathstr, name, ext] = fileparts(MATFilePrefix);

if(strcmp(ext, '.mat'))
	FileName = name;
else
	FileName = [name ext];
end

ManEdit = [];

if(exist(fullfile(pathstr, [FileName '_seg.mat']), 'file') == 2)
	Seg = load(fullfile(pathstr, [FileName '_seg.mat']));
	if(exist(fullfile(pathstr, [FileName '_seg_manedit.mat']), 'file') == 2)
		disp('man edit exists, using');
		ManEdit = load(fullfile(pathstr, [FileName '_seg_manedit.mat']));
	end
else
	disp('Seg file doesnt exist');
	return;
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
	[MATFilePath, name, ext] = fileparts(MATFilePrefix);
	
	[Output.Contours, ...
		Output.SolvedImage, ...
		Output.InitialPoints, ...
		Output.ArcLengthsLeft, ...
		Output.ArcLengthsRight, ...
		Output.LeftOffsetIDXStart, ...
		Output.CurLeftOffsetIDX, ...
		Output.RightOffsetIDXStart, ...
		Output.CurRightOffsetIDX] = endpoint_find_one_subject(FinalSeg, fullfile(MATFilePath, 'paper_figures', [name '-011-ThicknessStartPoints.png']));

	[Output.MaskInnerBoundary, ...
		Output.MaskOuterBoundary, ...
		Output.MaskFree, ...
		Output.SolvedImage, ...
		Output.XY, ...
		Output.X, ...
		Output.Y, ...
		Output.NormFX, ...
		Output.NormFY, ...
		Output.ValidStreamlines, ...
		Output.StartV] = laplace_get_points_2d_auto_mw(Output.Contours.xi, Output.Contours.yi, Output.Contours.xo(end:-1:1), Output.Contours.yo(end:-1:1), 1, 1, 1, 100);

	XYForArcLengths = cell(size(Output.XY));

	for z = 1:length(Output.XY)
		XYForArcLengths{z} = bsxfun(@times, Output.XY{z}, Seg.TemplatePixdims(:)');
	end

	T = cellfun(@arc_length, XYForArcLengths(Output.ValidStreamlines));
	Output.Thickness = zeros(size(Output.XY));
	Output.Thickness(Output.ValidStreamlines) = T;
	clear T;
catch ME
	disp('There was a problem with the thickness calculations');
	ThicknessCalculationError = true;
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
			Output.GND.XY, ...
			Output.GND.X, ...
			Output.GND.Y, ...
			Output.GND.NormFX, ...
			Output.GND.NormFY, ...
			Output.GND.ValidStreamlines, ...
			Output.GND.StartV] = laplace_get_points_2d_auto_mw(Output.GND.Contours.xi, Output.GND.Contours.yi, Output.GND.Contours.xo(end:-1:1), Output.GND.Contours.yo(end:-1:1), 1, 1, 1, 100);
			XYForArcLengths = cell(size(Output.GND.XY));
		%keyboard;
		for z = 1:length(Output.GND.XY)
			XYForArcLengths{z} = bsxfun(@times, Output.GND.XY{z}, Seg.TemplatePixdims(:)');
		end

		T = cellfun(@arc_length, XYForArcLengths(Output.GND.ValidStreamlines));
		Output.GND.Thickness = zeros(size(Output.GND.XY));
		Output.GND.Thickness(Output.GND.ValidStreamlines) = T;
		clear T;
	catch ME
		disp('There was a problem with the thickness calculations, ground truth image');
		ThicknessCalculationErrorGround = true;
		disp(ME.message);
		disp('Stack:');
		for z = 1:length(ME.stack)
			disp([ME.stack(z).file '>' ME.stack(z).name ': line ' num2str(ME.stack(z).line)]);
		end
		Output = rmfield(Output, 'GND');
		ReturnCode = 1;
		return;
	end
end

if(~ThicknessCalculationError)
	save(fullfile(pathstr, [FileName '_thickness.mat']), '-struct', 'Output');
end

PlotGraphics = true;

if(PlotGraphics)
	
	%%
% 	clf;
% 	SR = 2;
% 	SC = 2;
% 	
% 	if(~ThicknessCalculationError)
% 		subplot(SR, SC, 1);
% 		imshow(Seg.IMG, []);
% 		hold on;
% 		%keyboard;
% 		S = streamline(Output.XY(Output.ValidStreamlines));
% 		set(S, 'Color', 'r');
% 		plot(Output.Contours.xo, Output.Contours.yo, Output.Contours.xi, Output.Contours.yi);
% 		title('Streamlines');
% 		
% 		subplot(SR, SC, 3);
% 		plot(Output.Thickness);
% 		xlabel('Node');
% 		ylabel('Thickness (mm)');
% 		title('Thickness');
% 	end
% 	
% 	if(GroundTruthGiven)
% 		if(~ThicknessCalculationErrorGround)
% 			subplot(SR, SC, 2);
% 			imshow(Seg.IMG, []);
% 			hold on;
% 			S = streamline(Output.XY(Output.ValidStreamlines));
% 			set(S, 'Color', 'r');
% 			S = streamline(Output.GND.XY(Output.GND.ValidStreamlines));
% 			set(S, 'Color', 'b');
% 			title('Streamlines + Ground Truth');
% 
% 			subplot(SR, SC, 4);
% 			plot([Output.Thickness, Output.GND.Thickness(:)]);
% 			xlabel('Node');
% 			ylabel('Thickness (mm)');
% 			title('Thickness (blue) + Ground Truth (green)');
% 		end
% 	end
% 	if(~ThicknessCalculationError)
% 		OutputPNG = fullfile(pathstr, 'paper_figures', [FileName '-012-Thickness.png']);
% 		[~, ~, ~] = mkdir(fullfile(pathstr, 'paper_figures'));
% % 		FigPos = fullscreen_fig_pos;
% % 		set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% % 		exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% % 
% % 		IMG = imread(OutputPNG);
% % 		IMG = imautocropwhite(IMG);
% % 		imwrite(IMG, OutputPNG);
% 		cc_seg_save_figure_paper(OutputPNG);
% 	end
	if(GroundTruthGiven)
		if(~ThicknessCalculationErrorGround)
			%%
			clf;
			SR = 2;
			SC = 2;
			TitleProps = {'FontSize', 20};
			TextProps = {'FontSize', 16};
			AX = zeros(1, 3);
			AX(1) = subplot(SR, SC, 1);
			imshow(Seg.IMG, []);
			%keyboard;
			hold on;
			SEST = plot([Output.Contours.xi; Output.Contours.xo(end:-1:1); Output.Contours.xi(1)], [Output.Contours.yi; Output.Contours.yo(end:-1:1); Output.Contours.yi(1)], 'Color', 'b', 'LineWidth', 3);
			SGND = plot([Output.GND.Contours.xi; Output.GND.Contours.xo(end:-1:1); Output.GND.Contours.xi(1)], [Output.GND.Contours.yi; Output.GND.Contours.yo(end:-1:1); Output.GND.Contours.yi(1)], 'Color', [0 0.5 0], 'LineWidth', 3);
			
			legend([SEST(1), SGND(1)]', 'Estimated', 'Ground Truth', 'Location', 'South');
			set(gca, TextProps{:}, 'XTick', [], 'YTick', []);
			title([roman_label(1) ' Segmentations'], TitleProps{:});
			
			AX(2) = subplot(SR, SC, 2);
			imshow(Seg.IMG, []);
			hold on;
			SEST = streamline(Output.XY(Output.ValidStreamlines));
			set(SEST, 'Color', 'b', 'LineWidth', 3);
			SGND = streamline(Output.GND.XY(Output.GND.ValidStreamlines));
			set(SGND, 'Color', [0 0.5 0], 'LineWidth', 3);
			
			IDX = [1, 20:20:100];
			T = Output.GND.XY;
			GNDDIFF = cellfun(@diff, T, 'UniformOutput', false);
			
			FirstPoints = zeros(length(IDX), 2);
			Vectors = zeros(length(IDX), 2);
			for z = 1:length(IDX)
				FirstPoints(z, :) = T{IDX(z)}(1, :);
				Vectors(z, :) = -GNDDIFF{IDX(z)}(1, :);
			end
			IDXString = num2cell(IDX);
			IDXString = cellfun(@num2str, IDXString, 'UniformOutput', false);
			Vectors = bsxfun(@rdivide, Vectors, sqrt(sum(Vectors .* Vectors, 2)));
			Vectors = Vectors * 5;
			text(FirstPoints(:, 1) + Vectors(:, 1), FirstPoints(:, 2) + Vectors(:, 2), IDXString, ...
				'HorizontalAlignment', 'center', ...
				'VerticalAlignment', 'middle', ...
				'BackgroundColor', 'w');
			title([roman_label(2) ' Streamlines'], TitleProps{:});
			legend([SEST(1), SGND(1)]', 'Estimated', 'Ground Truth', 'Location', 'South');
			fix_figure_aspect(gcf, gca, size(Seg.IMG));
			set(gca, TextProps{:}, 'XTick', [], 'YTick', []);
			
			AX(3) = subplot(SR, SC, [3 4]);
			plot([Output.Thickness, Output.GND.Thickness(:)], 'LineWidth', 3);
			xlabel('Node', TextProps{:});
			ylabel('Thickness (mm)', TextProps{:});
			title([roman_label(3) ' Thickness Profiles'], TitleProps{:});
			legend({'Estimated', 'Ground Truth'}, 'Location', 'North');
			AXPosLeft = get(AX(1), 'Position');
			AXPosRight = get(AX(2), 'Position');
			AXPosRight([2 4]) = AXPosLeft([2 4]);
			AXPosRight(1) = AXPosRight(1) - 0.04;
			set(AX(2), 'Position', AXPosRight);
			set(gca, TextProps{:}, 'XTick', IDX);
			%%
			%keyboard;
		end
		if(~ThicknessCalculationError)
			OutputPNG = fullfile(pathstr, 'paper_figures', [FileName '-013-ThicknessGroundTruth.png']);
% 			[~, ~, ~] = mkdir(fullfile(pathstr, 'paper_figures'));
% 			FigPos = fullscreen_fig_pos;
% 			set(gcf, 'Position', FigPos, 'PaperPosition', FigPos, 'PaperUnits', 'points');
% 			exportfig(gcf, OutputPNG, 'Format', 'png', 'Width', FigPos(3), 'Height', FigPos(4), 'Color', 'rgb');
% 
% 			IMG = imread(OutputPNG);
% 			IMG = imautocropwhite(IMG);
% 			imwrite(IMG, OutputPNG);
			cc_seg_save_figure_paper_png(OutputPNG);
		end
	end
	%delete(gcf);
end

%save(OutputMAT, '-struct', 'OutputS');
