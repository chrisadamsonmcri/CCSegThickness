function varargout = cc_seg_manedit(varargin)
% CC_SEG_MANEDIT MATLAB code for cc_seg_manedit.fig
%      CC_SEG_MANEDIT, by itself, creates a new CC_SEG_MANEDIT or raises the existing
%      singleton*.
%
%      H = CC_SEG_MANEDIT returns the handle to a new CC_SEG_MANEDIT or the handle to
%      the existing singleton*.
%
%      CC_SEG_MANEDIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CC_SEG_MANEDIT.M with the given input arguments.
%
%      CC_SEG_MANEDIT('Property','Value',...) creates a new CC_SEG_MANEDIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cc_seg_manedit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cc_seg_manedit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Copyright 2013, Chris Adamson, Murdoch Childrens Research Institute
% See LICENSE for full license information.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cc_seg_manedit

% Last Modified by GUIDE v2.5 10-Oct-2013 13:36:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cc_seg_manedit_OpeningFcn, ...
                   'gui_OutputFcn',  @cc_seg_manedit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cc_seg_manedit is made visible.
function cc_seg_manedit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cc_seg_manedit (see VARARGIN)

% Choose default command line output for cc_seg_manedit
handles.output = hObject;
% 
% AllObjects = [handles.figure1,  handles.sld_opacity, ...
%          handles.text_opacity, handles.text_strelsize, ...
%         handles.sld_strelsize, handles.text3,...
%             handles.btn_close, handles.btn_open,...
%             handles.btn_erode, handles.btn_dilate,...
%       handles.text_cursubject, handles.btn_save,...
%        handles.sld_cursubject, handles.sld_brushsize,...
%        handles.text_brushsize, handles.radio_brush_circle,...
%    handles.radio_brush_square];
% 
% set(AllObjects, 'KeyPressFcn', @KeyPressFcn);

if(length(varargin) < 1)
	error('Give a directory');
end
handles.Directory = varargin{1};

if(exist(handles.Directory, 'dir') ~= 7)
	error('Argument given is not a directory');
end

D = dir(fullfile(handles.Directory, '*_seg.mat'));
[F{1:length(D)}] = deal(D.name);

[match] = regexp(F, '.+_manedit\.mat$', 'match');
NotManEdit = cellfun(@isempty, match);
clear match;
F = F(NotManEdit);

[match] = regexp(F, '.+_thickness\.mat$', 'match');
NotManEdit = cellfun(@isempty, match);
clear match;
F = F(NotManEdit);

handles.SegMATs = cell(1, length(F));
handles.SegManEdits = cell(1, length(F));
handles.ImageVariableName = 'IMG';
handles.InitialSegVariableName = 'FinalSeg';
handles.ManEditVariableName = 'FinalSegManEdit';

for z = 1:length(F)
	%disp(fullfile(handles.Directory, F{z}));
	handles.SegMATs{z} = load(fullfile(handles.Directory, F{z}), handles.ImageVariableName, handles.InitialSegVariableName);
	[~, name] = fileparts(F{z});
	%clear ext;
	if(exist(fullfile(handles.Directory, [name '_manedit.mat']), 'file') == 2)
		T = load(fullfile(handles.Directory, [name '_manedit.mat']));
		if(isfield(T, handles.ManEditVariableName))
			handles.SegManEdits{z} = T.(handles.ManEditVariableName);
			CC = bwconncomp(handles.SegManEdits{z});
			%keyboard;
			if(CC.NumObjects > 1)
				disp([name ', ' num2str(z) ': multiple regions']);
			end
		else
			disp(['Warning: the manual edit file does not have the right variable: ' F{z}]);
		end
	end
end

set(handles.sld_brushsize, 'Min', 1, 'Max', 10, 'Value', 1);
set(handles.sld_brushsize, 'SliderStep', [1 ./ (get(handles.sld_brushsize, 'Max') - get(handles.sld_brushsize, 'Min')), 1 ./ (get(handles.sld_brushsize, 'Max') - get(handles.sld_brushsize, 'Min'))]);
set(handles.text_brushsize, 'String', ['Brush size: ' num2str(get(handles.sld_brushsize, 'Value'))]);
%F
set(handles.sld_cursubject, 'Min', 1, 'Max', length(F), 'Value', 1);
set(handles.sld_cursubject, 'SliderStep', [1 ./ (get(handles.sld_cursubject, 'Max') - get(handles.sld_cursubject, 'Min')), 1 ./ (get(handles.sld_cursubject, 'Max') - get(handles.sld_cursubject, 'Min'))]);

set(handles.sld_strelsize, 'Min', 1, 'Max', 10, 'Value', 2);
set(handles.sld_strelsize, 'SliderStep', [1 ./ (get(handles.sld_strelsize, 'Max') - get(handles.sld_strelsize, 'Min')), 1 ./ (get(handles.sld_strelsize, 'Max') - get(handles.sld_strelsize, 'Min'))]);
set(handles.text_strelsize, 'String', ['Strel size: ' num2str(get(handles.sld_strelsize, 'Value'))]);

set(handles.sld_opacity, 'Min', 0, 'Max', 1, 'Value', 0.7);
set(handles.sld_opacity, 'SliderStep', [0.1 0.3]);
set(handles.text_opacity, 'String', ['Opacity: ' num2str(get(handles.sld_opacity, 'Value'), '%.1f')]);
handles.OriginalAxisPosition = get(handles.axes_main, 'Position');

axes(handles.axes_main);

handles.NeedToSave = false(1, length(F));
handles.MouseDown = false;
handles.FileNames = F;

handles.Undo.IMG = [];
handles.Undo.Subject = [];

handles = subject_changed(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cc_seg_manedit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [handles] = copy_image_to_edited(handles)

CurSubject = get(handles.sld_cursubject, 'Value');

if(isempty(handles.SegManEdits{CurSubject}))
	handles.SegManEdits{CurSubject} = handles.SegMATs{CurSubject}.(handles.InitialSegVariableName);
end
handles.NeedToSave(CurSubject) = true;

function [handles] = subject_changed(handles)

CurSubject = get(handles.sld_cursubject, 'Value');
set(handles.axes_main, 'Position', handles.OriginalAxisPosition);
hold off;
IMG = handles.SegMATs{CurSubject}.(handles.ImageVariableName);
T = sort(IMG(:));

%imagesc(IMG, [T(floor(length(T) * 0.05)), T(ceil(length(T) * 0.95))]);
imagesc(IMG);
axis equal;
hold on;
if(isempty(handles.SegManEdits{CurSubject}))
	CurSeg = handles.SegMATs{CurSubject}.(handles.InitialSegVariableName);
else
	CurSeg = handles.SegManEdits{CurSubject};
end
RGB = cat(3, double(CurSeg), zeros([size(CurSeg), 2]));
handles.SegHandle = imagesc(RGB, 'AlphaData', double(CurSeg) * get(handles.sld_opacity, 'Value'));
axis off;
colormap gray;
%fix_figure_aspect(handles.figure1, handles.axes_main, size(IMG));

set(handles.text_cursubject, 'String', ['Subject ' num2str(CurSubject) ', ' handles.FileNames{CurSubject}]);
%guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = cc_seg_manedit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sld_brushsize_Callback(hObject, eventdata, handles)
% hObject    handle to sld_brushsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject, 'Value', round(get(hObject, 'Value')));
set(handles.text_brushsize, 'String', ['Brush size: ' num2str(get(hObject, 'Value'))]);

% --- Executes during object creation, after setting all properties.
function sld_brushsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sld_brushsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sld_cursubject_Callback(hObject, eventdata, handles)
% hObject    handle to sld_cursubject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject, 'Value', round(get(hObject, 'Value')));
handles = subject_changed(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sld_cursubject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sld_cursubject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btn_save.
function btn_save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(any(handles.NeedToSave))
	I = find(handles.NeedToSave);
	for z = 1:length(I)
		[pathstr, name, ext] = fileparts(handles.FileNames{I(z)});
		T.(handles.ManEditVariableName) = handles.SegManEdits{I(z)};
		save(fullfile(handles.Directory, [name '_manedit.mat']), '-struct', 'T');
		clear T;
	end
end
handles.NeedToSave = false(size(handles.FileNames));
set(handles.btn_save, 'Enable', 'off');
guidata(hObject, handles);


% --- Executes on button press in btn_undo.
function btn_undo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [STR] = circular_mask(R)

xx = -R:R;
yy = -R:R;

[X, Y] = meshgrid(xx, yy);

MAG = sqrt(X .* X + Y .* Y);
STR = (MAG <= R);

function [M] = mask_circular_dilate(I, R)

M = conv2(single(I), single(circular_mask(R)), 'same');
M = (M > 0);

function [M] = mask_circular_erode(I, R)

M = ~mask_circular_dilate(~I, R);

function [M] = mask_circular_open(I, R)

M = mask_circular_erode(I, R);
M = mask_circular_dilate(M, R);

function [M] = mask_circular_close(I, R)

M = mask_circular_dilate(I, R);
M = mask_circular_erode(M, R);


% --- Executes on button press in btn_dilate.
function btn_dilate_Callback(hObject, eventdata, handles)
% hObject    handle to btn_dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);
handles.Undo.IMG = handles.SegManEdits{CurSubject};
handles.Undo.Subject = CurSubject;

%handles.SegManEdits{CurSubject} = imdilate(handles.SegManEdits{CurSubject}, strel('disk', get(handles.sld_strelsize, 'Value')));
handles.SegManEdits{CurSubject} = mask_circular_dilate(handles.SegManEdits{CurSubject}, get(handles.sld_strelsize, 'Value'));
handles = subject_changed(handles);
handles.NeedToSave(CurSubject) = 1;
set(handles.btn_save, 'Enable', 'on');

guidata(hObject, handles);



% --- Executes on button press in btn_erode.
function btn_erode_Callback(hObject, eventdata, handles)
% hObject    handle to btn_erode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);

handles.Undo.IMG = handles.SegManEdits{CurSubject};
handles.Undo.Subject = CurSubject;

%handles.SegManEdits{CurSubject} = imerode(handles.SegManEdits{CurSubject}, strel('disk', get(handles.sld_strelsize, 'Value')));
handles.SegManEdits{CurSubject} = mask_circular_erode(handles.SegManEdits{CurSubject}, get(handles.sld_strelsize, 'Value'));
handles = subject_changed(handles);
handles.NeedToSave(CurSubject) = 1;
set(handles.btn_save, 'Enable', 'on');

guidata(hObject, handles);


% --- Executes on button press in btn_open.
function btn_open_Callback(hObject, eventdata, handles)
% hObject    handle to btn_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);
handles.Undo.IMG = handles.SegManEdits{CurSubject};
handles.Undo.Subject = CurSubject;

%handles.SegManEdits{CurSubject} = imopen(handles.SegManEdits{CurSubject}, strel('disk', get(handles.sld_strelsize, 'Value')));
handles.SegManEdits{CurSubject} = mask_circular_open(handles.SegManEdits{CurSubject}, get(handles.sld_strelsize, 'Value'));
handles = subject_changed(handles);
handles.NeedToSave(CurSubject) = 1;
set(handles.btn_save, 'Enable', 'on');

guidata(hObject, handles);

% --- Executes on button press in btn_close.
function btn_close_Callback(hObject, eventdata, handles)
% hObject    handle to btn_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);
handles.Undo.IMG = handles.SegManEdits{CurSubject};
handles.Undo.Subject = CurSubject;

%handles.SegManEdits{CurSubject} = imclose(handles.SegManEdits{CurSubject}, strel('disk', get(handles.sld_strelsize, 'Value')));
handles.SegManEdits{CurSubject} = mask_circular_close(handles.SegManEdits{CurSubject}, get(handles.sld_strelsize, 'Value'));
handles = subject_changed(handles);
handles.NeedToSave(CurSubject) = 1;
set(handles.btn_save, 'Enable', 'on');
guidata(hObject, handles);


% --- Executes on slider movement.
function sld_strelsize_Callback(hObject, eventdata, handles)
% hObject    handle to sld_strelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject, 'Value', round(get(hObject, 'Value')));
set(handles.text_strelsize, 'String', ['Strel size: ' num2str(get(hObject, 'Value'))]);

% --- Executes during object creation, after setting all properties.
function sld_strelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sld_strelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function [P] = get_axis_point_from_figure_pos(handles, CurrentPoint)

CurSubject = get(handles.sld_cursubject, 'Value');

SZ = size(handles.SegMATs{CurSubject}.(handles.ImageVariableName));

AxesPos = get(handles.axes_main, 'Position');

XL = get(handles.axes_main, 'XLim');
YL = get(handles.axes_main, 'YLim');

AxesXFrac = (CurrentPoint(1) - AxesPos(1)) ./ AxesPos(3);
AxesYFrac = 1 - (CurrentPoint(2) - AxesPos(2)) ./ AxesPos(4);

P(1) = XL(1) + AxesXFrac * (XL(2) - XL(1));
P(2) = YL(1) + AxesYFrac * (YL(2) - YL(1));
%P(1) = YL(1) + DataYFrac * (YL(2) - YL(1))
%P
%P(1) = AxesXFrac * SZ(2);

P = round(P);
%SZ
P(P < 1 | P > SZ([2 1])) = NaN;


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp('figure1_WindowButtonMotionFcn');
if(isfield(handles, 'MouseDown'))
	if(handles.MouseDown)
        G = get(hObject, 'Units');
        set(hObject, 'Units', 'normalized');
		
        F = get(hObject, 'CurrentPoint');
        set(hObject, 'Units', G);
		T = get_axis_point_from_figure_pos(handles, F);

		if(all(~isnan(T)))
			handles.CurPointsClicked = unique(cat(1, handles.CurPointsClicked, T), 'rows');
		end
	end
end
guidata(hObject, handles);

%get(hObject, 'CurrentPoint')
%guidata(hObject, handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.MouseDown = true;
handles.SelectionType = get(hObject, 'SelectionType');
%disp('figure1_WindowButtonDownFcn');
G = get(hObject, 'Units');
set(hObject, 'Units', 'normalized');
%get(hObject, 'CurrentPoint')
T = get_axis_point_from_figure_pos(handles, get(hObject, 'CurrentPoint'));
set(hObject, 'Units', G);
if(all(~isnan(T)))
	handles.CurPointsClicked = T;
else
	handles.CurPointsClicked = [];
end
%handles.CurPointsClicked
guidata(hObject, handles);

%function [I] = get_indices_to_edit(handles)

function [I] = distance_function_square(P, BrushSize, IMGSZ)

%[X, Y] = meshgrid(1:IMGSZ(2), 1:IMGSZ(1));

I = zeros(IMGSZ, 'single');

ID = sub2ind(IMGSZ, P(:, 2), P(:, 1));
I(ID) = 1;

I = conv2(I, ones(1, BrushSize, 'single'), 'same');
I = conv2(I, ones(BrushSize, 1, 'single'), 'same');
I = (I > 0);
%for z = 1:size(P, 1)
%	I = I | ((abs(X - P(z, 1)) < BrushSize / 2) & (abs(Y - P(z, 2)) < BrushSize / 2));
%end

function [I] = distance_function_circle(P, BrushSize, IMGSZ)

[X, Y] = meshgrid(1:IMGSZ(2), 1:IMGSZ(1));

I = false(IMGSZ);

for z = 1:size(P, 1)
	XC = X - P(z, 1);
	YC = Y - P(z, 2);
	I = I | (sqrt(XC .* XC + YC .* YC) < BrushSize);
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
if(get(handles.radio_brush_circle, 'Value') == get(handles.radio_brush_circle, 'Max'))
	DistFunc = @distance_function_circle;
else
	DistFunc = @distance_function_square;
end
%DistFunc
CurBrushSize = get(handles.sld_brushsize, 'Value');

handles.MouseDown = false;
if(~isempty(handles.CurPointsClicked))
	I = feval(DistFunc, handles.CurPointsClicked, CurBrushSize, size(handles.SegMATs{CurSubject}.(handles.ImageVariableName)));
	handles = copy_image_to_edited(handles);
	handles.Undo.IMG = handles.SegManEdits{CurSubject};
	handles.Undo.Subject = CurSubject;

	switch(handles.SelectionType)
		case 'normal' % left mouse button
			handles.SegManEdits{CurSubject} = handles.SegManEdits{CurSubject} | I;
		case 'alt' % right mouse button
			handles.SegManEdits{CurSubject} = handles.SegManEdits{CurSubject} & ~I;
	end
	handles = subject_changed(handles);
	handles.NeedToSave(CurSubject) = 1;
	set(handles.btn_save, 'Enable', 'on');
end

guidata(hObject, handles);


% --- Executes on slider movement.
function sld_opacity_Callback(hObject, eventdata, handles)
% hObject    handle to sld_opacity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.text_opacity, 'String', ['Opacity: ' num2str(get(handles.sld_opacity, 'Value'), '%.1f')]);
handles = subject_changed(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sld_opacity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sld_opacity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on btn_dilate and none of its controls.
function btn_dilate_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to btn_dilate (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(eventdata.Key, 'z') && strcmp(eventdata.Modifier, 'control'))
	if(~isempty(handles.Undo.IMG))
		handles.SegManEdits{handles.Undo.Subject} = handles.Undo.IMG;
		handles.Undo.IMG = [];
		handles.Undo.Subject = [];
	end
	handles = subject_changed(handles);
end
guidata(hObject, handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in btn_largest_component.
function btn_largest_component_Callback(hObject, eventdata, handles)
% hObject    handle to btn_largest_component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);
if(any(handles.SegManEdits{CurSubject}(:)))
	handles.Undo.IMG = handles.SegManEdits{CurSubject};
	handles.Undo.Subject = CurSubject;

	SegCC = bwconncomp(handles.SegManEdits{CurSubject});
	Areas = cellfun('length', SegCC.PixelIdxList);
	[~, I] = max(Areas);
	handles.SegManEdits{CurSubject}(:) = false;
	handles.SegManEdits{CurSubject}(SegCC.PixelIdxList{I}) = 1;
	handles = subject_changed(handles);
	handles.NeedToSave(CurSubject) = 1;
	set(handles.btn_save, 'Enable', 'on');
end
guidata(hObject, handles);


% --- Executes on button press in btn_fill_holes.
function btn_fill_holes_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fill_holes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurSubject = get(handles.sld_cursubject, 'Value');
handles = copy_image_to_edited(handles);
if(any(handles.SegManEdits{CurSubject}(:)))
	handles.Undo.IMG = handles.SegManEdits{CurSubject};
	handles.Undo.Subject = CurSubject;

	handles.SegManEdits{CurSubject}(:) = bwfill(handles.SegManEdits{CurSubject}, 'holes');
	handles = subject_changed(handles);
	handles.NeedToSave(CurSubject) = 1;
	set(handles.btn_save, 'Enable', 'on');
end
guidata(hObject, handles);
