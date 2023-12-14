function varargout = manual_refine_registration(varargin)
% manual_refine_registration MATLAB code for manual_refine_registration.fig
%      manual_refine_registration, by itself, creates a new manual_refine_registration or raises the existing
%      singleton*.
%
%      H = manual_refine_registration returns the handle to a new manual_refine_registration or the handle to
%      the existing singleton*.
%
%      manual_refine_registration('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in manual_refine_registration.M with the given input arguments.
%
%      manual_refine_registration('Property','Value',...) creates a new manual_refine_registration or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manual_refine_registration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manual_refine_registration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manual_refine_registration

% Last Modified by GUIDE v2.5 19-Dec-2023 15:02:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manual_refine_registration_OpeningFcn, ...
                   'gui_OutputFcn',  @manual_refine_registration_OutputFcn, ...
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


% --- Executes just before manual_refine_registration is made visible.
function manual_refine_registration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manual_refine_registration (see VARARGIN)

% Choose default command line output for manual_refine_registration
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manual_refine_registration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manual_refine_registration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PATHS TO DATA %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function datapath_Callback(hObject, eventdata, handles)
% hObject    handle to datapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datapath as text
%        str2double(get(hObject,'String')) returns contents of datapath as a double


% --- Executes during object creation, after setting all properties.
function datapath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datapath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in datapath_button.
function datapath_button_Callback(hObject, eventdata, handles)
% hObject    handle to datapath_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile('*.tif');
set(handles.datapath,'String',[pathname filename])




% --- Executes on button press in atlaspath_button.
function atlaspath_button_Callback(hObject, eventdata, handles)
% hObject    handle to atlaspath_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = uigetdir;
set(handles.atlaspath,'String',pathname);


function atlaspath_Callback(hObject, eventdata, handles)
% hObject    handle to atlaspath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atlaspath as text
%        str2double(get(hObject,'String')) returns contents of atlaspath as a double


% --- Executes during object creation, after setting all properties.
function atlaspath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlaspath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% INITIALIZE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% turn off some warnings
warning('off','imageio:tiffmexutils:libtiffWarning')
warning('off','imageio:tiffmexutils:libtiffErrorAsWarning')
warning('off','imageio:tifftagsread:expectedTagDataFormatMultiple')
warning('off','MATLAB:polyshape:repairedBySimplify')

% check that file exists
datapath = get(handles.datapath,'String');
while ~exist(datapath,'file')
    set(handles.datapath,'String','file does not exist. reselect.')
end
% once we've verified both files exist
disp(['loading atlas and ' datapath])    
% load atlas & note size
MRIdir = get(handles.atlaspath,'String');
% record some useful information
tmp.atlas = load_tiffs_fast(fullfile(MRIdir,'CCF','average_template_10_coronal.tif'));
tmp.slice_limits = containers.Map({'coronal','axial','sagittal'},...
    {size(tmp.atlas,3),size(tmp.atlas,1),size(tmp.atlas,2)});
tmp.flat_dims = containers.Map({'coronal','axial','sagittal'},{[1 2],[3 2],[1 3]});
set(handles.coronal,'String',['cor (' num2str(tmp.slice_limits('coronal')) ')'])
set(handles.axial,'String',['ax (' num2str(tmp.slice_limits('axial')) ')'])
set(handles.sagittal,'String',['sag (' num2str(tmp.slice_limits('sagittal')) ')'])

% load data  
tmp.ct_path = datapath;    
tifinfo = imfinfo(tmp.ct_path);
% if a 3d tif
if size(tifinfo,1)>1 || (isfield(tifinfo,'ImageDescription') && startsWith(tifinfo.ImageDescription,'ImageJ'))
    tmp.ct = load_tiffs_fast(tmp.ct_path);  
else % if a single slice only
    tmp.ct = imread(tmp.ct_path);
    if size(ct,3)==3            
        tmp.ct = rgb2gray(ct);
    end            
end    
tmp.ct_working = tmp.ct;

% set to middle coronal slice
set(handles.slice_number,'String',num2str(round(tmp.slice_limits('coronal')/2)));

% update GUI: turn off things
set(handles.datapath_button,'Enable','off')
set(handles.datapath,'Enable','off')
set(handles.atlaspath_button,'Enable','off')
set(handles.atlaspath,'Enable','off')
set(handles.ok,'Enable','off')
% update GUI: turn on things
set(handles.coronal,'Enable','on')
set(handles.axial,'Enable','on')
set(handles.sagittal,'Enable','on')
set(handles.slice_number,'Enable','on')
set(handles.prev_slice,'Enable','on')
set(handles.next_slice,'Enable','on')
set(handles.raw,'Enable','on')
set(handles.imadjust,'Enable','on')
set(handles.histeq,'Enable','on')
set(handles.adapthisteq,'Enable','on')
set(handles.show_ct,'Enable','on')
set(handles.show_atlas,'Enable','on')
set(handles.shift_up,'Enable','on')
set(handles.shift_down,'Enable','on')
set(handles.shift_left,'Enable','on')
set(handles.shift_right,'Enable','on')
set(handles.shift_magnitude,'Enable','on')
set(handles.rotate_cw,'Enable','on')
set(handles.rotate_ccw,'Enable','on')
set(handles.rotate_angle,'Enable','on')
set(handles.scale_vert_up,'Enable','on')
set(handles.scale_vert_down,'Enable','on')
set(handles.scale_horiz_up,'Enable','on')
set(handles.scale_horiz_down,'Enable','on')
set(handles.scale_both_up,'Enable','on')
set(handles.scale_both_down,'Enable','on')
set(handles.scale_magnitude,'Enable','on')
set(handles.reset,'Enable','on')
set(handles.save,'Enable','on')
set(handles.quit,'Enable','on')
% let's get started: initialize to middle coronal slice
set(handles.slice_number,'String',num2str(round(tmp.slice_limits('coronal')/2)));
tmp.slice_tforms = [];
view_ct(handles,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% SHIFT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in shift_down.
function shift_down_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
apply_shift(handles,2,tmp.shift);

% --- Executes on button press in shift_up.
function shift_up_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
apply_shift(handles,2,-tmp.shift);

% --- Executes on button press in shift_left.
function shift_left_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
apply_shift(handles,1,-tmp.shift);

% --- Executes on button press in shift_right.
function shift_right_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
apply_shift(handles,1,tmp.shift);

function shift_magnitude_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shift_magnitude as text
%        str2double(get(hObject,'String')) returns contents of shift_magnitude as a double
tmp.shift = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function shift_magnitude_CreateFcn(hObject, eventdata, handles)
global tmp
% hObject    handle to shift_magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
tmp.shift = 2;

function apply_shift(handles,dimension,magnitude) 
global tmp
tform = affine2d;
tform.T(3,dimension) = magnitude;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% ROTATE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in rotate_cw.
function rotate_cw_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to rotate_cw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta = deg2rad(tmp.angle);
tform = affine2d;
% rotation around image center: translate center to origin, rotate, translate back
T1 = tform.T;
T1(3,1:2) = [-floor(size(tmp.current_slice,2)/2) -floor(size(tmp.current_slice,1)/2)];
T2 = tform.T;
T2(1:2,1:2) = [cos(theta) sin(theta); -sin(theta) cos(theta)];
T3 = tform.T;
T3(3,1:2) = -T1(3,1:2);
tform.T = T1*T2*T3;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in rotate_ccw.
function rotate_ccw_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to rotate_ccw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta = -deg2rad(tmp.angle);
tform = affine2d;
% rotation around image center: translate center to origin, rotate, translate back
T1 = tform.T;
T1(3,1:2) = [-floor(size(tmp.current_slice,2)/2) -floor(size(tmp.current_slice,1)/2)];
T2 = tform.T;
T2(1:2,1:2) = [cos(theta) sin(theta); -sin(theta) cos(theta)];
T3 = tform.T;
T3(3,1:2) = -T1(3,1:2);
tform.T = T1*T2*T3;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

function rotate_angle_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to rotate_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotate_angle as text
%        str2double(get(hObject,'String')) returns contents of rotate_angle as a double
tmp.angle = str2double(get(handles.rotate_angle,'String'));

% --- Executes during object creation, after setting all properties.
function rotate_angle_CreateFcn(hObject, eventdata, handles)
global tmp
% hObject    handle to rotate_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
tmp.angle = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% SCALE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in scale_both_up.
function scale_both_up_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_both_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = 1 + tmp.scale;
tform = affine2d;
tform.T(1,1) = s;
tform.T(2,2) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in scale_both_down.
function scale_both_down_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_both_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = 1 - tmp.scale;
tform = affine2d;
tform.T(1,1) = s;
tform.T(2,2) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in scale_horiz_up.
function scale_horiz_up_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_horiz_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('scaling...')
s = 1 + tmp.scale;
tform = affine2d;
tform.T(1,1) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in scale_horiz_down.
function scale_horiz_down_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_horiz_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = 1 - tmp.scale;
tform = affine2d;
tform.T(1,1) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in scale_vert_up.
function scale_vert_up_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_vert_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = 1 + tmp.scale;
tform = affine2d;
tform.T(2,2) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

% --- Executes on button press in scale_vert_down.
function scale_vert_down_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to scale_vert_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = 1 - tmp.scale;
tform = affine2d;
tform.T(2,2) = s;
tmp.slice_tforms = cat(3,tmp.slice_tforms,tform.T); % add this to our chain of transforms
view_ct(handles,0)

function scale_magnitude_Callback(hObject, eventdata, handles)
% hObject    handle to scale_magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale_magnitude as text
%        str2double(get(hObject,'String')) returns contents of scale_magnitude as a double
global tmp
tmp.scale = str2double(get(hObject,'String'))/100;

% --- Executes during object creation, after setting all properties.
function scale_magnitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global tmp
tmp.scale = .01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DISPLAY OPTIONS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global tmp
tmp.contrast = 'raw';

% --- Executes when selected object is changed in contrast.
function contrast_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in contrast 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tmp
tmp.contrast = get(eventdata.NewValue,'Tag');
view_ct(handles,0)

% --- Executes on button press in show_ct.
function show_ct_Callback(hObject, eventdata, handles)
% hObject    handle to show_ct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_ct
view_ct(handles,0)

% --- Executes on button press in show_atlas.
function show_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to show_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_atlas
view_ct(handles,0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ORIENTATION & SLICE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes when selected object is changed in orientation.
function orientation_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in orientation 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tmp
tmp.orientation = get(eventdata.NewValue,'Tag');
% automatically send it to the middle slice
set(handles.slice_number,'String',num2str(round(tmp.slice_limits(tmp.orientation)/2)));
view_ct(handles,1)

% --- Executes during object creation, after setting all properties.
function orientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global tmp
tmp.orientation = 'coronal';


function slice_number_Callback(hObject, eventdata, handles)
% hObject    handle to slice_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slice_number as text
%        str2double(get(hObject,'String')) returns contents of slice_number as a double
global tmp
% turn off next slice if we're at end, and don't let user pick a slice
% bigger than that
if str2double(get(hObject,'String'))>=tmp.slice_limits(tmp.orientation)
    set(hObject,'String',num2str(tmp.slice_limits(tmp.orientation)));
    set(handles.next_slice,'Enable','off');
else
    set(handles.next_slice,'Enable','on');
end 
% turn off prev slice if we're at beginning, and don't let user pick a
% smaller than that
if str2double(get(hObject,'String'))<=1
    set(hObject,'String','1');
    set(handles.prev_slice,'Enable','off');
else
    set(handles.prev_slice,'Enable','on');
end 
view_ct(handles,1)

% --- Executes during object creation, after setting all properties.
function slice_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in prev_slice.
function prev_slice_Callback(hObject, eventdata, handles)
% hObject    handle to prev_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slice_number,'String',num2str(max([str2double(get(handles.slice_number,'String'))-1,1])));
if str2double(get(handles.slice_number,'String'))==1
    set(handles.prev_slice,'Enable','off');
else
    set(handles.prev_slice,'Enable','on');
end 
view_ct(handles,1)

% --- Executes on button press in next_slice.
function next_slice_Callback(hObject, eventdata, handles)
% hObject    handle to next_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tmp
set(handles.slice_number,'String',num2str(min([str2double(get(handles.slice_number,'String'))+1,...
    tmp.slice_limits(tmp.orientation)])));
if str2double(get(handles.slice_number,'String'))==tmp.slice_limits(tmp.orientation)
    set(handles.next_slice,'Enable','off');
else
    set(handles.next_slice,'Enable','on');
end 
view_ct(handles,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% VISUALIZE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function view_ct(handles,apply3d) % apply3d = 0 means just apply to slice, else apply to 3D CT (for changing views, etc)
global tmp

if apply3d == 0 % if we're just adjusting the current slice
    if ~isempty(tmp.slice_tforms)
        apply_tform(2);
    end
    if isfield(tmp,'current_slice_to_view')
        ct_slice_to_view = tmp.current_slice_to_view;
    else
        ct_slice_to_view = tmp.current_slice;    
    end
else % if we're going to a new view
    if ~isempty(tmp.slice_tforms)
        disp('applying transformation to whole CT... ')
        apply_tform(3)
        disp('    ... CT updated')
    end
    switch tmp.orientation
        case 'coronal'
            ct_slice_to_view = tmp.ct_working(:,:,str2double(get(handles.slice_number,'String')));
        case 'axial'
            ct_slice_to_view = transpose(squeeze(tmp.ct_working(str2double(get(handles.slice_number,'String')),:,:)));
        case 'sagittal'
            ct_slice_to_view = squeeze(tmp.ct_working(:,str2double(get(handles.slice_number,'String')),:));
    end
    tmp.current_slice = ct_slice_to_view;
    tmp.current_slice_to_view = ct_slice_to_view;
    tmp.current_slice_orientation = tmp.orientation;
end
% atlas orientation
switch tmp.orientation
    case 'coronal'
         atlas_slice_to_view = tmp.atlas(:,:,str2double(get(handles.slice_number,'String')));
    case 'axial'
        atlas_slice_to_view = transpose(squeeze(tmp.atlas(str2double(get(handles.slice_number,'String')),:,:)));
    case 'sagittal'
        atlas_slice_to_view = squeeze(tmp.atlas(:,str2double(get(handles.slice_number,'String')),:));
end

% contrast
switch tmp.contrast    
    case 'imadjust'
        ct_slice_to_view = imadjust(ct_slice_to_view);
    case 'histeq'
        ct_slice_to_view = histeq(ct_slice_to_view);
    case 'adapthisteq'
        ct_slice_to_view = adapthisteq(ct_slice_to_view);
end

% now show
if get(handles.show_ct,'Value') == 1 && get(handles.show_atlas,'Value')==1
    imshowpair(ct_slice_to_view,atlas_slice_to_view)
elseif get(handles.show_ct,'Value') == 1
    imshow(ct_slice_to_view);
else
    imshow(imadjust(atlas_slice_to_view));
end



function apply_tform(dim) % set dim to 2 to do just the slice, 3 to do the whole CT
global tmp

% combine transformations
sequential_tform = tmp.slice_tforms(:,:,1);
for i = 2:size(tmp.slice_tforms,3)
    sequential_tform = sequential_tform * tmp.slice_tforms(:,:,i);
end
tform = affine2d;
tform.T = sequential_tform;

if dim == 2 % apply transformation to slice only
    tmp.current_slice_to_view = imwarp(tmp.current_slice,tform,'OutputView',imref2d(size(tmp.current_slice),1,1));   
elseif dim == 3 % apply transform to 3D, taking into consideration which view we're in
    switch tmp.current_slice_orientation
        case 'coronal' % x = ML, y = DV        
            tmp.ct_working = imwarp(tmp.ct_working,tform,'OutputView',imref2d(size(tmp.current_slice),1,1));
        case 'axial' % x = ML, y = AP
            tmp.ct_working = permute(imwarp(permute(tmp.ct_working,[3 2 1]),tform,...
                'OutputView',imref2d(size(tmp.current_slice),1,1)),[3 2 1]);

        case 'sagittal' % x = AP, y = DV
            tmp.ct_working = permute(imwarp(permute(tmp.ct_working,[1 3 2]),tform,...
                'OutputView',imref2d(size(tmp.current_slice),1,1)),[1 3 2]);
    end
    tmp.slice_tforms = []; % reset once you've applied it
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% DONE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% reset to original ROIs and radius
tmp.ct_working = tmp.ct;
tmp.slice_tforms = [];
disp('reset')
view_ct(handles,1)

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
global tmp
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('saving tif...')

if contains(path,'Fast_Tiff_Write')
    fast_save_tif([tmp.ct_path(1:end-4) '_manually_adjusted.tif'],tmp.ct_working)
else
    disp('cannot find Fast_Tiff_Write in path. using bfsave, though this is slower.')
    bfsave(tmp.ct_working,[tmp.ct_path(1:end-4) '_manually_adjusted.tif'],...
        'dimensionOrder','XYTZC','BigTiff',true,'Compression', 'LZW');
end
close(gcbf) % close GUI
clear global tmp   % clear our global tmp variable

% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('quitting (no save) ...')
close(gcbf) % close GUI
clear global tmp   % clear our global tmp variable



