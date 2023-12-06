function varargout = orient_CT(varargin)
% ORIENT_CT MATLAB code for orient_CT.fig
%      ORIENT_CT, by itself, creates a new ORIENT_CT or raises the existing
%      singleton*.
%
%      H = ORIENT_CT returns the handle to a new ORIENT_CT or the handle to
%      the existing singleton*.
%
%      ORIENT_CT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ORIENT_CT.M with the given input arguments.
%
%      ORIENT_CT('Property','Value',...) creates a new ORIENT_CT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before orient_CT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to orient_CT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help orient_CT

% Last Modified by GUIDE v2.5 05-Dec-2023 10:14:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @orient_CT_OpeningFcn, ...
                   'gui_OutputFcn',  @orient_CT_OutputFcn, ...
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


% --- Executes just before orient_CT is made visible.
function orient_CT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to orient_CT (see VARARGIN)

% Choose default command line output for orient_CT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes orient_CT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = orient_CT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%
%%%% PATH SETUP & %%%%
%%%% INSTRUCTIONS %%%%
%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in dicom_dir_button.
function dicom_dir_button_Callback(hObject, eventdata, handles)
% hObject    handle to dicom_dir_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dicom_dir=uigetdir;
set(handles.dicom_dir_text,'String',dicom_dir);
set(handles.output_filepath_text,'String',[dicom_dir '.tif']);
set(handles.ok_go_button,'Enable','on')


function dicom_dir_text_Callback(hObject, eventdata, handles)
% hObject    handle to dicom_dir_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dicom_dir_text as text
%        str2double(get(hObject,'String')) returns contents of dicom_dir_text as a double


% --- Executes during object creation, after setting all properties.
function dicom_dir_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dicom_dir_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output_filepath_text_Callback(hObject, eventdata, handles)
% hObject    handle to output_filepath_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_filepath_text as text
%        str2double(get(hObject,'String')) returns contents of output_filepath_text as a double


% --- Executes during object creation, after setting all properties.
function output_filepath_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_filepath_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function file_prefix_text_Callback(hObject, eventdata, handles)
% hObject    handle to file_prefix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_prefix_text as text
%        str2double(get(hObject,'String')) returns contents of file_prefix_text as a double


% --- Executes during object creation, after setting all properties.
function file_prefix_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_prefix_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_prefix_text as text
%        str2double(get(hObject,'String')) returns contents of file_prefix_text as a double

%%%%%%%%%%%%%%%%%%%%%%
%%%% BUTTONS AND %%%%%
%%%%%%% STUFF %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in ok_go_button.
function ok_go_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_go_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global temp

% paths & filenames
temp.dicom_dir = get(handles.dicom_dir_text,'String');
temp.output_filepath = get(handles.output_filepath_text,'String');
temp.dicom_prefix = get(handles.file_prefix_text,'String');
temp.im_guide = imread('CT_dicom_align_and_convert_guide.png');
% load and crop dicoms
load_and_crop;

% enable some things
set(handles.reset_button,'Enable','on')
set(handles.next_button,'Enable','on')
set(handles.imadjust_radio,'Enable','on')
set(handles.histeq_radio,'Enable','on')
set(handles.adapthisteq_radio,'Enable','on')
set(handles.gray_radio,'Enable','on')
set(handles.parula_radio,'Enable','on')
set(handles.autumn_radio,'Enable','on')

% preallocate some things
temp.temp_points = [];
temp.midsag_points = [];
temp.current_slice = 50;
temp.current_step = 1;

% now select midsagittal points
select_points(handles, 5);

% --- Executes on button press in next_button.
function next_button_Callback(hObject, eventdata, handles)
% hObject    handle to next_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global temp
% if we're on step 1 (midsag alignment)
if temp.current_step == 1
    % add the points onto our list, reset temp_points
    if ~isempty(temp.temp_points)
        temp.midsag_points = cat(1,temp.midsag_points,[temp.temp_points ones(size(temp.temp_points,1),1)*temp.current_slice]);
        temp.temp_points = [];
    end
    % increment the slice number if we're not at the end
    if temp.current_slice < size(temp.data,3)-20
        temp.current_slice = temp.current_slice+20;
        % set this label to "Done" if we're on the last slice
        if temp.current_slice > size(temp.data,3)-20
            set(handles.next_button,'String','Done Step 1')
        end
        select_points(handles,5);
    else
        align_midsag_points;
        % update stage and some GUI things
        temp.current_step = 2;
        instr_str = {'DETERMINE VOLUME ML ORIENTATION (Step 2/4)';'Please click the right hemisphere (the implant should help you figure this out).'};
        set(handles.instructions,'String',instr_str);
        set(handles.next_button,'String','Done Step 2');
        temp.temp_points = [];
        temp.current_slice = round(size(temp.data,3)*3/4);
        select_points(handles,1);
    end 
% if we're doing ML orientation alignment
elseif temp.current_step == 2 
    % do step 2
    orient_volume_ML;
    % move on to step 3
    temp.current_step = 3;
    instr_str = {'DETERMINE VOLUME AP ORIENTATION (Step 3/4)';'Please click the front of the brain.'};
    set(handles.instructions,'String',instr_str);
    set(handles.next_button,'String','Done Step 3');
    temp.temp_points = [];
    temp.current_slice = round(size(temp.data,3)/2);
    select_points(handles,1);
% if we're doing AP orientatation alignment
elseif temp.current_step == 3
    % do step 3 and then move on
    orient_volume_AP;
    % move on to step 4
    temp.current_step = 4;
    instr_str = {'ADJUST VOLUME PITCH ORIENTATION (Step 4/4)';'On the midsagittal slice shown, please click the point where the cerebral aqueduct meets the 4th ventricle, and the anterior ventral point of the Islands of Calleja (see pop-up guide). Guess if necessary.'};
    set(handles.instructions,'String',instr_str);
    set(handles.next_button,'String','Done Step 4');
    temp.gui_fig = gcf;
    temp.guide_fig = figure;
    imshow(temp.im_guide)
    temp.temp_points = [];
    temp.current_slice = round(size(temp.data,2)/2);
    select_points2(handles,2);
% if we're doing the AP orientation alignment
elseif temp.current_step == 4
    % do step 4
    orient_volume_pitch;
    % save out
    close(temp.guide_fig);    
    %disp('saving 2 tif files:')
    %disp(' 1. _coronal.tif the CT in coronal orientation')
    %disp(' 2. _coronal_histeq.tif: the CT in coronal orientation, histeq for visualization')    
    %save_tiff(temp.data,temp.output_filepath);
    save_tiff(permute(temp.data,[3 2 1]),[temp.output_filepath(1:end-4) '_coronal.tif'],1);
    %save_tiff(histeq(permute(temp.data,[3 2 1])),[temp.output_filepath(1:end-4) '_coronal_histeq.tif'],0);    
end
    
    

% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global temp
if temp.current_step == 1
    temp.temp_points = [];
    if temp.current_step == 1
        select_points(handles,5)
    elseif temp.current_step == 2 
        select_points(handles,1)
    elseif temp.current_step == 3
        select_points(handles,1)
    elseif temp.current_step == 4
        select_points2(handles,2)
    end
else
end


% --- Executes on button press in imadjust_radio.
function imadjust_radio_Callback(hObject, eventdata, handles)
% hObject    handle to imadjust_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of imadjust_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end

% --- Executes on button press in histeq_radio.
function histeq_radio_Callback(hObject, eventdata, handles)
% hObject    handle to histeq_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of histeq_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end

% --- Executes on button press in adapthisteq_radio.
function adapthisteq_radio_Callback(hObject, eventdata, handles)
% hObject    handle to adapthisteq_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adapthisteq_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end


% --- Executes on button press in gray_radio.
function gray_radio_Callback(hObject, eventdata, handles)
% hObject    handle to gray_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gray_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end

% --- Executes on button press in parula_radio.
function parula_radio_Callback(hObject, eventdata, handles)
% hObject    handle to parula_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parula_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end

% --- Executes on button press in autumn_radio.
function autumn_radio_Callback(hObject, eventdata, handles)
% hObject    handle to autumn_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autumn_radio
global temp
if temp.current_step == 1 && get(hObject,'Value')==1
    select_points(handles,5)
elseif temp.current_step == 2 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 3 && get(hObject,'Value')==1
    select_points(handles,1)
elseif temp.current_step == 4 && get(hObject,'Value')==1
    select_points2(handles,2)
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%% INTERNAL %%%%%%%
%%%%% FUNCTIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

function load_and_crop
global temp
% let's figure out all the files in there
disp('loading dicoms...')
dicom_files = dir(temp.dicom_dir);
dicom_files = {dicom_files.name}';
dicom_files = dicom_files(startsWith(dicom_files,temp.dicom_prefix));

% now let's read it in to a 3D matrix
% read in the first one
data = dicomread(fullfile(temp.dicom_dir,dicom_files{1}));
% preallocate the rest
data = cat(3,data,nan(size(data,1),size(data,2),numel(dicom_files)-1));
% now read the rest in
for i = 2:numel(dicom_files)
    tmp = dicomread(fullfile(temp.dicom_dir,dicom_files{i}));
    data(:,:,i) = tmp;
end

% now we've gotta figure out which z-slices have the data so we can crop
z_means = squeeze(mean(double(data),[1:2]));
z_means_max = max(z_means);
z_means_min = min(z_means);
% figure out if it goes ventral-to-dorsal or dorsal-to-ventral
if z_means(1)<z_means(end)
    z_start = find(z_means>z_means_min+.05*(z_means_max-z_means_min),1,'first');
    z_end = numel(z_means);
    temp.vd = 1;
else
    z_start = 1;
    z_end = find(z_means>z_means_min+.05*(z_means_max-z_means_min),1,'last');
    temp.vd = 0;
end
data = data(:,:,z_start:z_end);

% let's also trim cols and rows with no data
colsum = sum(data,[1 3]);
rowsum = sum(data,[2 3]);
temp.data = data(rowsum>0,colsum>0,:);

disp('dicoms loaded.')


function select_points(handles,npts)
global temp

% display the image
cla; % clear axis first just in case
hold on; % hol don
% color options
if get(handles.gray_radio,'Value')==1
    colormap gray
elseif get(handles.parula_radio,'Value')==1
    colormap parula
else
    colormap autumn  
end
% display option
if get(handles.imadjust_radio,'Value')==1
    imagesc(handles.image_axes,imadjust(temp.data(:,:,temp.current_slice)))
elseif get(handles.histeq_radio,'Value')==1
    imagesc(handles.image_axes,histeq(temp.data(:,:,temp.current_slice)))
else
    imagesc(handles.image_axes,adapthisteq(temp.data(:,:,temp.current_slice)))    
end

% plot points if there are any
if ~isempty(temp.temp_points)
    plot(handles.image_axes,temp.temp_points(:,1),temp.temp_points(:,2),'o','Color','k','MarkerSize',6,'MarkerFaceColor',lines(1));
end
set(gca,'YDir','reverse','YLim',[.5 size(temp.data,1)+.5],'XLim',[.5 size(temp.data,2)+.5])

% now get points
try
    while size(temp.temp_points,1)<npts
        upoint=drawpoint(handles.image_axes,'InteractionsAllowed','none');
        x = upoint.Position(1);
        y = upoint.Position(2);
        temp.temp_points = cat(1,temp.temp_points,[x y]);
    end
catch exception
end

function select_points2(handles,npts)
global temp

% display the image
figure(temp.gui_fig);
cla; % clear axis first just in case
hold on; 
% color options
if get(handles.gray_radio,'Value')==1
    colormap gray
elseif get(handles.parula_radio,'Value')==1
    colormap parula
else
    colormap autumn  
end
% display option
if get(handles.imadjust_radio,'Value')==1
    imagesc(handles.image_axes,imadjust(transpose(squeeze(temp.data(:,temp.current_slice,:)))))
elseif get(handles.histeq_radio,'Value')==1
    imagesc(handles.image_axes,histeq(transpose(squeeze(temp.data(:,temp.current_slice,:)))))
else
    imagesc(handles.image_axes,adapthisteq(transpose(squeeze(temp.data(:,temp.current_slice,:)))))    
end

% plot points if there are any
if ~isempty(temp.temp_points)
    plot(handles.image_axes,temp.temp_points(:,1),temp.temp_points(:,2),'o','Color','k','MarkerSize',6,'MarkerFaceColor',lines(1));
end
set(handles.image_axes,'YDir','reverse','YLim',[.5 size(temp.data,3)+.5],'XLim',[.5 size(temp.data,1)+.5])

% now get points
try
    while size(temp.temp_points,1)<npts
        upoint=drawpoint(handles.image_axes,'InteractionsAllowed','none');
        x = upoint.Position(1);
        y = upoint.Position(2);
        temp.temp_points = cat(1,temp.temp_points,[x y]);
    end
catch exception 
end


function align_midsag_points
global temp
disp('aligning midsagittal plane...')

% let's make a black box with the midsag pts marked in white
% this will help us figure out the misag plane later
block = zeros(size(temp.data));
for i = 1:size(temp.midsag_points,1)
    block(round(temp.midsag_points(i,2)),round(temp.midsag_points(i,1)),round(temp.midsag_points(i,3)))=1;
end

% figure out the transformation
x = temp.midsag_points(:,1);
y = temp.midsag_points(:,2);
z = temp.midsag_points(:,3);
v = table(x,y,z,'VariableNames',{'x','y','z'});
mdlspec = 'y ~ x + z';
mdl = fitglm(v,mdlspec);
coef = [mdl.Coefficients.Estimate];
nrm = [coef(2) -1 coef(3)]; 
ref = [1 0 0];
% the rotation vector needed to get there
r = vrrotvec(ref,nrm); 

% turn it into a transform
tform = affine3d;
tform.T = rotm2tform(vrrotvec2mat(r));

% apply the transform
temp.data = imwarp(temp.data,tform);
block = imwarp(block,tform);

% do some trimming
colsum = sum(temp.data,[1 3]);
rowsum = sum(temp.data,[2 3]);
block = block(rowsum>0,colsum>0,:);
temp.data = temp.data(rowsum>0,colsum>0,:);

% now find the midsag plane
xsums = sum(block,[1 3]);
[~,midsag_plane]=max(xsums);
% now pad the image with zeros to make it symmetrical
hemi_width = max(abs([midsag_plane size(block,2)-midsag_plane]));
new_width = 2*hemi_width-1;
pad_width = new_width-size(block,2);
padding = zeros(size(block,1),pad_width,size(block,3));
if midsag_plane < hemi_width % pad left
    temp.data = cat(2,padding,temp.data);
else % pad right
    temp.data = cat(2,temp.data,padding);
end
disp('midsagittal plane aligned.')

function orient_volume_ML
global temp
% if the right is on the left
if temp.temp_points(1) < size(temp.data,2)/2 
    temp.data = imrotate(temp.data,180);
end

function orient_volume_AP
global temp
% if the front is at the back
if temp.temp_points(2) > size(temp.data,1)/2 
    temp.data = permute(temp.data,[1 3 2]);
    temp.data = imrotate(temp.data,180);
    temp.data = permute(temp.data,[1 3 2]);
end

function orient_volume_pitch
global temp
ang_data = atand(diff(temp.temp_points(:,2))/diff(temp.temp_points(:,1)));
ang_atlas = -12;
ang_diff = ang_atlas-ang_data;
temp.data = permute(temp.data,[1 3 2]);
temp.data = imrotate(temp.data,ang_diff);
temp.data = permute(temp.data,[1 3 2]);
% trim AP and DV edges (leave ML alone for symmetry)
rowsum = sum(temp.data,[2 3]);
zsum = sum(temp.data,[1 2]);
temp.data = temp.data(rowsum>0,:,zsum>0);



function save_tiff(data,filepath,finish)
disp('saving tiff...')

if contains(path,'Fast_Tiff_Write')
    fast_save_tif(filepath,data)
else
    disp('cannot find Fast_Tiff_Write in path. using bfsave, though this is slower.')
    bfsave(data,filepath,'dimensionOrder','XYTZC','BigTiff',true,'Compression', 'LZW');
end


disp('saved.')
if finish == 1
    closereq()
    clear global
end
