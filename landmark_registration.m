% function landmark_registration(ct_path,ct_landmarks_path,varargin)
%
% This function aligns the CT to the Allen Brain Institute CCF Atlas based
% on landmarks as referenced, marked, and detailed in the
% REGISTER_CT_README.pdf document
%
% Make sure you have the MRIAtlas folder downloaded and put in the same
% folder as CTProcessing. For example:
%     C:\Users\maianhvu\Documents\MATLAB\Scripts\CTProcessing
%     C:\Users\maianhvu\Documents\MATLAB\Scripts\MRIAtlas
%
% inputs:
% ct_path           - the path to the ct_coronal.tif file output from
%                     CT_dicom_align_and_convert
% ct_landmarks_path - the path to the .points file that has the landmarks
%                     you  identified using the FIJI plugin: Name Landmarks
%                     and Register
%
% optional inputs:
% output_name       - if left blank, the output name will be the same as
%                     the input name with _reg.tif and _reg_info.mat
%                     appended. default = blank
% reg_approach      - 2 possible approaches here:
%                       1.  regid registration (w/isotropic scaling), then 
%                           separate DV and AP scaling, then another rigid 
%                           registration (w/isotropic scaling)
%                       2.  rigid registration (w/isotropic scaling), then 
%                           affine registration
%                     Approach #1 works better if there aren't many
%                     landmarks, otherwise approach #2 works fine. When
%                     left blank, the approach will be determined by the
%                     number of landmarks identified and which ones.
%                     default = blank
%
% outputs:
% it outputs 2 files
% 1. a _____reg.tif that is your CT registered to the atlas
% 2. a _____reg_info.mat file that includes some info about the transforms
%
% Briefly, this algorithm registers the brain in 2 parts:
% 1. the ML scale and midsagittal alignment refinement (midsagittal plane 
%    was already established via CT_dicom_align_and_convert). This steps 
%    uses the non-midsagittal landmarks to estimate the scale, and an
%    average of the midsagittal landmarks and original midsagittal plane
%    estimation to re-assign the midsagittal plane.
%
% 2. the AP and DV pitch, and scale. This happens one of two ways:
%
%    The first way, which performs better when there aren't very many
%    landmarks identified, occurs in 3 steps. Note that I picked these 
%    steps over a straight affine registration, because the affine 
%    registration was overfitting the landmarks in some CTs and causing 
%    skewing and warping. This step uses the midsagittal landmarks.
%       i.      rigid registration of landmarks (basically to get the right
%               orientation)
%       ii.     estimation of scaling in DV and AP separately
%       iii.    another rigid registration of the resecaled landmarks
%
%    The second way, which performs better when there are enough landmarks
%    (at least 2 each of the more posterior nad more anterior midsagittal
%    landmarks, occurs in 2 steps. 
%       i.      rigid registration of landmarks (basically to get the right
%               orientation, with uniform/isotropic scaling)
%       ii.     affine registration of landmarks to refine the
%               registration, allowing nonuniform/anisotropic scaling
%
% Mai-Anh Vu
% 3/18/2022
% updated 4/7/22
% updated 4/10/22: documentation and final version  that doesn't overfit
% updated 9/2/2022: to make use of Fast_Tiff_Write if available
% updated 12/15/2022: to allow input of output name, as well as two
%   different approaches

function landmark_registration(ct_path,ct_landmarks_path,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('reg_approach',[]);
ip.addParameter('output_name',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end


%%% display a start message
disp('Registering your CT to the Allen Mouse Brain CCF Atlas. Will update when finished.')

%%% turn off some warnings
% turn off some warnings
warning('off','imageio:tiffmexutils:libtiffWarning')
warning('off','imageio:tiffmexutils:libtiffErrorAsWarning')
warning('off','imageio:tifftagsread:expectedTagDataFormatMultiple')

%%% load things
% CCF
MRIdir = strsplit(fileparts(mfilename('fullpath')),filesep);
MRIdir = fullfile(strjoin(MRIdir(1:end-1),filesep),'MRIAtlas');
atlas_landmarks = load_fiji_landmarks(fullfile(MRIdir,'CCF','landmarks.points'));
key = load(fullfile(MRIdir,'CCF','ccf_key.mat'),'landmark_notes');
key = key.landmark_notes;
tmp = imfinfo(fullfile(MRIdir,'CCF','average_template_10_coronal.tif'));
atlas_info = tmp(1);
atlas_info.Depth = numel(tmp);
atlas_midsag_ref2d = imref2d([atlas_info.Height,atlas_info.Depth],1,1); % midsag
atlas_midsag_ref3d = imref3d([atlas_info.Height,atlas_info.Width,atlas_info.Depth],1,1,1); %3D

clear tmp
% CT 
ct_landmarks = load_fiji_landmarks(ct_landmarks_path);
ct = load_tiffs_fast(ct_path);
ct_midsag_plane = ceil(size(ct,2)/2);

%%%%%%%%%%%%%
%%% sort ROIs
% midsagittal landmarks: AP/DV scaling, pitch
midsag_landmarks_all = key.name(key.midsag==1);
[~,midsag_ct,~] = intersect(ct_landmarks.name,midsag_landmarks_all);
midsag_ct = ct_landmarks(midsag_ct,:);
[~,midsag_atlas,~] = intersect(atlas_landmarks.name,midsag_landmarks_all);
midsag_atlas = atlas_landmarks(midsag_atlas,:);
% only keep the ones in common
[~,midsag_ct_idx,midsag_atlas_idx] = intersect(midsag_ct.name,midsag_atlas.name);
midsag_ct = midsag_ct(midsag_ct_idx,:);
midsag_atlas = midsag_atlas(midsag_atlas_idx,:);
% lateral extent landmarks: ML scaling
lat_landmarks_all = key.name(key.ML_extent==1);
[~,lat_ct,~] = intersect(ct_landmarks.name,lat_landmarks_all);
lat_ct = ct_landmarks(lat_ct,:);
[~,lat_atlas,~] = intersect(atlas_landmarks.name,lat_landmarks_all);
lat_atlas = atlas_landmarks(lat_atlas,:);
atlas_midsag_plane = imfinfo(fullfile(MRIdir,'CCF','average_template_10_coronal.tif'));
atlas_midsag_plane = atlas_midsag_plane(1).Width/2+.5;
% only keep the ones in common
[~,lat_ct_idx,lat_atlas_idx] = intersect(lat_ct.name,lat_atlas.name);
lat_ct = lat_ct(lat_ct_idx,:);
lat_atlas = lat_atlas(lat_atlas_idx,:);
% see if there's enough to do an affine transformation
if isempty(reg_approach)
    posterior_midsag = {'KM','PM','HM1','IP','CCM'};
    anterior_midsag = {'VM1','FM1','AC'};
    post_ant = [numel(intersect(midsag_ct.name,posterior_midsag)) numel(intersect(midsag_ct.name,anterior_midsag))];
    if sum(post_ant>=2) == 2
        reg_approach = 2;
    else
        reg_approach = 1;
    end
end
%%%%%%%%%%%%%%%%%%%
%%% transformations
% see here https://www.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html

% tform1: ML scale & shift
% ML centering: re-calculate midsagittal plane on CT
if size(ct,2)/2 < ct_midsag_plane % pad right
    pad_width = (2*ct_midsag_plane-1)-size(ct,2);
    ct = cat(2,ct,zeros(size(ct,1),pad_width,size(ct,3)));
elseif size(ct,2)/2 > ct_midsag_plane % pad left
    pad_width = size(ct,2)-(2*ct_midsag_plane-1);
    ct = cat(2,zeros(size(ct,1),pad_width,size(ct,3)),ct);
end
% ML rescale
ct_ML_rescale = mean(abs(lat_atlas.x-atlas_midsag_plane)./abs(lat_ct.x-ct_midsag_plane));
tform1 = affine3d;
tform1.T(1,1) = ct_ML_rescale;
ct = imwarp(ct,tform1);
% now trim to size symmetrically
MLmarg = abs(atlas_info.Width-size(ct,2))/2;
MLmarg = [ceil(MLmarg) floor(MLmarg)];
if size(ct,2)>atlas_info.Width
    ct = ct(:,MLmarg(2)+1:end-MLmarg(1),:);
elseif size(ct,2)<atlas_info.Width
    ct = cat(2,zeros(size(ct,1),MLmarg(1),size(ct,3)),...
        ct,zeros(size(ct,1),MLmarg(2),size(ct,3)));
end

% tform2: AP & DV scale, rotation, shift: rigid
tform2 = fitgeotrans([midsag_ct.z midsag_ct.y],[midsag_atlas.z midsag_atlas.y],'NonreflectiveSimilarity');
ct_midsag_tmp = transformPointsForward(tform2,[midsag_ct.z midsag_ct.y]);
ct = permute(imwarp(permute(ct,[1 3 2]),tform2,'OutputView',atlas_midsag_ref2d),[1 3 2]); % rotate CT pitch

if reg_approach == 2
    
    % tform3: AP & DV scale, rotation, shift: affine
    tform3 = fitgeotrans(ct_midsag_tmp,[midsag_atlas.z midsag_atlas.y],'affine');
    ct_midsag_tmp = transformPointsForward(tform3,ct_midsag_tmp);
    ct = permute(imwarp(permute(ct,[1 3 2]),tform3,'OutputView',atlas_midsag_ref2d),[1 3 2]); % rotate CT pitch
   
else % approach # 1
    
    % tform3: independent DV/AP scaling
    ct_DV_rescale = nanmean((midsag_atlas.y-min(midsag_atlas.y))./(ct_midsag_tmp(:,2)-min(ct_midsag_tmp(:,2))));
    ct_AP_rescale = nanmean((midsag_atlas.z-min(midsag_atlas.z))./(ct_midsag_tmp(:,1)-min(ct_midsag_tmp(:,1))));
    tform3 = affine2d;
    tform3.T(1,1) = ct_AP_rescale;
    tform3.T(2,2) = ct_DV_rescale;
    ct_midsag_tmp = transformPointsForward(tform3,ct_midsag_tmp);
    ct = permute(imwarp(permute(ct,[1 3 2]),tform3,'OutputView',atlas_midsag_ref2d),[1 3 2]); % rotate CT pitch

    % tform4: AP & DV scale, rotation, shift one more time, rigid again
    tform4 = fitgeotrans(ct_midsag_tmp,[midsag_atlas.z midsag_atlas.y],'NonreflectiveSimilarity');
    ct_midsag_tmp = transformPointsForward(tform4,ct_midsag_tmp);
    ct = permute(imwarp(permute(ct,[1 3 2]),tform4,'OutputView',atlas_midsag_ref2d),[1 3 2]); % rotate CT pitch
   
end
%%%%%%%%%%%%
%%% save out
% save out .tif and some info
if isempty(output_name)
    output_name = ct_path(1:end-4);
end
fTIF = Fast_BigTiff_Write([output_name '_reg.tif'],0,Tiff.Compression.LZW);    
for i = 1:size(ct,3)
    fTIF.WriteIMG(ct(:,:,i));
end
fTIF.close;
output = struct;
output.tform1 = tform1.T;
output.tform2 = tform2.T;
output.tform3 = tform3.T;
if reg_approach == 1
    output.tform4 = tform4.T;
else
    output.tform4 = [];
end
output.midsag_pts_tformed = ct_midsag_tmp;
output.midsag_atlas_pts = [midsag_atlas.z midsag_atlas.y];
save([output_name '_reg_info.mat'],'-struct','output')
%%% display an end message
disp('Finished.')

