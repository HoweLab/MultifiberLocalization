% function output = atlas_labels(varargin)
% 
% This function returns atlas labels for a set of coordinates or atlas 
% indices. For example,
%
%       output = atlas_labels('coord',[0.8, 2.3, -3.4; 0.5, 2, -2]);
%       output = atlas_labels('idx',[445 800 413; 475 770 254]);
%
% Note that input coord or idx can have multiple rows, where each row is a
% point.
%
% Specifically, this function takes as input a m x 3 matrix of coordinates, 
% where m = the number of points, and the columns, in order, are AP, ML,
% DV. It uses the labels from the Allen Mouse Brain CCF Atlas and the Kim
% Lab Atlas, and returns a struct with the following fields:
% 
% coords        the input coordinates
% atlas_coords  a table with columns CCF_AP, CCF_AP_idx, CCF_ML,
%               CCF_ML_idx, CCF_DV, CCF_DV_idx, Chon_AP, Chon_AP_idx. The
%               columns labeled "idx" represent the index of the point in 
%               the atlas 3D matrix. The columns without "idx" represent
%               the coordinate in the atlas closest to the given
%               coordinate.
% ccf_labels    a table with all the labels from the Allen Mouse Brain CCF
%               Atlas corresponding to the given points
% chon_labels   a table with all the labels from the Kim Lab (Chon) Atlas
%               Atlas corresponding to the given points
%
%
%
% NOTE #1: This assumes that those atlases are in a folder labeled 
% MRIAtlas,located in the same parent directory as your DataAnalysis folder 
% (which contains this function). 
%
% go to https://github.com/HoweLab/MultifiberLocalization
% and follow the instructions to download atlas files and run the
% GENERATE_ATLAS_FILES App
%
%
%
% NOTE #2: all atlases are slightly different from one another, so
% this may not match exactly to the Franklin-Paxinos Atlas, for example.
% 
%
%
% References:
%
% Allen Mouse CCF Atlas:
% publication: https://pubmed.ncbi.nlm.nih.gov/32386544/
% more info: https://help.brain-map.org/display/mouseconnectivity/API#API-DownloadAtlas3-DReferenceModels
%
% Kim Lab (Chon) Atlas:
% publication: https://pubmed.ncbi.nlm.nih.gov/31699990/
% more info: https://kimlab.io/brain-map/atlas/
%
% Mai-Anh Vu, 12/13/2023
%

function output = atlas_labels(varargin)

% input parser
ip = inputParser;
ip.addParameter('coord',[]);
ip.addParameter('idx',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% turn off some warnings
warning('off','imageio:tiffmexutils:libtiffWarning')
warning('off','imageio:tiffmexutils:libtiffErrorAsWarning')
warning('off','imageio:tifftagsread:expectedTagDataFormatMultiple')


% MRI Atlas Stuff
mridir = strsplit(fileparts(mfilename('fullpath')),filesep);
mridir = fullfile(strjoin(mridir(1:end-1),filesep),'MRIAtlas');
ccf_key = load(fullfile(mridir,'CCF','ccf_key.mat'),'AP','ML','DV','labels','updated');
chon_key = load(fullfile(mridir,'Chon','chon_key.mat'),'AP','labels');

coord_str = {'AP','ML','DV'};
% atlas indices
if ~isempty(idx)
    AP_ML_DV = idx;     

% atlas_coordinates    
elseif ~isempty(coord)
    AP_ML_DV = coord;
    % now translate from coordinates to indices
    for j = 1:size(AP_ML_DV,1)
        for i = 1:3
            if ~isnan(AP_ML_DV(j,i))
                [~,this_idx] = min(abs(ccf_key.(coord_str{i}).(coord_str{i})-AP_ML_DV(j,i)));
                AP_ML_DV(j,i) = ccf_key.(coord_str{i}).slice(this_idx);
            end
        end
    end  

% missing inputs
else
    disp('You need to input coordinates or indices.')
    disp('Type "help atlas_slice" for more info.')
end


% output
output = struct;
output.coord = coord;
output.idx = idx;
% preallocate atlas coordinates
output.atlas_coord = table;
atlas_fields = {'CCF_AP','CCF_AP_idx',...
    'CCF_ML','CCF_ML_idx',...
    'CCF_DV','CCF_DV_idx',...
    'Chon_AP','Chon_AP_idx'};
for a = 1:numel(atlas_fields)
    output.atlas_coord.(atlas_fields{a}) = ...
        zeros(max([size(coord,1) size(idx,1)]),1);
end
output.ccf_labels = table;
output.chon_labels = table;

% figure out coord first
ccf_key_fields = {'AP','ML','DV'};
for c = 1:size(AP_ML_DV,1)
    for i = 1:numel(ccf_key_fields)
        % CCF coordinates
        %[~,idx] = min(abs(table2array(ccf_key.(ccf_key_fields{i})(:,2))-coord(c,i)));   
        idx = AP_ML_DV(c,i);
        output.atlas_coord.(['CCF_' ccf_key_fields{i}])(c) = ccf_key.(ccf_key_fields{i}){idx,2};
        output.atlas_coord.(['CCF_' ccf_key_fields{i} '_idx'])(c) = ccf_key.(ccf_key_fields{i}){idx,1};
        if i == 1 % Chon atlas has diff AP            
            [~,idx2] = min(abs(table2array(chon_key.(ccf_key_fields{i})(:,2))-output.atlas_coord.(['CCF_' ccf_key_fields{i}])(c)));
            output.atlas_coord.(['Chon_' ccf_key_fields{i}])(c) = chon_key.(ccf_key_fields{i}){idx2,2};
            output.atlas_coord.(['Chon_' ccf_key_fields{i} '_idx'])(c) = chon_key.(ccf_key_fields{i}){idx2,1};
        end
    end
end

% now loop over and get labels
output.chon_labels = chon_key.labels(1,:);
% first preallocate in case there are blank sections
f = fieldnames(output.chon_labels);
for i = 1:size(output.chon_labels,2)    
    if isnumeric(output.chon_labels.(f{i}))
        output.chon_labels.(f{i}) = -1;
    else 
        output.chon_labels.(f{i}) = {''};
    end
end
output.ccf_labels = ccf_key.labels(1,:);
f = fieldnames(output.ccf_labels);
for i = 1:size(output.ccf_labels,2)
    if isnumeric(output.ccf_labels.(f{i}))
        output.ccf_labels.(f{i}) = -1;
    else 
        output.ccf_labels.(f{i}) = {''};
    end
end
if size(AP_ML_DV,1)>1
    for c = 2:size(AP_ML_DV,1)
        output.chon_labels = vertcat(output.chon_labels,output.chon_labels(1,:));
        output.ccf_labels = vertcat(output.ccf_labels,output.ccf_labels(1,:));
    end
end
for c = 1:size(AP_ML_DV,1)
    labels_chon = tiffreadVolume(fullfile(mridir,'Chon','Chon_labels_coronal.tif'),...
        'pixelRegion',{[1 Inf],[1 Inf],[output.atlas_coord.Chon_AP_idx(c) output.atlas_coord.Chon_AP_idx(c)]});
    val_chon = labels_chon(output.atlas_coord.CCF_DV_idx(c),output.atlas_coord.CCF_ML_idx(c));        
    idx_chon = find(chon_key.labels.StructuralID==val_chon);
    %output.chon_labels = vertcat(output.chon_labels,chon_key.labels(idx_chon,:));
    if ~isempty(val_chon) && ~isempty(idx_chon)
        output.chon_labels(c,:) = chon_key.labels(idx_chon,:);
    end
    % ccf
    labels_ccf = tiffreadVolume(fullfile(mridir,'CCF','annotation_10_coronal.tif'),...
        'pixelRegion',{[1 Inf],[1 Inf],[output.atlas_coord.CCF_AP_idx(c) output.atlas_coord.CCF_AP_idx(c)]});
    val_ccf = labels_ccf(output.atlas_coord.CCF_DV_idx(c),output.atlas_coord.CCF_ML_idx(c));
    idx_ccf = find(ccf_key.labels.value==val_ccf);
    %output.ccf_labels = vertcat(output.ccf_labels,ccf_key.labels(idx_ccf,:));    
    if ~isempty(val_ccf) && ~isempty(idx_ccf)
        output.ccf_labels(c,:) = ccf_key.labels(idx_ccf,:);
    end
end


