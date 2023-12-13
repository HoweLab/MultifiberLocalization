% function output = atlas_labels(coords)
% 
% This function returns atlas labels for a set of coordinates.
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
% Mai-Anh Vu, 11/30/2023
%

function output = atlas_labels(coords)


% MRI Atlas Stuff
mridir = strsplit(fileparts(mfilename('fullpath')),filesep);
mridir = fullfile(strjoin(mridir(1:end-1),filesep),'MRIAtlas');
ccf_key = load(fullfile(mridir,'CCF','ccf_key.mat'),'AP','ML','DV','labels','updated');
chon_key = load(fullfile(mridir,'Chon','chon_key.mat'),'AP','labels');

% output
output = struct;
output.coords = coords;
% preallocate atlas coordinates
output.atlas_coords = table;
atlas_fields = {'CCF_AP','CCF_AP_idx',...
    'CCF_ML','CCF_ML_idx',...
    'CCF_DV','CCF_DV_idx',...
    'Chon_AP','Chon_AP_idx'};
for a = 1:numel(atlas_fields)
    output.atlas_coords.(atlas_fields{a}) = zeros(size(coords(:,1)));
end
output.ccf_labels = table;
output.chon_labels = table;

% figure out coords first
ccf_key_fields = {'AP','ML','DV'};
for c = 1:size(coords,1)
    for i = 1:numel(ccf_key_fields)
        % CCF coordinates
        [~,idx] = min(abs(table2array(ccf_key.(ccf_key_fields{i})(:,2))-coords(c,i)));   
        output.atlas_coords.(['CCF_' ccf_key_fields{i}])(c) = ccf_key.(ccf_key_fields{i}){idx,2};
        output.atlas_coords.(['CCF_' ccf_key_fields{i} '_idx'])(c) = ccf_key.(ccf_key_fields{i}){idx,1};
        if i == 1 % Chon atlas has diff AP
            [~,idx2] = min(abs(table2array(chon_key.(ccf_key_fields{i})(:,2))-coords(c,i)));
            output.atlas_coords.(['Chon_' ccf_key_fields{i}])(c) = chon_key.(ccf_key_fields{i}){idx2,2};
            output.atlas_coords.(['Chon_' ccf_key_fields{i} '_idx'])(c) = chon_key.(ccf_key_fields{i}){idx2,1};
        end
    end
end

% now loop over and get labels
for c = 1:size(coords,1)
    labels_chon = tiffreadVolume(fullfile(mridir,'Chon','Chon_labels_coronal.tif'),...
        'pixelRegion',{[1 Inf],[1 Inf],[output.atlas_coords.Chon_AP_idx(c) output.atlas_coords.Chon_AP_idx(c)]});
    val_chon = labels_chon(output.atlas_coords.CCF_DV_idx(c),output.atlas_coords.CCF_ML_idx(c));        
    idx_chon = find(chon_key.labels.StructuralID==val_chon);
    output.chon_labels = vertcat(output.chon_labels,chon_key.labels(idx_chon,:));
    % ccf
    labels_ccf = tiffreadVolume(fullfile(mridir,'CCF','annotation_10_coronal.tif'),...
        'pixelRegion',{[1 Inf],[1 Inf],[output.atlas_coords.CCF_AP_idx(c) output.atlas_coords.CCF_AP_idx(c)]});
    val_ccf = labels_ccf(output.atlas_coords.CCF_DV_idx(c),output.atlas_coords.CCF_ML_idx(c));
    idx_ccf = find(ccf_key.labels.value==val_ccf);
    output.ccf_labels = vertcat(output.ccf_labels,ccf_key.labels(idx_ccf,:));    
end


