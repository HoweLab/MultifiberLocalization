function output = howelab_table_labels(input_table,varargin)
% Change the fieldnames of the table to ones matching Howe Lab legacy code 
% this returns a struct "output" with a single field "table" that has the
% table in it.
%
% This requires a single input, which can be the new fiber table itself,
% the struct with the table in it as the field "fiber_table", or the path
% to a struct with the table in it as the field "fiber_table".
%
% Optional inputs:
% save_path:    if supplied, the output struct will be saved there
%
% Mai-Anh 4/12/2024



%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('save_path',[]);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

%%% make something like a dictionary for fieldnames
fieldmap = struct;
fieldmap.ID = 'ROI';
fieldmap.top_ID = 'top_ID';
fieldmap.bottom_ID = 'bottom_ID';
fieldmap.grid_row = 'calib_row';
fieldmap.grid_col = 'calib_col';
fieldmap.in_striatum = 'in_striatum';
fieldmap.is_multifiber = 'is_multifiber';
fieldmap.ccf_ID = 'ccf_ID';
fieldmap.chon_ID = 'chon_ID';
fieldmap.ccf_label = 'ccf_label';
fieldmap.chon_allen_label = 'allen_label';
fieldmap.chon_allen_label_abbrev = 'allen_label_abbrev';
fieldmap.chon_fp_label = 'fp_label';
fieldmap.chon_fp_label_abbrev = 'fp_label_abbrev';
coords = {'AP','ML','DV'};
which_end = {'bottom','top'};
which_coordinate_system = {'','_idx'};
for c = 1:numel(coords)
    for e = 1:numel(which_end)
        for s = 1:numel(which_coordinate_system)            
            fieldmap.([which_end{e} '_' coords{c} which_coordinate_system{s}]) = ...
                ['fiber_' which_end{e} '_' coords{c} which_coordinate_system{s}];
        end
    end
end
fieldmap.bottom_DV2 = 'fiber_bottom_DV2'; % then add this one
    

%%% load input table 
% if is char, let's assume it's a path
if ischar(input_table) 
    input_table = load(input_table);
end
% if it's a struct, let's grab the fiber_table field
if isstruct(input_table) 
    input_table = input_table.fiber_table;
end
% get the fieldnames
input_table_fields = fieldnames(input_table);
fields_ignore = {'Properties','Row','Variables'};
input_table_fields = input_table_fields(~ismember(input_table_fields,fields_ignore));

%%% now make output table
output_table = table;
for f = 1:numel(input_table_fields)
    if isfield(fieldmap, input_table_fields{f})
        output_table.(fieldmap.(input_table_fields{f})) = input_table.(input_table_fields{f});
    else
        output_table.(input_table_fields{f}) = input_table.(input_table_fields{f});
    end
end

%%% output struct
output = struct;
output.table = output_table;
% save if save_path is supplied
if ~isempty(save_path)
    save(save_path,'-struct','output')
end



    
    
    
    