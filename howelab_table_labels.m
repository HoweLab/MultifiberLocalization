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
% updated 7/9/2024 by Mai-Anh to add in_striatum if it isn't already there



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

%%% add in_striatum if it isn't in there
if ~isfield(output_table,'in_striatum')
%     % NOTE: this is just here for reference on how the values were
%     % generated
%     % NOTE: ccf_key and chon_key are the loaded files from
%     % the atlas folders
%     ccf_is_str = ccf_key.labels.parent_value == 477 | ...
%         ccf_key.labels.parent_value == 485 | ...
%         ccf_key.labels.parent_value == 493 | ...
%         ccf_key.labels.value == 477 | ...
%         ccf_key.labels.value == 485 |...
%         ccf_key.labels.value == 493;
%     ccf_str_values = ccf_key.labels.value(ccf_is_str);
%     chon_is_str = chon_key.labels.ParentID == 477 | ...
%         chon_key.labels.ParentID == 485 | ...
%         chon_key.labels.ParentID == 493 | ...
%         chon_key.labels.StructuralID == 477 | ...
%         chon_key.labels.StructuralID == 485 |...
%         chon_key.labels.StructuralID == 493 | ...
%         contains(chon_key.labels.Franklin_PaxinosFullName,'riatum') |...
%         contains(chon_key.labels.Franklin_PaxinosFullName,'ccumbens') | ...
%         contains(chon_key.labels.Franklin_PaxinosFullName,'utamen');
%     chon_str_values = chon_key.labels.StructuralID(chon_is_str);
    
    % let's make a striatum filter (these are hard-coded here for brevity,
    % but see commented section above for how they were derived
    ccf_str_values = [477 485 672 493 56 998 754 549009199 275 278];
    
    chon_str_values = [477 485 672 2376 2491 2294 2295 2296 2497 2395 ,...
        2297 2492 2498 2299 2298 2374 2380 2500 2302 2480 2483 2499 ,...
        2300 2301 2501 2479 2482 2481 2370 2496 2493 2484 2485 2486 ,...
        2494 2487 2489 2488 2490 2495 2001 2050 493 56 2074 2006 2007 ,...
        998 754 275 278 2013];
    
    output_table.in_striatum = zeros(size(output_table,1),1);
    for f = 1:size(output_table.in_striatum,1)
        output_table.in_striatum(f) = ...
            ismember(output_table.ccf_ID(f),ccf_str_values) |...
            ismember(output_table.chon_ID(f),chon_str_values);
    end
end

%%% output struct
output = struct;
output.table = output_table;
% save if save_path is supplied
if ~isempty(save_path)
    save(save_path,'-struct','output')
end



    
    
    
    