
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAP FIBER TIPS - MAIN %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = map_fiber_tips(output,varargin)

%%%  parse optional inputs %%%
ip = inputParser;
ip.addParameter('watershed_sep',1); % 0 for GMM, 1 for watershed separation
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end


% slide our detected fiber tops into a zero matrix: this will help us line
% it up with fiber blobs
fiber_tops = zeros(output.info.ct_size);
fiber_tops(output.fiber_tops.axial_slice_num,:,:) = ...
    permute(output.fiber_tops.tops_mask,[3 2 1]);

% make an adaptive mask
fibers_mask = zeros(output.info.ct_size);
for i = output.fiber_tops.axial_slice_num:output.info.ct_size(1)    
    % top 3 slices: mask = fiber tops
    if i<output.fiber_tops.axial_slice_num+3
        bw = fiber_tops(output.fiber_tops.axial_slice_num,:,:)>0;        
    else
        this_slice = output.ct(i,:,:);
        this_slice(~output.fiber_bottoms.axial_proj.polymask(i,:,:)) = 0;
        t = adaptthresh(this_slice); % adaptive threshold    
        bw = imbinarize(this_slice,max(t(:))); % use the max
    end
    fibers_mask(i,:,:) = bw;
end

% blob detection and properties
detected_fibers = bwconncomp(fibers_mask,26); % default connectivity
fiber_props = regionprops(detected_fibers);
output.fiber_bottoms.boxes = round(transpose(cell2mat(cellfun(@(x) transpose(x),{fiber_props.BoundingBox},'UniformOutput',false))));
output.fiber_bottoms.voxels = transpose(detected_fibers.PixelIdxList);
output.fiber_bottoms.bottom_ID = transpose(1:size(output.fiber_bottoms.boxes,1));

% now line them up with fiber tops
fiber_topID = cell(numel(output.fiber_bottoms.bottom_ID),1);
for i = 1:numel(output.fiber_bottoms.bottom_ID)
    blobValues = unique(fiber_tops(detected_fibers.PixelIdxList{i}));
    fiber_topID{i} = blobValues(blobValues>0);
end
output.fiber_bottoms.top_ID = fiber_topID;

% now sort our fibers: first get all the top info in one place
fiber_table = output.fiber_tops.tops_table;
[~,idx] = sort(fiber_table.top_ID);
fiber_table = fiber_table(idx,:);
tmp_nums = 1:size(output.fiber_tops.centers_AP);
fiber_table.top_AP_idx = round(output.fiber_tops.centers_AP(ismember(tmp_nums,fiber_table.top_ID)));
fiber_table.top_ML_idx = round(output.fiber_tops.centers_ML(ismember(tmp_nums,fiber_table.top_ID)));
fiber_table.top_DV_idx = ones(size(fiber_table.top_AP_idx))*output.fiber_tops.axial_slice_num;
% now start to fill in the bottom
fiber_table.bottom_ID = zeros(size(fiber_table.top_ID));
fiber_table.bottom_AP_idx = zeros(size(fiber_table.top_ID));
fiber_table.bottom_ML_idx = zeros(size(fiber_table.top_ID));
fiber_table.bottom_DV_idx = zeros(size(fiber_table.top_ID));

% fibers that line up with a single top (these are the best!)
fiber_lists = struct;
fiber_lists.fiber_single_top = find(cellfun(@(x) numel(x)==1,fiber_topID));
% map these to the table
for f = 1:numel(fiber_lists.fiber_single_top)
    fibnum = fiber_lists.fiber_single_top(f);
    fibtop = fiber_topID{fibnum};
    fiber_table.bottom_ID(fiber_table.top_ID==fibtop) = fibnum;
    % now add in the location of the bottom (recording location)
    % this is the center of mass of the mask bottom slice
    [dv,ml,ap] = ind2sub(size(output.ct),output.fiber_bottoms.voxels{fibnum});
    dv_end = max(dv);
    fiber_table.bottom_ML_idx(fiber_table.top_ID==fibtop) = round(mean(ml(dv==dv_end)));
    fiber_table.bottom_AP_idx(fiber_table.top_ID==fibtop) = round(mean(ap(dv==dv_end)));
    fiber_table.bottom_DV_idx(fiber_table.top_ID==fibtop) = dv_end;
end

% fibers that line up with multiple tops (these will get separated)
fiber_lists.fiber_multiple_tops = find(cellfun(@(x) numel(x)>1,fiber_topID));
accounted_tops = [cell2mat(fiber_topID(fiber_lists.fiber_single_top));...
    cell2mat(fiber_topID(fiber_lists.fiber_multiple_tops))];
[~,idx,~] = unique(accounted_tops);
% save this in our table
fiber_table.is_multifiber = zeros(size(fiber_table.top_ID));
fiber_table.is_multifiber(ismember(fiber_table.bottom_ID,fiber_lists.fiber_multiple_tops)) = 1;

% fibers that line up with no top (either non-fibers, or disconnected fibers)
fiber_lists.fiber_no_top = find(cellfun(@isempty,fiber_topID));

% tops not associated with any fibers(disconnected thresholding)
fiber_lists.top_no_fibers = setdiff(fiber_table.top_ID,accounted_tops);

% tops that line up with more than one fiber (VERY unlikely)
fiber_lists.top_multiple_fibers = accounted_tops; 
fiber_lists.top_multiple_fibers(idx) = [];

% add to our output struct
output.fiber_bottoms.fiber_lists = fiber_lists;
output.fiber_table = fiber_table;

% separate
if ~isempty(fiber_lists.fiber_multiple_tops)
    for i = 1:numel(fiber_lists.fiber_multiple_tops)
        current_multifiber = fiber_lists.fiber_multiple_tops(i);        
        output = separate_multifibers(current_multifiber,output,'watershed_sep',watershed_sep);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MULTIFIBER SEPARATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = separate_multifibers(bottom_fibnum,output,varargin)

%%%  parse optional inputs %%%
% watershed_sep = 1: uses watershed to separate multifibers
% watershed_sep = 0: uses gaussian mixture model to separate multifibers

ip = inputParser;
ip.addParameter('watershed_sep',0);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end
% suppress GMM warnings
warning('off','stats:gmdistribution:FailedToConverge')

% how many fibers are we separating 
%(i.e., how many tops are included in this this multifiber)?
n_tops = numel(output.fiber_bottoms.top_ID{bottom_fibnum});

% make a small box to work with, to be memory-efficient
[dv,ml,ap] = ind2sub(output.info.ct_size,output.fiber_bottoms.voxels{bottom_fibnum});
this_box = output.fiber_bottoms.boxes(bottom_fibnum,:);
dv_box = dv-this_box(2)+1;
ml_box = ml-this_box(1)+1;
ap_box = ap-this_box(3)+1;
dv_slices = unique(dv_box); 

% make a fiber mask in the small box (where the fiber blob is)
this_fiber_mask = zeros(this_box(5), this_box(4), this_box(6));
this_fiber_mask(sub2ind(size(this_fiber_mask),dv_box,ml_box,ap_box)) = 1;
% now our CT, but with the non-blob masked out
% output.ct: rows = DV, cols = ML, slices = AP
this_fiber_ct = output.ct(...
    this_box(2):(this_box(2)+this_box(5)-1),...
    this_box(1):(this_box(1)+this_box(4)-1),...
    this_box(3):(this_box(3)+this_box(6)-1));
this_fiber_ct(this_fiber_mask==0) = 0;

% map of the relevant fiber tops: and make it match the top of the small box
% otput.fiber_tops.tops_mask: rows == AP, cols = ML
this_fiber_top = output.fiber_tops.tops_mask(...
    this_box(3):(this_box(3)+this_box(6)-1),...
    this_box(1):(this_box(1)+this_box(4)-1));
this_fiber_top(~ismember(this_fiber_top,output.fiber_bottoms.top_ID{bottom_fibnum})) = 0;

if watershed_sep == 1
    % go slice by slice and use the watershed approach to separate blobs and
    % extract centers
    new_mask = zeros(size(this_fiber_mask));
    new_mask(1:10,:,:) = this_fiber_mask(1:10,:,:);
    blob_info = table(transpose(1:numel(dv_slices)),'VariableNames',{'dv_box'});
    max_blob_area = 0;
    for d = 1:max(dv_slices)
        bw = squeeze(this_fiber_mask(d,:,:)); % orientation: ML (rows,y) AP (cols,x)        
        if d<=10 % don't need to separate first 10 slices
            blobs = bwconncomp(bw);
            blobprops = regionprops(blobs);
            if max(cellfun(@(x) x,{blobprops.Area})) > max_blob_area
                max_blob_area = max(cellfun(@(x) x,{blobprops.Area}));
            end
        else % watershed for slices 11 onwards
            D = -bwdist(~bw);
            Ld = watershed(D,8); 
            bw(Ld==0)=0;
            blobs = bwconncomp(bw);
            blob_sizes = cellfun(@numel,blobs.PixelIdxList);    
            if sum(blob_sizes>max_blob_area*1.5)>0 % if any blobs are too big, split again
                D =-bwdist(~bw);
                Ld = watershed(D,8); 
                bw(Ld==0)=0;
                blobs = bwconncomp(bw);                
            end
            blobprops = regionprops(blobs);
        end
        blob_info.n(d) = blobs.NumObjects;
        blob_info.voxels{d} = blobs.PixelIdxList;
        blob_info.centers(d) = {transpose(cell2mat(cellfun(@transpose,{blobprops.Centroid},'UniformOutput',false)))};  
        new_mask(d,:,:) = bw;
    end
elseif watershed_sep == 0
    % go slice by slice and fit a 2d gaussian mixture model to the blobs
    new_mask = zeros(size(this_fiber_mask));
    new_mask(1:10,:,:) = this_fiber_mask(1:10,:,:);
    blob_info = table(transpose(1:numel(dv_slices)),'VariableNames',{'dv_box'});
    for d = 1:max(dv_slices)
        if d<=10 % first 10 slices just leave
            bw = squeeze(this_fiber_mask(d,:,:)); % orientation: ML (rows,y) AP (cols,x)        
            blobs = bwconncomp(bw);
            blobprops = regionprops(blobs);
            blob_info.n(d) = blobs.NumObjects;
            blob_info.voxels{d} = blobs.PixelIdxList;
            blob_info.centers(d) = {transpose(cell2mat(cellfun(@transpose,{blobprops.Centroid},'UniformOutput',false)))};  
            new_mask(d,:,:) = bw;
        else
            this_slice = this_fiber_ct(d,:,:);
            % generate something like weighted point clouds to fit gaussian
            % mixture model
            [~,these_x,these_z] = ind2sub(size(this_slice),find(this_slice>0));
            these_y = double(this_slice(this_slice>0));
            these_y = round(100*these_y/max(these_y));
            these_y = these_y-min(these_y)+1;
            rep_y = cell2mat(arrayfun(@(x) [ones(1,x) nan(1,max(these_y)-x)],these_y,'UniformOutput',false));        
            rep_x = repmat(these_x,1,max(these_y)).*rep_y;
            rep_x = rep_x(:);
            rep_x = rep_x(rep_x>0);
            rep_z = repmat(these_z,1,max(these_y)).*rep_y;
            rep_z = rep_z(:);
            rep_z = rep_z(rep_z>0);
            if numel(rep_x) > 2 % needs at least 3 points
                for i = 1:n_tops
                    %this_gmm = fitgmdist([these_x, these_z],i);
                    %this_gmm = fitgmdist([these_x, these_z, these_y],i);
                    this_gmm = fitgmdist([rep_x, rep_z],i);
                    best_nll = inf;
                    if this_gmm.NegativeLogLikelihood < best_nll
                        best_model = this_gmm;
                        best_n = i;
                        best_nll = this_gmm.NegativeLogLikelihood;
                        best_centroids = this_gmm.mu;
                    end
                end
                % get the separated fiber identities
                gmm_clust = best_model.cluster([these_x these_z]);
                gmm_clust_vox = cell(1,best_n);
                bw = zeros(size(this_slice));
                for i = 1:best_n
                    clust_x = these_x(gmm_clust==i);
                    clust_z = these_z(gmm_clust==i);
                    clust_idx = sub2ind(size(squeeze(this_slice)),clust_x,clust_z);
                    gmm_clust_vox{i} = clust_idx;
                    tmp_bw = zeros(size(squeeze(bw)));
                    tmp_bw(clust_idx) = 1;
                    tmp_bw = imerode(tmp_bw,strel('disk',1));
                    bw(tmp_bw == 1) = 1;
                end
                blob_info.n(d) = best_n;
                blob_info.voxels{d} = gmm_clust_vox;
                blob_info.centers{d} = fliplr(best_centroids);
                new_mask(d,:,:) = bw;
            else % just assign them based on distance 
                prev_centroids = blob_info.centers{d-1};
                tmp_pts = fliplr([rep_x rep_z]);
                tmp_dist = pdist2(tmp_pts, prev_centroids);
                [~,tmp_idx] = min(tmp_dist,[],2);
                blob_info.n(d) = numel(unique(tmp_idx));
                if numel(unique(tmp_idx)) > 1
                    blob_info.centers{d} = tmp_pts;                   
                    tmp_pts = find(this_slice>0); 
                    blob_info.voxels{d}{1} = tmp_pts(1);
                    blob_info.voxels{d}{2} = tmp_pts(2);
                else
                    blob_info.centers{d} = round(mean(tmp_pts));                                        
                    blob_info.voxels{d}{1} = find(this_slice>0); 
                end
                
            end
        end
    end
end


%%%% BUILD UP FIBERS SLICE BY SLICE %%%%
% fit a line to the data
fiber_lines = struct;
last_vox = struct;
% loop over dorsal-ventral slices
for d = 1:size(blob_info,1)
    % if it's the first one, just assign and note the voxel mask
    if  d == 1     
        for f = 1:n_tops
            fiber_lines.(['line' num2str(f)]) = [blob_info.centers{d}(f,:) d];
            last_vox.(['line' num2str(f)]) = blob_info.voxels{d}{f};
        end

    else
        available_lines = fieldnames(fiber_lines);
        current_points = [blob_info.centers{d} repmat(d,size(blob_info.centers{d},1),1)];
        current_points_idx = sub2ind([size(this_fiber_mask,2) size(this_fiber_mask,3)],round(current_points(:,2)),round(current_points(:,1)));
        % see if any new points overlap with the blob masks from the last slice
        for i_line = 1:numel(available_lines)
            % make a convex hull around last slice's mask just in case
            % there's a split in it or something
            these_vox = last_vox.(['line' num2str(i_line)]);[these_r,these_c] = ind2sub(size(this_fiber_mask,2,3),these_vox);
            try
                this_hull = convhull(these_c,these_r);
                these_vox = find(roipoly(zeros(size(this_fiber_mask,2,3)),these_c(this_hull),these_r(this_hull))==1);               
            catch exception
            end
            [~,~,i_pt] = intersect(these_vox,current_points_idx);
            if ~isempty(i_pt)
                % if more than one point fit
                if numel(i_pt)>1
                    % combine all relevant blobs for next voxel mask
                    last_vox.(['line' num2str(i_line)]) = [];
                    for i = 1:numel(i_pt)
                        last_vox.(['line' num2str(i_line)]) = [...
                            last_vox.(['line' num2str(i_line)]); blob_info.voxels{d}{i_pt(i)}];
                    end
                    if d >= 10
                        % if we have enough of a line going, 
                        % see what fits better, one of them, or the average?
                        this_line = fiber_lines.(available_lines{i_line});
                        exy = prctile(this_line,[25 75]);
                        % fix by Alice in case fiber is too short
                        if size(this_line,1)<=3
                            this_line = [this_line; this_line];
                        end
                        this_fit = polyfitn(this_line(:,1:2),this_line(:,3),{'constant', 'x', 'y'});
                        exy(:,3) = this_fit.Coefficients(1)+this_fit.Coefficients(2)*exy(:,1)+this_fit.Coefficients(3)*exy(:,2);                  
                        candidate_points = [mean(current_points(i_pt,:));...
                            current_points(i_pt,:)];
                        fit_check = nan(size(candidate_points,1),1);
                        for i = 1:numel(fit_check)
                            fit_check(i) = point_to_line_distance(candidate_points(i,:),exy(1,:),exy(2,:));
                        end
                        [~,i_pt] = min(fit_check);
                        current_point = candidate_points(i_pt,:);
                    else
                        current_point = mean(current_points(i_pt,:));
                    end
                else
                    last_vox.(['line' num2str(i_line)]) = blob_info.voxels{d}{i_pt};
                    current_point = current_points(i_pt,:);
                end
                fiber_lines.(['line' num2str(i_line)]) = [...
                    fiber_lines.(['line' num2str(i_line)]);...
                    current_point];
            end
        end                       
    end
end
available_lines = fieldnames(fiber_lines);

%%%% SORT SHARED POINTS BETW LINES %%%
for i = 1:numel(available_lines)
    line_i = available_lines{i};
    for j = (i+1):numel(available_lines)
        line_j = available_lines{j};
        shared_points = intersect(fiber_lines.(line_i),fiber_lines.(line_j),'rows');
        if ~isempty(shared_points)
            % the unshared points: fit a line and get 2 example points
            line_i_unshared = setdiff(fiber_lines.(line_i),shared_points,'rows');
            exy_i = prctile(line_i_unshared,[30 60]);
            fit_i = polyfitn(line_i_unshared(:,1:2),line_i_unshared(:,3),{'constant', 'x', 'y'});
            exy_i(:,3) = fit_i.Coefficients(1)+fit_i.Coefficients(2)*exy_i(:,1)+fit_i.Coefficients(3)*exy_i(:,2);
            % the unshared points: fit a line and get to example points
            line_j_unshared = setdiff(fiber_lines.(line_j),shared_points,'rows');
            exy_j = prctile(line_j_unshared,[30 60]);
            fit_j = polyfitn(line_j_unshared(:,1:2),line_j_unshared(:,3),{'constant', 'x', 'y'});
            exy_j(:,3) = fit_j.Coefficients(1)+fit_j.Coefficients(2)*exy_j(:,1)+fit_j.Coefficients(3)*exy_j(:,2);               
            for k = 1:size(shared_points,1)
                % distance to line i
                dist_i = point_to_line_distance(shared_points(k,:),exy_i(1,:),exy_i(2,:));
                % distance to line j
                dist_j = point_to_line_distance(shared_points(k,:),exy_j(1,:),exy_j(2,:));
                % add it to the closer line
                if dist_i<dist_j
                    line_i_unshared = [line_i_unshared; shared_points(k,:)];
                else
                    line_j_unshared = [line_j_unshared; shared_points(k,:)];
                end
            end
            [~,sort_idx] = sort(line_i_unshared(:,3));
            fiber_lines.(line_i) = line_i_unshared(sort_idx,:);
            [~,sort_idx] = sort(line_j_unshared(:,3));
            fiber_lines.(line_j) = line_j_unshared(sort_idx,:);
        end
    end
end

%%%% IDENTIFY FIBER TOPS AND TIPS %%%%
fiber_tops = struct;
fiber_ends = struct;
for i = 1:numel(available_lines)
    this_top = this_fiber_top(...
        round(fiber_lines.(available_lines{i})(1,1)),... % ML = X
        round(fiber_lines.(available_lines{i})(1,2))); % AP = Y
    this_tip = round(fiber_lines.(available_lines{i})(end,:));
    fiber_tops.(available_lines{i}) = this_top; % store this info
    if sum(output.fiber_table.top_ID==this_top) > 0 && this_top > 0
        fiber_ends.(available_lines{i}) = [this_tip;...
            output.fiber_table.top_AP_idx(output.fiber_table.top_ID==this_top)-this_box(3)+1,...
            output.fiber_table.top_ML_idx(output.fiber_table.top_ID==this_top)-this_box(1)+1 1];
        % now add it to our table
        output.fiber_table.bottom_ID(output.fiber_table.top_ID==this_top) = bottom_fibnum;
        output.fiber_table.bottom_ML_idx(output.fiber_table.top_ID==this_top) = round(this_tip(2)+this_box(1)-1);
        output.fiber_table.bottom_AP_idx(output.fiber_table.top_ID==this_top) = round(this_tip(1)+this_box(3)-1);
        output.fiber_table.bottom_DV_idx(output.fiber_table.top_ID==this_top) = round(this_tip(3)+this_box(2)-1);
    end
end

output.fiber_bottoms.multi_fiber_checks.(['fibnum' num2str(bottom_fibnum)]).fibers = fiber_lines;
output.fiber_bottoms.multi_fiber_checks.(['fibnum' num2str(bottom_fibnum)]).ends = fiber_ends;
output.fiber_bottoms.multi_fiber_checks.(['fibnum' num2str(bottom_fibnum)]).top_nums = fiber_tops;
output.fiber_bottoms.multi_fiber_checks.(['fibnum' num2str(bottom_fibnum)]).new_mask = find(new_mask==1);
