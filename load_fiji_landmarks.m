% landmark_table = load_fiji_landmarks(points_file_path)
%
% This function takes in the path to a .points file created using the Fiji
% plugin Landmarks > Name Landmarks and Register and outputs a table
% containing the names and coordinates of those points.
%
% Mai-Anh Vu
% 3/16/2022

function landmark_table = load_fiji_landmarks(points_file_path)

%%% load in contents
fid = fopen(points_file_path); % open the file
content = textscan(fid,'%s'); % read contents in as a cell array of strings
content = content{1,1};
fclose(fid); % close the file

%%% extract info
% points_set
point_set = content(startsWith(content,'set='));
point_set = cellfun(@(s) s(strfind(s,'="')+2:end-1),point_set,'UniformOutput',false);
% names
point_names = content(startsWith(content,'name='));
point_names = cellfun(@(s) s(strfind(s,'="')+2:end-1),point_names,'UniformOutput',false);
point_names = point_names(contains(point_set,'true'));
% x
x_coords = content(startsWith(content,'x='));
x_coords = cell2mat(cellfun(@(s) str2double(s(strfind(s,'="')+2:end-1)),x_coords,'UniformOutput',false));
% y
y_coords = content(startsWith(content,'y='));
y_coords = cell2mat(cellfun(@(s) str2double(s(strfind(s,'="')+2:end-1)),y_coords,'UniformOutput',false));
% z
z_coords = content(startsWith(content,'z='));
z_coords = cell2mat(cellfun(@(s) str2double(s(strfind(s,'="')+2:strfind(s,'"/')-1)),z_coords,'UniformOutput',false));

%%% fill in table
landmark_table = table(point_names,'VariableNames',{'name'});
landmark_table.x = x_coords;
landmark_table.y = y_coords;
landmark_table.z = z_coords;




