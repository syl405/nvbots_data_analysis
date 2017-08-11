%clear

%% Import Metrology Data R1
% list directory contents
csv_dir_path = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAPFIXCORR\parsed_csv\';
dir_struct = dir(csv_dir_path);

% initialize data structure for imported CMM data
C_R1 = struct('data',{},'run',{}, 'coupling_pts', {});

% iterate through directory list to import data
for i = 1:numel(dir_struct)
    if dir_struct(i).isdir %skip directories
        continue
    elseif ~endsWith(dir_struct(i).name,'.csv') %handle non-text files
        error('Non-CSV file encountered in directory.')
    end
    
    filename = [dir_struct(i).folder '\' dir_struct(i).name];
    run_index = regexp(dir_struct(i).name,'MFC\d?_R(\d+)_parsed.csv','tokenExtents');
    run = str2double(dir_struct(i).name(run_index{1}(1):run_index{1}(2)));
    
    delimiter = ',';
    startRow = 2;
    
    % Format for each line of text:
    %   column1: categorical (%C)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%C%f%f%f%f%[^\n\r]';
    
    % Open the text file.
    fileID = fopen(filename,'r');
    
    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    % Close the text file.
    fclose(fileID);
    
    % Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.
    
    % Create output variable
    metrology_data = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Flatness','Z','XYParallelism'});
    
    % Renumber feature numbers to be in sequence
    metrology_data = sortrows(metrology_data,{'FeatureNumber'},{'ascend'}); %sort by feature number
    metrology_data.FeatureNumber = (1:84)';
    
    % Separate column data coupling error measurements
    coupling_point_data = metrology_data(size(metrology_data,1)-2:end,:); % coupling data (3 planes)
   	coupling_point_data.FeatureNumber = (1:3)'; % renumber kc points 1-3
    
    metrology_data = metrology_data(1:size(metrology_data,1)-3,:); % column data (81 planes)
    
    % Apply nominal offset to reference height against bottom of part
    metrology_data.Z = metrology_data.Z - 1.83;
    
    % Clear temporary variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    % Create new entry in C with newly imported measurements
    C_R1 = [C_R1; struct('data', metrology_data, 'coupling_pts', coupling_point_data, 'run', run)];
end

%% Import Metrology Data R2
% list directory contents
csv_dir_path = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAPFIXCORR2\parsed_csv\';
dir_struct = dir(csv_dir_path);

% initialize data structure for imported CMM data
C_R2 = struct('data',{},'run',{}, 'coupling_pts', {});

% iterate through directory list to import data
for i = 1:numel(dir_struct)
    if dir_struct(i).isdir %skip directories
        continue
    elseif ~endsWith(dir_struct(i).name,'.csv') %handle non-text files
        error('Non-CSV file encountered in directory.')
    end
    
    filename = [dir_struct(i).folder '\' dir_struct(i).name];
    run_index = regexp(dir_struct(i).name,'MFC\d?_R(\d+)_parsed.csv','tokenExtents');
    run = str2double(dir_struct(i).name(run_index{1}(1):run_index{1}(2)));
    
    delimiter = ',';
    startRow = 2;
    
    % Format for each line of text:
    %   column1: categorical (%C)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%C%f%f%f%f%[^\n\r]';
    
    % Open the text file.
    fileID = fopen(filename,'r');
    
    % Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    % Close the text file.
    fclose(fileID);
    
    % Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.
    
    % Create output variable
    metrology_data = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Flatness','Z','XYParallelism'});
    
    % Renumber feature numbers to be in sequence
    metrology_data = sortrows(metrology_data,{'FeatureNumber'},{'ascend'}); %sort by feature number
    metrology_data.FeatureNumber = (1:84)';
    
    % Separate column data coupling error measurements
    coupling_point_data = metrology_data(size(metrology_data,1)-2:end,:); % coupling data (3 planes)
   	coupling_point_data.FeatureNumber = (1:3)'; % renumber kc points 1-3
    
    metrology_data = metrology_data(1:size(metrology_data,1)-3,:); % column data (81 planes)
    
    % Apply nominal offset to reference height against bottom of part
    metrology_data.Z = metrology_data.Z - 1.83;
    
    % Clear temporary variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    
    % Create new entry in C with newly imported measurements
    C_R2 = [C_R2; struct('data', metrology_data, 'coupling_pts', coupling_point_data, 'run', run)];
end

%% Import Design Data
% Script for importing data from the following spreadsheet:
%
%    Workbook: D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\artifact_R1 - Variable Pattern Table - VarPattern.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/21 14:04:43

% Import the data
[~, ~, raw] = xlsread('D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\artifact_R1 - Variable Pattern Table - VarPattern.xlsx','Sheet1','A3:E83');

% Create output variable
data = reshape([raw{:}],size(raw));

% Create table
design_data = table;

% Allocate imported array to column variable names
design_data.FeatureNumber = (1:81)';
design_data.DesignX = data(:,3);
design_data.DesignY = 200 - data(:,4);
design_data.DesignZ = data(:,5);

% Offset to reference off of bottom face of part
design_data.DesignZ = design_data.DesignZ + 3.9;

% Clear temporary variables
clearvars data raw;

rows = {[1 1 1 1 1 1 1 1 1];...
        [1 0 1 0 1 0 1 0 1];...
        [1 0 0 0 1 0 0 0 1];...
        [1 0 0 0 0 0 0 0 1]};
masks = cell(4,1);
for i = 1:numel(rows)
    mask = nan(9,9);
    mask(1,:) = rows{i}; %initialise as first row
    for j = 2:size(mask,2) %iterate over first row
        mask(j,:) = mask(1,j) * mask(1,:);
    end
    masks{i} = mask;
end

%% Calculate Averaged Metrology Data R1
for i = 1:numel(C_R1)
    % naive height error wrt bottom of part assuming perfect kinematic coupling
    C_R1(i).data.error = C_R1(i).data.Z - design_data.DesignZ;
    
    % calculate local base error at each contact point
    C_R1(i).coupling_pts.error = C_R1(i).coupling_pts.Z - 5.73; %nominal dist between center of balls and top surface of base
    
    % fit plane to coupling error
    C_R1(i).coupling_lm = fitlm([90 15; 16.3878 142.5; 163.6121 142.5], C_R1(i).coupling_pts.error, 'linear');
    
    % predict coupling error
    C_R1(i).predicted_coupling_error = predict(C_R1(i).coupling_lm, [design_data.DesignX design_data.DesignY]);
    
    % subtract predicted coupling error from measured height error (wrt center of balls)
    C_R1(i).data.error  = C_R1(i).data.error - C_R1(i).predicted_coupling_error;
end

for i = 1:numel(C_R1)
    y1(:,i) = C_R1(i).data.error;
end

y1 = mean(y1,2);

for i = 1:numel(C_R2)
    % naive height error wrt bottom of part assuming perfect kinematic coupling
    C_R2(i).data.error = C_R2(i).data.Z - design_data.DesignZ;
    
    % calculate local base error at each contact point
    C_R2(i).coupling_pts.error = C_R2(i).coupling_pts.Z - 5.73; %nominal dist between center of balls and top surface of base
    
    % fit plane to coupling error
    C_R2(i).coupling_lm = fitlm([90 15; 16.3878 142.5; 163.6121 142.5], C_R2(i).coupling_pts.error, 'linear');
    
    % predict coupling error
    C_R2(i).predicted_coupling_error = predict(C_R2(i).coupling_lm, [design_data.DesignX design_data.DesignY]);
    
    % subtract predicted coupling error from measured height error (wrt center of balls)
    C_R2(i).data.error  = C_R2(i).data.error - C_R2(i).predicted_coupling_error;
end

for i = 1:numel(C_R2)
    y2(:,i) = C_R2(i).data.error;
end

%% Fit models to error
lms = cell(4,1);

for i = 1:4
    lin_mask = reshape(logical(masks{i}),81,1);
    selected_design_data = repmat(design_data(lin_mask,:),2,1);
    selected_z_err = [y1(lin_mask,:); y2(lin_mask,:)];
    predictors = [selected_design_data.DesignX selected_design_data.DesignY selected_design_data.DesignZ];
    lms{i} = fitlm(predictors, selected_z_err ,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);
end

%% Compare estimated coefficients
coeffs = nan(4,4);
for i = 1:4
    coeffs(:,i) = lms{i}.Coefficients.Estimate;
end

%% Nested variance analysis
%[p,anova_tbl] = anovan(error,{design_data.DesignX,design_data.DesignY,design_data.DesignZ},'nested',[0 0 1; 0 0 1; 0 0 0]);