clear

%% Import Metrology Data from Measurement 1

% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data_auto.csv';
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
metrology_data = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','Flatness','Z','FeatureNumber','XYParallelism'});

% Renumber feature numbers to be in sequence
metrology_data = sortrows(metrology_data,{'FeatureNumber'},{'ascend'}); %sort by feature number
metrology_data.FeatureNumber = (1:81)';

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import Metrology Data from Measurement 2

% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data_auto_verify.csv';
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
metrology_data_ver = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Flatness','Z','XYParallelism'});

% Renumber feature numbers to be in sequence
metrology_data_ver = sortrows(metrology_data_ver,{'FeatureNumber'},{'ascend'}); %sort by feature number
metrology_data_ver.FeatureNumber = (1:81)';

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import Metrology Data from Measurement 3
% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data_auto_verify_2.csv';
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
metrology_data_ver2 = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Flatness','Z','XYParallelism'});

% Renumber feature numbers to be in sequence
metrology_data_ver2 = sortrows(metrology_data_ver2,{'FeatureNumber'},{'ascend'}); %sort by feature number
metrology_data_ver2.FeatureNumber = (1:81)';

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import Metrology Data from Measurement 4
% Script for importing data from the following text file:
%
%    D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/21 14:00:59

% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data_auto_verify_3.csv';
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

% Create output variable
metrology_data_ver3 = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Flatness','Z','XYParallelism'});

% Renumber feature numbers to be in sequence
metrology_data_ver3 = sortrows(metrology_data_ver3,{'FeatureNumber'},{'ascend'}); %sort by feature number
metrology_data_ver3.FeatureNumber = (1:81)';

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import Part Design Data

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

% Clear temporary variables
clearvars data raw;

%% Create error, flatness, and parallelism surface plots
error = metrology_data.Z - design_data.DesignZ;
x = reshape(design_data.DesignX,9,9);
y = reshape(design_data.DesignY,9,9);
z_err = reshape(error,9,9);
figure
surf(x,y,z_err)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Height Error (mm)')
title('Height Error')
savefig(gcf, 'figures/error_surf_plot.fig')
close(gcf)

parallelism = metrology_data.XYParallelism;
figure
surf(x,y,reshape(parallelism,9,9))
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Parallelism (mm)')
title('Parallelism')
savefig(gcf, 'figures/parallelism_surf_plot.fig')
close(gcf)

flatness = metrology_data.Flatness;
figure
surf(x,y,reshape(flatness,9,9))
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Flatness (mm)')
title('Flatness')
savefig(gcf, 'figures/flatness_surf_plot.fig')
close(gcf)

%% Fit XYZ error models
%calculate measured height errors
errors = [error metrology_data_ver.Z-design_data.DesignZ metrology_data_ver2.Z-design_data.DesignZ metrology_data_ver3.Z-design_data.DesignZ];
    
lm_lin = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],errors(:,1),'interactions', 'VarNames', {'X','Y','Z','Error'},'Intercept',false);
lm_lin_ver = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],errors(:,2),'interactions', 'VarNames', {'X','Y','Z','Error'},'Intercept',false);
lm_lin_ver2 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],errors(:,3),'interactions', 'VarNames', {'X','Y','Z','Error'},'Intercept',false);
lm_lin_ver3 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],errors(:,4),'interactions', 'VarNames', {'X','Y','Z','Error'},'Intercept',false);

%% Export model coefficients to csv file for compensation use
model_coefficients_path = './error_model_coefficients.csv';
csvwrite(model_coefficients_path, lm_lin.Coefficients.Estimate)

%% Analyze measurement uncertainty
repeat_error = metrology_data.Z - metrology_data_ver.Z;
repeat_error_2 = metrology_data.Z - metrology_data_ver2.Z;
repeat_error_3 = metrology_data.Z - metrology_data_ver3.Z;

%=======================================================================
% Check whether there is spatial dependence in measurement repeatability
%========================================================================
lm_rep = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);
lm_rep_2 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error_2,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);
lm_rep_3 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error_3,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);

%% Formal Gauge R&R Study
%comparing individual columns
grr_data = [metrology_data.Z metrology_data_ver.Z metrology_data_ver2.Z metrology_data_ver3.Z];
x_bar = mean(grr_data,2);
s = std(grr_data,0,2);
r = range(grr_data,2);
sigma_gauge_all = mean(r)/2.059; %d2 factor from montgomery
sigma_total_all = std(reshape(grr_data,[],1));
sigma_part_all = sqrt(sigma_total_all^2 - sigma_gauge_all^2);
rho_p_all = sigma_part_all^2/sigma_total_all^2;
snr_all = sqrt(2*rho_p_all/(1-rho_p_all));

%comparing columns of same nominal height
heights = unique(design_data.DesignZ);
grr_data_by_height = nan(9,4,numel(heights));
for i = 1:numel(heights)
    grr_data_by_height(:,:,i) = grr_data(design_data.DesignZ==heights(i),:);
end
r_by_height = range(grr_data_by_height,2);
sigma_gauge_by_height = squeeze(mean(r_by_height,1)/2.059);
sigma_total_by_height = nan(size(grr_data_by_height,3),1);
for i = 1:size(grr_data_by_height,3)
    sigma_total_by_height(i) = std(reshape(grr_data_by_height(:,:,i),[],1));
end
sigma_part_by_height = sqrt(sigma_total_by_height.^2 - sigma_gauge_by_height.^2);
rho_p_by_height = sigma_part_by_height.^2./sigma_total_by_height.^2;
snr_by_height = (2*rho_p_by_height./(1-rho_p_by_height)).^0.5;

%% Quantifying Iatrogenics of Model Coefficients Uncertainty
%============================================
%generate grid of points across work envelope
%============================================
[probe_x,probe_y,probe_z] = meshgrid(0:20:180,0:20:200,0:20:240);

%reshape into column vectors
probe_x = reshape(probe_x,[],1);
probe_y = reshape(probe_y,[],1);
probe_z = reshape(probe_z,[],1);
probe_xyz = [probe_x,probe_y,probe_z]; %nx3 matrix

%generate predicted error vectors with various error models
predicted_offset = [predict(lm_lin,probe_xyz) predict(lm_lin_ver,probe_xyz) predict(lm_lin_ver2,probe_xyz) predict(lm_lin_ver3,probe_xyz)];

%calculate dispersion in predicted offset
offset_range = range(predicted_offset,2);
offset_std = std(predicted_offset,0,2);

%plot absolute range of predicted offsets across work envelope
colormap('jet')
scatter3(probe_x,probe_y,probe_z,200,offset_range,'.')
colorbar()
savefig(gcf, 'figures/pred_offset_range_work_envelope.fig')
close(gcf)

%==================================
%probe models at measured positions
%==================================
%generate probe points
probe_arti = [design_data.DesignX design_data.DesignY design_data.DesignZ];

%predict offsets at measured positions
predicted_arti_offset = [predict(lm_lin,probe_arti) predict(lm_lin_ver,probe_arti) predict(lm_lin_ver2,probe_arti) predict(lm_lin_ver3,probe_arti)];

%calculate range of predicted offsets at measured positions
arti_offset_range = range(predicted_arti_offset,2);

%NEED TO CALCULATE RESIDUAL DEVIATIONS AFTER APPLYING PREDICTED OFFSET AS
%COMPENSATION

%express predicted offset range as fraction of actual measured deviation,
%averaged across measurements
arti_rel_offset_range = arti_offset_range./mean(errors,2);

colormap('jet')
scatter3(design_data.DesignX,design_data.DesignY,design_data.DesignZ,500,arti_rel_offset_range,'.')
colorbar()
savefig(gcf, 'figures/rel_pred_offset_range_measured_pts.fig')
close(gcf)
