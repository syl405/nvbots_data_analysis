clear

% Script for importing data from the following text file:
%
%    D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP1\parsed_raw_data.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/21 14:00:59

%% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP2\parsed_raw_data.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: categorical (%C)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%C%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
metrology_data = table(dataArray{1:end-1}, 'VariableNames', {'SpecimenCode','FeatureNumber','Z','Flatness','XYParallelism'});

%% Renumber feature numbers to be in sequence
metrology_data = sortrows(metrology_data,{'FeatureNumber'},{'ascend'}); %sort by feature number
%metrology_data.FeatureNumber = (1:9)';

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\artifact_R1 - Variable Pattern Table - VarPattern.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/21 14:04:43

%% Import the data
[~, ~, raw] = xlsread('D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\artifact_R2 - Variable Pattern Table - VarPattern.xlsx','Sheet1','A3:E83');

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
design_data = table;

%% Allocate imported array to column variable names
design_data.DesignX = data(:,3);
design_data.DesignY = 200 - data(:,4);
design_data.DesignZ = data(:,5);
design_data = design_data(data(:,2)==0,:);
design_data.FeatureNumber = (1:9)';

%% Clear temporary variables
clearvars data raw;

%% Create error surface plot
error = metrology_data.Z - design_data.DesignZ;
x = reshape(design_data.DesignX,3,3);
y = reshape(design_data.DesignY,3,3);
z_err = reshape(error,3,3);
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
surf(x,y,reshape(parallelism,3,3))
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Parallelism (mm)')
title('Parallelism')
savefig(gcf, 'figures/parallelism_surf_plot.fig')
close(gcf)

flatness = metrology_data.Flatness;
figure
surf(x,y,reshape(flatness,3,3))
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Flatness (mm)')
title('Flatness')
savefig(gcf, 'figures/flatness_surf_plot.fig')
close(gcf)

%% Fit model to error
lm_lin = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],error,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',false,'Exclude',3); %exclude outlier point (visibly bad column)

%% Save model coefficients to csv file
model_coefficients_path = './error_model_coefficients.csv';
csvwrite(model_coefficients_path, lm_lin.Coefficients.Estimate)

%% Nested variance analysis
%[p,anova_tbl] = anovan(error,{design_data.DesignX,design_data.DesignY,design_data.DesignZ},'nested',[0 0 1; 0 0 1; 0 0 0]);