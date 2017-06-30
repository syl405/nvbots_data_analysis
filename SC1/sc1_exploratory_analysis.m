clear
%% Import data from text file.
% Script for importing data from the following text file:
%
%    D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\input_screening_experiment\metrology_data\sc1\parsed_raw_data.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/12 15:31:17

%% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\input_screening_experiment\metrology_data\sc1\parsed_raw_data.csv';
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

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
raw_data = table(dataArray{1:end-1}, 'VariableNames', {'specimen_code','datum_flatness','col_1_dist','col_1_flatness','col_1_parallelism','col_2_dist','col_2_flatness','col_2_parallelism','col_3_dist','col_3_flatness','col_3_parallelism','col_4_dist','col_4_flatness','col_4_parallelism'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% Generate DOE
design = fracfact('W L M V WLM LMV');
design(design == -1) = 0; % convert to zeros and ones

%% sort and reshape raw data into separate vectors
datum_flatness = []; %nx1
col_flatness = []; %nx4
mean_layer_thickness = []; %nx4
col_parallelism = []; %nx4
nominal_slice_thickness = []; %nx1
for i=1:size(design,1)
    if design(i,2) == 0
        cur_nominal_slice_thickness = 0.1;
    else
        cur_nominal_slice_thickness = 0.35;
    end
    
    specimen_code = ['1-'...
                     num2str(design(i,1))...
                     num2str(design(i,2))...
                     num2str(design(i,3))...
                     num2str(design(i,4))...
                     num2str(design(i,5))...
                     num2str(design(i,6))...
                     '-1'];
    row = raw_data(strcmp(raw_data.specimen_code,specimen_code),:);
    datum_flatness = [datum_flatness; row.datum_flatness];
    col_flatness = [col_flatness; row.col_1_flatness row.col_2_flatness row.col_3_flatness row.col_4_flatness];
    col_parallelism = [col_parallelism; row.col_1_parallelism row.col_2_parallelism row.col_3_parallelism row.col_4_parallelism];
    mean_layer_thickness = [mean_layer_thickness; row.col_1_dist/4 row.col_2_dist/4 row.col_3_dist/4 row.col_4_dist/4];
    nominal_slice_thickness = [nominal_slice_thickness; cur_nominal_slice_thickness];
end

mean_layer_thickness_deviation = (mean_layer_thickness - nominal_slice_thickness)./nominal_slice_thickness;
mean_col_flatness = mean(col_flatness,2);
sd_col_flatness = std(col_flatness,0,2);
mean_col_parallelism = mean(col_parallelism,2);
sd_col_parallelism = std(col_parallelism,0,2);

%% Fit regression models
var_names = {'extrusion_width','slice_thickness','extrusion_multiplier','print_speed','cooling','temperature'};
% column-level models
lm_indiv_deviation = fitlm([design;design;design;design],reshape(mean_layer_thickness_deviation,[],1),'linear','VarNames',[var_names,'col_level_mean_layer_thickness_deviation']);
lm_indiv_abs = fitlm([design;design;design;design],reshape(mean_layer_thickness,[],1),'linear','VarNames',[var_names,'col_level_mean_layer_thickness']);

% part-level models
lm_part_deviation = fitlm(design,mean(mean_layer_thickness_deviation,2),'linear','VarNames',[var_names,'part_level_mean_layer_thickness_deviation']); %layer thickness deviation, part-level mean
lm_std_deviation = fitlm(design,std(mean_layer_thickness_deviation,0,2),'linear','VarNames',[var_names,'stdev_of_layer_thickness_deviation']); %std of layer thickness deviation within parts

lm_part_abs = fitlm(design,mean(mean_layer_thickness,2),'linear','VarNames',[var_names,'part_level_mean_layer_thickness']); %layer thickness deviation, part-level mean
lm_std_abs = fitlm(design,std(mean_layer_thickness,0,2),'linear','VarNames',[var_names,'stdev_of_layer_thickness']); %std of layer thickness deviation within parts

%% Residual Plots
% Individual column model
fh = figure(); %figure handle
title('Residual Plots - Individual Column Layer Thickness Model'); %overall figure title
subplot(2,2,1); %subplot handle [TODO: make number of subplots dynamic to number of model forms]

plotResiduals(lm_indiv_deviation,'histogram')
subplot(2,2,2)
plotResiduals(lm_indiv_deviation,'probability')
subplot(2,2,3)
plotResiduals(lm_indiv_deviation,'fitted')

savefig(gcf,'exploratory_analysis/indiv_col_model/residual_analysis/Residual Plots.fig')
saveas(gcf,'exploratory_analysis/indiv_col_model/residual_analysis/Residual Plots.eps','epsc')
saveas(gcf,'exploratory_analysis/indiv_col_model/residual_analysis/Residual Plots.png')
close(gcf)

% Part level mean model
fh = figure(); %figure handle
title('Residual Plots - Individual Column Layer Thickness Model'); %overall figure title
subplot(2,2,1); %subplot handle [TODO: make number of subplots dynamic to number of model forms]

plotResiduals(lm_part_deviation,'histogram')
subplot(2,2,2)
plotResiduals(lm_part_deviation,'probability')
subplot(2,2,3)
plotResiduals(lm_part_deviation,'fitted')

savefig(gcf,'exploratory_analysis/part_mean_model/residual_analysis/Residual Plots.fig')
saveas(gcf,'exploratory_analysis/part_mean_model/residual_analysis/Residual Plots.eps','epsc')
saveas(gcf,'exploratory_analysis/part_mean_model/residual_analysis/Residual Plots.png')
close(gcf)

% Part level stdev model
fh = figure(); %figure handle
title('Residual Plots - Individual Column Layer Thickness Model'); %overall figure title
subplot(2,2,1); %subplot handle [TODO: make number of subplots dynamic to number of model forms]

plotResiduals(lm_std_deviation,'histogram')
subplot(2,2,2)
plotResiduals(lm_std_deviation,'probability')
subplot(2,2,3)
plotResiduals(lm_std_deviation,'fitted')

savefig(gcf,'exploratory_analysis/part_std_model/residual_analysis/Residual Plots.fig')
saveas(gcf,'exploratory_analysis/part_std_model/residual_analysis/Residual Plots.eps','epsc')
saveas(gcf,'exploratory_analysis/part_std_model/residual_analysis/Residual Plots.png')
close(gcf)

%% Plot run chart
plot([raw_data.col_1_dist raw_data.col_2_dist raw_data.col_3_dist raw_data.col_4_dist]/4,'x')
axis([0 17 0 0.5])
xlabel('Print sequence')
ylabel('Mean Layer Height (mm)')
legend('Column 1','Column 2','Column 3','Column 4','Location','northwest')
savefig(gcf,'exploratory_analysis/run_chart.fig')
saveas(gcf,'exploratory_analysis/run_chart.eps','epsc')
saveas(gcf,'exploratory_analysis/run_chart.png')
close(gcf)