function [C] = mapfix_analysis()

%% Import Metrology Data
% list directory contents
csv_dir_path = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAP_ST\parsed_csv\';
dir_struct = dir(csv_dir_path);

% initialize data structure for imported CMM data
C = struct('data',{},'run',{}, 'coupling_pts', {});

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
    C = [C; struct('data', metrology_data, 'coupling_pts', coupling_point_data, 'run', run)];
end

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

% Apply offset to height data to reference against base of part
design_data.DesignZ = design_data.DesignZ + 3.9;

% Clear temporary variables
clearvars data raw;

%% Import Experimental Run Data
% Import data from spreadsheet
%========================

% Script for importing data from the following spreadsheet:
%
%    Workbook: D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAPFIX1\run_info\expt_log.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/07/12 10:29:27

% Import the data, extracting spreadsheet dates in Excel serial date format
%========================

[~, ~, raw, dates] = xlsread('D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\error_mapping_artifact\metrology_data\MAPFIX1\run_info\expt_log.xlsx','Sheet1','A2:I21','',@convertSpreadsheetExcelDates);
stringVectors = string(raw(:,[1,2,3,5,6,8,9]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,7);
dates = dates(:,4);

% Create output variable
%========================

data = reshape([raw{:}],size(raw));

% Create table
%========================

exptlog = table;

% Allocate imported array to column variable names
%========================

exptlog.Run = stringVectors(:,1);
exptlog.Part = categorical(stringVectors(:,2));
exptlog.Printer = categorical(stringVectors(:,3));
exptlog.Date = categorical(datetime([dates{:,1}].', 'ConvertFrom', 'Excel'));
exptlog.Probe = categorical(stringVectors(:,4));
exptlog.Pallet = categorical(stringVectors(:,5));
exptlog.Part = data(:,1);
exptlog.Done = categorical(stringVectors(:,6));
exptlog.Remarks = categorical(stringVectors(:,7));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% exptlog.Date=datenum(exptlog.Date);

% Clear temporary variables
%========================
clearvars data raw dates stringVectors;

%% Calculate Height Errors

for i = 1:numel(C)
    % naive height error wrt bottom of part assuming perfect kinematic coupling
    C(i).data.error = C(i).data.Z - design_data.DesignZ;
    
    % calculate local base error at each contact point
    C(i).coupling_pts.error = C(i).coupling_pts.Z - 5.73; %nominal dist between center of balls and top surface of base
    
    % fit plane to coupling error
    C(i).coupling_lm = fitlm([90 15; 16.3878 142.5; 163.6121 142.5], C(i).coupling_pts.error, 'linear');
    
    % predict coupling error
    C(i).predicted_coupling_error = predict(C(i).coupling_lm, [design_data.DesignX design_data.DesignY]);
    
    % subtract predicted coupling error from measured height error (wrt center of balls)
    C(i).data.error  = C(i).data.error - C(i).predicted_coupling_error;
end

%% Nested Variance Analysis
% Nesting levels: probe setup > pallet setup > part setup > column number

%reshape data into vector form for MATLAB built-in anovan()
y = nan(81,numel(C)); %matrix containing all height measurements

for i = 1:numel(C)
    y(:,i) = C(i).data.error;
end

y = reshape(y,[],1); %reshape into column vector

%generate grouping vectors
column = repmat(categorical((1:81)'),numel(C),1);
probe =[];
pallet = [];
part = [];
date = [];

for i = 1:numel(C)
    cur_run = C(i).run;
    cur_probe = exptlog(exptlog.Part == cur_run,:).Probe;
    probe = [probe; repmat(cur_probe,81,1)];
    cur_pallet = exptlog(exptlog.Part == cur_run,:).Pallet;
    pallet = [pallet; repmat(cur_pallet,81,1)];
    cur_part = categorical(exptlog(exptlog.Part == cur_run,:).Part);
    part = [part; repmat(cur_part,81,1)];
    cur_date = exptlog(exptlog.Part == cur_run,:).Date;
    date = [date; repmat(cur_date,81,1)];
    clearvars cur_run cur_probe cur_pallet cur_part cur_date
end

%with measurement date
% [P,tbl,stats] = anovan(y,{date, probe, pallet, part, column},...
%                 'nested',[0 0 0 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0],...
%                 'random', [1 2 3 4],...
%                 'varnames',{'Date' 'Probe' 'Pallet' 'Part' 'Column'});
% 
% %without measurement date
% [P,tbl,stats] = anovan(y,{probe, pallet, part, column},...
%             'nested',[0 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0],...
%             'random', [1 2 3],...
%             'varnames',{'Probe' 'Pallet' 'Part' 'Column'});
            
%lumping all measurement-related factors
 [P,tbl,stats] = anovan(y,{part, column},...
                 'random', [1],...
                 'varnames',{'PartSetup' 'Column'},...
                 'display',false);
            
% calculate percentage variance components
nested_var_comps = table(stats.rtnames,'VariableNames',{'Source'});
nested_var_comps.VarComp = stats.varest;
if ~isempty(nested_var_comps(nested_var_comps.VarComp < 0,:))
    nested_var_comps(nested_var_comps.VarComp < 0,:).VarComp = 0;
end
total_var = sum(nested_var_comps.VarComp);
nested_var_comps.VarPercentage = 100 * nested_var_comps.VarComp / total_var;

% calculate signal to noise ratio
rho_part = nested_var_comps.VarPercentage(strcmp(nested_var_comps.Source,'Error'))/100;
gage_snr = sqrt((2*rho_part)/(1-rho_part));

%% Create error, flatness, and parallelism surface plots
% error = metrology_data.Z - design_data.DesignZ;
% x = reshape(design_data.DesignX,9,9);
% y = reshape(design_data.DesignY,9,9);
% z_err = reshape(error,9,9);
% figure
% surf(x,y,z_err)
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Height Error (mm)')
% title('Height Error')
% savefig(gcf, 'figures/error_surf_plot.fig')
% close(gcf)
% 
% parallelism = metrology_data.XYParallelism;
% figure
% surf(x,y,reshape(parallelism,9,9))
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Parallelism (mm)')
% title('Parallelism')
% savefig(gcf, 'figures/parallelism_surf_plot.fig')
% close(gcf)
% 
% flatness = metrology_data.Flatness;
% figure
% surf(x,y,reshape(flatness,9,9))
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Flatness (mm)')
% title('Flatness')
% savefig(gcf, 'figures/flatness_surf_plot.fig')
% close(gcf)

%% Create XYZ Height Table
% x = reshape(design_data.DesignX,9,9);
% y = reshape(design_data.DesignY,9,9);
% z = reshape(design_data.DesignZ,9,9);
% 
% input.data = [[NaN; y(1,:)'] [x(:,1)'; z]];
% input.tableBorders = false;
% input.dataFormat = {'%.1f'};
% latexTable(input)

%% Fit XYZ error models
% fit model for each measurement (for iatrogenics analysis)
for i = 1:numel(C)
    C(i).error_lm = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],C(i).data.error, 'linear', 'VarNames', {'X','Y','Z','HeightError'});
end

% fit model using repeated measurements as replicates ("Error" componenent refers to measurement error + inherent randomness.
predictor_mat = repmat([design_data.DesignX design_data.DesignY design_data.DesignZ],numel(C),1);
lm_repeat = fitlm(predictor_mat, y, 'linear', 'VarNames', {'X','Y','Z','HeightError'});

% fit model using average measurements as replicates ("Error" componenent refers to inherent randomness only
predictor_mat_avg =[design_data.DesignX design_data.DesignY design_data.DesignZ];
lm_avg = fitlm(predictor_mat_avg, mean(reshape(y,size(predictor_mat_avg,1),[]),2), 'poly333', 'VarNames', {'X','Y','Z','HeightError'},'intercept',true);

%remove insignificant terms
insig_terms = '';
for i = 1:size(lm_avg.Coefficients,1)
    if lm_avg.Coefficients.pValue(i) > 0.001 % current term insignificant
        if isempty(insig_terms)
            if strcmp(lm_avg.CoefficientNames{i},'(Intercept)')
                insig_terms = '1';
            else
                insig_terms = lm_avg.CoefficientNames{i};
            end
        else
            insig_terms = [insig_terms ' + ' lm_avg.CoefficientNames{i}];
        end
    end
end
lm_avg = removeTerms(lm_avg,insig_terms); 

%% Plot predicted height offsets using aggregate model ("true error independent of measurement error")
[probe_x,probe_y,probe_z] = meshgrid(0:20:180,0:20:200,0:20:100);

%reshape into column vectors
probe_x_linear = reshape(probe_x,[],1);
probe_y_linear = reshape(probe_y,[],1);
probe_z_linear = reshape(probe_z,[],1);
probe_xyz = [probe_x_linear,probe_y_linear,probe_z_linear]; %nx3 matrix

%generate predicted error vectors with various error models
predicted_offset = predict(lm_avg,probe_xyz);

%height error predicted by aggregate model
%plot absolute range of predicted offsets across work envelope
colormap('jet')
scatter3(probe_x_linear,probe_y_linear,probe_z_linear,200,predicted_offset,'.')
axis([0 180 0 200 0 240])
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
c = colorbar();
ylabel(c,'Predicted Height Offset(mm)')
title('Predicted Height Offset From Repeated Measurements')
savefig(gcf,'figures/pred_total_height_error_repeated_measurements.fig')
close(gcf)

%generate predicted length offsets (step height) between top and 20.1 mm ref base (for arjun comparison)
[probe_x_ref,probe_y_ref,probe_z_ref] = meshgrid(0:20:180,0:20:200,20.1);
probe_x_ref = reshape(probe_x_ref,[],1);
probe_y_ref = reshape(probe_y_ref,[],1);
probe_z_ref = reshape(probe_z_ref,[],1);
probe_xyz_ref = [probe_x_ref,probe_y_ref,probe_z_ref]; %nx3 matrix

predicted_offset_ref_base_arjun = predict(lm_repeat, probe_xyz_ref);
predicted_step_length_error = predicted_offset - repmat(predicted_offset_ref_base_arjun,numel(0:20:100),1);

%height error predicted by aggregate model
%plot absolute range of predicted offsets across work envelope
colormap('jet')
scatter3(probe_x_linear,probe_y_linear,probe_z_linear,200,predicted_step_length_error,'.')
axis([0 180 0 200 0 240])
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
c = colorbar();
ylabel(c,'Predicted Step Length Offset(mm)')
title('Predicted Step Length Offset From Repeated Measurements')
savefig(gcf,'figures/pred_step_length_error_repeated_measurements.fig')
close(gcf)

%% Export model coefficients to csv file for compensation use
model_coefficients_path = './error_model_coefficients.csv';
predictor_names = {'(Intercept)','X','Y','Z','X^2','X:Y','Y^2','X:Z','Y:Z','Z^2','X^3','X^2:Y','X:Y^2','Y^3','X^2:Z','X:Y:Z','Y^2:Z','X:Z^2','Y:Z^2','Z^3'};
model_coefficients = zeros(numel(predictor_names),1); %initialize coefficient vector
for i = 1:numel(predictor_names)
    if lm_avg.Coefficients.Estimate(find(strcmp(lm_avg.CoefficientNames,predictor_names{i})))
        model_coefficients(i) = lm_avg.Coefficients.Estimate(find(strcmp(lm_avg.CoefficientNames,predictor_names{i})));
    end
end
dlmwrite(model_coefficients_path, model_coefficients, 'delimiter', ',', 'precision', 12); 

%% Analyze measurement uncertainty
%repeat_error = metrology_data.Z - metrology_data_ver.Z;
% repeat_error_2 = metrology_data.Z - metrology_data_ver2.Z;
% repeat_error_3 = metrology_data.Z - metrology_data_ver3.Z;

%=======================================================================
% Check whether there is spatial dependence in measurement repeatability
%========================================================================
%lm_rep = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);
% lm_rep_2 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error_2,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);
% lm_rep_3 = fitlm([design_data.DesignX design_data.DesignY design_data.DesignZ],repeat_error_3,'linear', 'VarNames', {'X','Y','Z','Error'},'Intercept',true);

%% Formal Gauge R&R Study
%comparing individual columns
err = nan(81,numel(C)); %matrix containing all height measurements
for i = 1:numel(C)
    err(:,i) = C(i).data.error;
end
err = reshape(err,[],1); %reshape into column vector

grr_data = reshape(err,81,[]);
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
grr_data_by_height = nan(9,size(grr_data,2),numel(heights));
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
[probe_x,probe_y,probe_z] = meshgrid(0:20:180,0:20:200,0:20:100);

%reshape into column vectors
probe_x = reshape(probe_x,[],1);
probe_y = reshape(probe_y,[],1);
probe_z = reshape(probe_z,[],1);
probe_xyz = [probe_x,probe_y,probe_z]; %nx3 matrix

%generate predicted error vectors with various error models
predicted_offset = nan(size(probe_xyz,1),numel(C));
for i = 1:numel(C)
    predicted_offset(:,i) = predict(C(i).error_lm,probe_xyz);
end
%calculate dispersion in predicted offset
offset_range = range(predicted_offset,2);
offset_std = std(predicted_offset,0,2);

%plot absolute range of predicted offsets across work envelope
colormap('jet')
scatter3(probe_x,probe_y,probe_z,200,offset_range,'.')
axis([0 180 0 200 0 240])
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
c = colorbar();
ylabel(c,'Range of Predicted Height Offset (mm)')
title(['Range of Predicted Height Offset Across' num2str(numel(C)) ' Measurements'])
savefig(gcf, 'figures/pred_offset_range_work_envelope.fig')
close(gcf)

%==================================
%probe models at measured positions
%==================================
%generate probe points
probe_arti = [design_data.DesignX design_data.DesignY design_data.DesignZ];

%predict offsets at measured positions
predicted_arti_offset = nan(size(probe_arti,1),numel(C));
for i = 1:numel(C)
    predicted_arti_offset(:,i) = predict(C(i).error_lm,probe_arti);
end

%calculate range of predicted offsets at measured positions
arti_offset_range = range(predicted_arti_offset,2);

%NEED TO CALCULATE RESIDUAL DEVIATIONS AFTER APPLYING PREDICTED OFFSET AS
%COMPENSATION

%express predicted offset range as fraction of actual measured deviation,
%averaged across measurements
errors = nan(81,numel(C));
for i = 1:numel(C)
    errors(:,i) = C(i).data.error;
end

arti_rel_offset_range = arti_offset_range./abs(mean(errors,2));

colormap('jet')
scatter3(design_data.DesignX,design_data.DesignY,design_data.DesignZ,500,arti_rel_offset_range,'.')
c=colorbar();
caxis([0 5])
savefig(gcf, 'figures/rel_pred_offset_range_measured_pts_mean.fig')
close(gcf)

%express predicted offset range as fraction of actual measured deviation,
%in each measurements
arti_rel_offset_range_indiv = repmat(arti_offset_range,1,numel(C))./abs(errors);

figure
for i = 1:numel(C)
    subplot(ceil(numel(C)/2),2,i)
    colormap('jet')
    scatter3(design_data.DesignX,design_data.DesignY,design_data.DesignZ,300,arti_rel_offset_range_indiv(:,i),'.')
    colorbar()
    caxis([0 5])
    title(i)
end
savefig(gcf, 'figures/rel_pred_offset_range_measured_pts_indiv.fig')
close(gcf)
%===============================================
% Probe averaged model for prediction intervals
%===============================================

[h_pred,ci_pred] = predict(lm_avg,probe_xyz,'Prediction','curve','Simultaneous',false);
ci_width = range(ci_pred,2);

%% Inspect prediction interval for future observations of height errors
%===============================================
% Probe averaged model for prediction intervals
%===============================================
% using work envelope xyz probe from previous section
[h_pred,ci_pred] = predict(lm_avg,probe_xyz,'Prediction','observation','Simultaneous',false);

%==========================================================
% Check whether existing observations fall in pred interval
%==========================================================
%data_in_int = and((lm_avg.Variables.HeightError > ci_pred(:,1)),(lm_avg.Variables.HeightError < ci_pred(:,2)));

%===============================================
% Summary statistics for prediction interval width
%===============================================
pred_int.data = ci_pred;
pred_int.widths = range(ci_pred,2);
pred_int.mean_width = mean(pred_int.widths);
pred_int.median_width = median(pred_int.widths);
pred_int.min_width = min(pred_int.widths);
pred_int.max_width = max(pred_int.widths);

%% Create morphed pointcloud for illustration
[pc_x,pc_y,pc_z] = meshgrid(0:45:180,0:50:200,0:25:100);

%reshape into column vectors
pc_x_linear = reshape(pc_x,[],1);
pc_y_linear = reshape(pc_y,[],1);
pc_z_linear = reshape(pc_z,[],1);
pc_xyz = [pc_x_linear,pc_y_linear,pc_z_linear]; %nx3 matrix

predicted_offset_pc = predict(lm_avg,pc_xyz);

%apply compensation
pc_z_comped = pc_z_linear - predicted_offset_pc*10;

%plot pointcloud
figure
hold on
hs1 = scatter3(pc_x_linear,pc_y_linear,pc_z_linear,300,'r.');
hs2 = scatter3(pc_x_linear,pc_y_linear,pc_z_comped,300,'g.');
axis([0 180 0 200 0 120])
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
title('Predicted Height Offset From Repeated Measurements')
close(gcf)
