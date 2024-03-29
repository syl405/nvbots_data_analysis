clear

%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\shien\Documents\NVBOTS\data_analysis\PGK Data\PLA\column_experiment_remeasured_v2.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/27 10:57:04

%% Import the data
[~, ~, raw] = xlsread('C:\Users\shien\Documents\NVBOTS\data_analysis\PGK Data\PLA\column_experiment_remeasured_v2.xlsx','Sheet1','A2:AB38');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,[1,2,4,5,6,9,11,14,15,16,17]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[3,7,8,10,12,13,18,19,20,21,22,23,24,25,26,27,28]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
R = table;

%% Allocate imported array to column variable names
R.name = categorical(stringVectors(:,1));
R.experiment = categorical(stringVectors(:,2));
R.number = data(:,1);
R.printer = categorical(stringVectors(:,3));
R.location = categorical(stringVectors(:,4));
R.date = stringVectors(:,5);
R.time = data(:,2);
R.infill = data(:,3);
%columnexperimentremeasuredv2.area = categorical(stringVectors(:,6));
R.Area = data(:,4);
R.material = categorical(stringVectors(:,7));
R.spool = data(:,5);
R.nozzle = data(:,6);
R.color = categorical(stringVectors(:,8));
R.quality = categorical(stringVectors(:,9));
R.structure = categorical(stringVectors(:,10));
R.operator = categorical(stringVectors(:,11));
%columnexperimentremeasuredv2.desired_height = data(:,7);
R.nominal_build_height = data(:,8);
R.layer_number = data(:,9);
R.nominal_layer_thickness = data(:,10);
R.DP = data(:,11);
R.PL = data(:,12);
R.FL_LOW = data(:,13);
R.FL_HIGH = data(:,14);
%columnexperimentremeasuredv2.quality_b = data(:,15);
R.xbar = data(:,16);
R.deltaxbar = data(:,17);

%% Clear temporary variables
clearvars data raw stringVectors;

%% Truncate data to honeycomb samples
D_honeycomb = R(R.structure=='Honeycomb',:);

%% Fit linear models
%fast quality
D_honeycomb_fast = D_honeycomb(D_honeycomb.quality == 'Quick',:);
lm_cubic_fast = fitlm(D_honeycomb_fast.nominal_build_height,D_honeycomb_fast.deltaxbar,'mean_layer_thickness ~ nominal_build_height^3','VarNames',{'nominal_build_height','mean_layer_thickness'});
lm_lin_fast = fitlm(D_honeycomb_fast.nominal_build_height,D_honeycomb_fast.deltaxbar,'linear','VarNames',{'nominal_build_height','mean_layer_thickness'});

%% Predict cumulative errors
probe = (0:0.3:240)';
predicted_lbl_error = predict(lm_cubic_fast,probe);
predicted_cumul_error = nan(numel(probe),1);
for i = 1:numel(predicted_cumul_error)
    predicted_cumul_error(i) = sum(predicted_lbl_error(1:i));
end

%% Fit model to predicted cumulative errors
lm_cumul_quartic = fitlm(probe,predicted_cumul_error,'y~x1^4-1'); %this is an "exact" fit, just a quick way to get coefficients - relationship between layer thickness model and this model is deterministic