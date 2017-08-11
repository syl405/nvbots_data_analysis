clear

%% Get metrology data
C_uncomped = mapfix_analysis;
C_comped = mapcomped_analysis_full;

%exclude broken columns
C_comped(1).data(C_comped(1).data.FeatureNumber == 37,:)=[];
C_comped(1).data(C_comped(1).data.FeatureNumber == 46,:)=[];
C_comped(1).data(C_comped(1).data.FeatureNumber == 55,:)=[];
C_comped(2).data(C_comped(2).data.FeatureNumber == 37,:)=[];
C_comped(2).data(C_comped(2).data.FeatureNumber == 46,:)=[];
C_comped(2).data(C_comped(2).data.FeatureNumber == 55,:)=[];
C_uncomped(1).data(C_uncomped(1).data.FeatureNumber == 37,:)=[];
C_uncomped(1).data(C_uncomped(1).data.FeatureNumber == 46,:)=[];
C_uncomped(1).data(C_uncomped(1).data.FeatureNumber == 55,:)=[];
C_uncomped(2).data(C_uncomped(2).data.FeatureNumber == 37,:)=[];
C_uncomped(2).data(C_uncomped(2).data.FeatureNumber == 46,:)=[];
C_uncomped(2).data(C_uncomped(2).data.FeatureNumber == 55,:)=[];

%% Get design data
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

% Exclude broken columns
design_data(design_data.FeatureNumber == 37,:)=[];
design_data(design_data.FeatureNumber == 46,:)=[];
design_data(design_data.FeatureNumber == 55,:)=[];

% Clear temporary variables
clearvars data raw;

%% Calculate averaged metrology data for comped part
%=========================================================
% average measurements for column tops and coupling points
%=========================================================
D_comped_avg = C_comped(1).data;
D_comped_avg.Flatness = mean([C_comped(1).data.Flatness,C_comped(2).data.Flatness],2);
D_comped_avg.Z = mean([C_comped(1).data.Z,C_comped(2).data.Z],2);
D_comped_avg.XYParallelism = mean([C_comped(1).data.XYParallelism,C_comped(2).data.XYParallelism],2);
CP_comped_avg = C_comped(1).coupling_pts;
CP_comped_avg.Flatness = mean([C_comped(1).coupling_pts.Flatness,C_comped(2).coupling_pts.Flatness],2);
CP_comped_avg.Z = mean([C_comped(1).coupling_pts.Z,C_comped(2).coupling_pts.Z],2);
CP_comped_avg.XYParallelism = mean([C_comped(1).coupling_pts.XYParallelism,C_comped(2).coupling_pts.XYParallelism],2);
CP_comped_avg.error = CP_comped_avg.Z - 5.73;

% D_comped_avg = C_comped(1).data;
% D_comped_avg.Flatness = C_comped(1).data.Flatness;
% D_comped_avg.Z = C_comped(1).data.Z;
% D_comped_avg.XYParallelism = C_comped(1).data.XYParallelism;
% CP_comped_avg = C_comped(1).coupling_pts;
% CP_comped_avg.Flatness = C_comped(1).coupling_pts.Flatness;
% CP_comped_avg.Z = C_comped(1).coupling_pts.Z;
% CP_comped_avg.XYParallelism = C_comped(1).coupling_pts.XYParallelism;
% CP_comped_avg.error = CP_comped_avg.Z - 5.73;


%==================
% recalculate error
%==================
% naive height error wrt bottom of part assuming perfect kinematic coupling
D_comped_avg.error = D_comped_avg.Z - design_data.DesignZ;
% fit plane to coupling error
comped_coupling_lm = fitlm([90 15; 16.3878 142.5; 163.6121 142.5], CP_comped_avg.error, 'linear');
% predict coupling error
comped_predicted_coupling_error = predict(comped_coupling_lm, [design_data.DesignX design_data.DesignY]);

% subtract predicted coupling error from measured height error (wrt center of balls)
D_comped_avg.error  = D_comped_avg.error - comped_predicted_coupling_error;

%% Calculate averaged metrology data for uncomped part
%=========================================================
% average measurements for column tops and coupling points
%=========================================================
D_uncomped_avg = C_uncomped(1).data;
D_uncomped_avg.Flatness = mean([C_uncomped(1).data.Flatness,C_uncomped(2).data.Flatness],2);
D_uncomped_avg.Z = mean([C_uncomped(1).data.Z,C_uncomped(2).data.Z],2);
D_uncomped_avg.XYParallelism = mean([C_uncomped(1).data.XYParallelism,C_uncomped(2).data.XYParallelism],2);
CP_uncomped_avg = C_uncomped(1).coupling_pts;
CP_uncomped_avg.Flatness = mean([C_uncomped(1).coupling_pts.Flatness,C_uncomped(2).coupling_pts.Flatness],2);
CP_uncomped_avg.Z = mean([C_uncomped(1).coupling_pts.Z,C_uncomped(2).coupling_pts.Z],2);
CP_uncomped_avg.XYParallelism = mean([C_uncomped(1).coupling_pts.XYParallelism,C_uncomped(2).coupling_pts.XYParallelism],2);
CP_uncomped_avg.error = CP_uncomped_avg.Z - 5.73;

%==================
% recalculate error
%==================
% naive height error wrt bottom of part assuming perfect kinematic coupling
D_uncomped_avg.error = D_uncomped_avg.Z - design_data.DesignZ;
% fit plane to coupling error
uncomped_coupling_lm = fitlm([90 15; 16.3878 142.5; 163.6121 142.5], CP_uncomped_avg.error, 'linear');
% predict coupling error
uncomped_predicted_coupling_error = predict(uncomped_coupling_lm, [design_data.DesignX design_data.DesignY]);
% subtract predicted coupling error from measured height error (wrt center of balls)
D_uncomped_avg.error  = D_uncomped_avg.error - uncomped_predicted_coupling_error;

%% Clear redundant variables
clear C_comped C_uncomped

%% Calculate RMS Height Error
rmse_comped = sqrt(mean(D_comped_avg.error.^2));
rmse_uncomped = sqrt(mean(D_uncomped_avg.error.^2));

%% Calculate change before and after compensation
abs_improvement = abs(D_uncomped_avg.error) - abs(D_comped_avg.error);
rel_improvement = abs_improvement./abs(D_uncomped_avg.error);

%% Plot histograms of error before and after compensation
figure
histogram(D_comped_avg.error,15,'facecolor','green','facealpha',.5,'edgecolor','none')
hold on
histogram(D_uncomped_avg.error,15,'facecolor','red','facealpha',.5,'edgecolor','none')

%% Calculate process capability metrics
S_comped = capability(D_comped_avg.error,[-0.3 0.3]);
S_uncomped = capability(D_uncomped_avg.error,[-0.3 0.3]);

%% Make process capability plots
figure
subplot(1,2,1)
P_comped = capaplot(D_comped_avg.error,[S_comped.mu-2*S_comped.sigma S_comped.mu+2*S_comped.sigma]);
title(['95% Height Error Interval with Compensation: (' num2str(S_comped.mu-2*S_comped.sigma) ',' num2str(S_comped.mu+2*S_comped.sigma) ')'])
hold on
subplot(1,2,2)
P_uncomped = capaplot(D_uncomped_avg.error,[S_uncomped.mu-2*S_uncomped.sigma S_uncomped.mu+2*S_uncomped.sigma]);
title(['95% Height Error Interval w/o Compensation: (' num2str(S_uncomped.mu-2*S_uncomped.sigma) ',' num2str(S_uncomped.mu+2*S_uncomped.sigma) ')'])
