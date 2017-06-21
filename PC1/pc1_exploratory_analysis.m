clear

%% Parse in raw uncompensated data from PAE1
raw_filename = 'C:/Users/shien/Documents/NVBOTS/piecewise_constant_compensation/error_data/printer_accuracy_arjun.xlsx';
[~,~,headers] = xlsread(raw_filename, 'Z run charts', 'A1:D1','basic'); 
headers = replace(headers,' ','_'); %strip whitespace
headers = replace(headers,'%','Pct'); %strip special characters
[~,~,R_uncomped] = xlsread(raw_filename, 'Z run charts', 'A2:D176','basic');
R_uncomped = cell2table(R_uncomped,'VariableNames',headers);
R_uncomped.Measured_Length = R_uncomped.Nominal_Length + R_uncomped.Raw_Deviation;

% correct for quantization error
D_uncomped = table(changem(R_uncomped.Nominal_Length,[30-20.1 60-20.1 99.9-20.1 180-20.1 240-20.1],[10 40 80 160 220]),'VariableNames',{'Nominal_Length'});
D_uncomped.Measured_Length = R_uncomped.Measured_Length;
D_uncomped.Absolute_Error = D_uncomped.Measured_Length - D_uncomped.Nominal_Length;

%% Calculate Summary Statistics for Uncomped Data
% S = table(unique(D.Nominal_Length+20.1),'VariableNames',{'Nominal_Length'},'RowNames',cellstr(num2str(unique(D.Nominal_Length))));
% nom_length_group = findgroups(D.Nominal_Length);
% S.Mean_Measured_Height = splitapply(@mean,D.Measured_Height,nom_length_group);
S_uncomped = table([20.1;30;60;99.9;180;240],'VariableNames',{'Target_Height'});
S_uncomped.Mean_Measured_Height(S_uncomped.Target_Height==20.1)= 20.1;
for i = 2:numel(S_uncomped.Target_Height)
    S_uncomped.Mean_Measured_Height(i) = 20.1 + mean(D_uncomped.Measured_Length(D_uncomped.Nominal_Length+20.1==S_uncomped.Target_Height(i)));
end

% Calculate mean absolute error
S_uncomped.Mean_Abs_Err = S_uncomped.Mean_Measured_Height-S_uncomped.Target_Height;

% Calculate total error in each block
S_uncomped.Error_This_Block(1) = 0;
S_uncomped.Num_Layers_This_Block(1) = 20.1/0.3;
for i = 2:size(S_uncomped,1)
    S_uncomped.Error_This_Block(i) = S_uncomped.Mean_Abs_Err(i)-S_uncomped.Mean_Abs_Err(i-1);
    S_uncomped.Num_Layers_This_Block(i) = (S_uncomped.Target_Height(i)-S_uncomped.Target_Height(i-1))/0.3;
end
S_uncomped.Error_Per_Layer_This_Block = S_uncomped.Error_This_Block./S_uncomped.Num_Layers_This_Block;

%% Import comped CMM data from text file.
% Script for importing data from the following text file:
%
%    D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\piecewise_compensation_experiment\metrology_data\pc1\parsed_raw_data.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/06/19 11:38:44

%% Initialize variables.
filename = 'D:\Dropbox (MIT)\Spring 2017\NVBOTS Dropbox\Summer\piecewise_compensation_experiment\metrology_data\pc1\parsed_raw_data.csv';
delimiter = ',';

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,3]);
rawStringColumns = string(raw(:, 1));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
R_comped = table;
R_comped.specimen_code = rawStringColumns(2:end, 1);
R_comped.datum_flatness = cell2mat(rawNumericColumns(2:end, 1));
R_comped.perp_dist = cell2mat(rawNumericColumns(2:end, 2));

%% Calculate nominal lengths
for i = 1:size(R_comped,1)
    split_array = R_comped.specimen_code(i).split('-');
    top_target = str2num(split_array{2});
    bot_target = str2num(split_array{3});
    %apply discretization compensation
    if mod(top_target,0.3) > 0.15
        top_target = top_target + (0.3 - mod(top_target,0.3));
    elseif mod(top_target,0.3) < 0.15
        top_target = top_target - mod(top_target,0.3);
    elseif mod(top_target,0.3) == 0
        top_target = top_target;
    end
    if mod(bot_target,0.3) > 0.15
        bot_target = bot_target + (0.3 - mod(bot_target,0.3));
    elseif mod(bot_target,0.3) < 0.15
        bot_target = bot_target - mod(bot_target,0.3);
    elseif mod(bot_target,0.3) == 0
        bot_target = bot_target;
    end
    R_comped.top_target(i) = top_target;
    R_comped.bot_target(i) = bot_target;
end

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;

%% Generate clean data table for comped measurements
D_comped = table(R_comped.top_target-R_comped.bot_target,'VariableNames',{'Nominal_Length'});
D_comped.Measured_Length = R_comped.perp_dist;
D_comped.Absolute_Error = D_comped.Measured_Length - D_comped.Nominal_Length;

%% Pairwise T-tests
grouped_data = cell(1,3); % columns: nom_length, uncomped_error_vector, comped_error_vector
comped_nom_lengths = unique(D_comped.Nominal_Length);
for i = 1:numel(comped_nom_lengths)
    nom_length = comped_nom_lengths(i);
    grouped_data{i,1} = nom_length;
    grouped_data{i,2} = D_uncomped.Absolute_Error(abs(D_uncomped.Nominal_Length-nom_length)<0.00001);
    grouped_data{i,3} = D_comped.Absolute_Error(abs(D_comped.Nominal_Length-nom_length)<0.00001);
end

t_test_output = cell(1,4); % columns: nom_length, uncomped_vs_zero, comped_vs_zero, comped_vs_uncomped
for i = 1:size(grouped_data,1)
    t_test_output{i,1} = grouped_data{i,1};
    [t_test_output{i,2}(1),  t_test_output{i,2}(2)] = ttest(grouped_data{i,2});
    [t_test_output{i,3}(1),  t_test_output{i,3}(2)] = ttest(grouped_data{i,3});
    [t_test_output{i,4}(1),  t_test_output{i,4}(2)] = ttest2(grouped_data{i,2},grouped_data{i,3});
end

%% Mean Abs Error Magnitude Comparisons
mean_abs_error_output = nan(1,4); % columns: nom_length, mean_abs_error_uncomped, mean_abs_error_comped, pct_improvement
for i = 1:size(grouped_data,1)
    mean_abs_error_output(i,1) = grouped_data{i,1};
    mean_abs_error_output(i,2) = mean(grouped_data{i,2});
    mean_abs_error_output(i,3) = mean(grouped_data{i,3});
    mean_abs_error_output(i,4) = (abs(mean_abs_error_output(i,2)) - abs(mean_abs_error_output(i,3)))/abs(mean_abs_error_output(i,2));
end

%% Calculate Summary Statistics for Comped Data
% S = table(unique(D.Nominal_Length+20.1),'VariableNames',{'Nominal_Length'},'RowNames',cellstr(num2str(unique(D.Nominal_Length))));
% nom_length_group = findgroups(D.Nominal_Length);
% S.Mean_Measured_Height = splitapply(@mean,D.Measured_Height,nom_length_group);
S_comped = table([20.1;30;99.9;240],'VariableNames',{'Target_Height'});
S_comped.Mean_Measured_Height(S_comped.Target_Height==20.1)= 20.1;
for i = 2:numel(S_comped.Target_Height)
    S_comped.Mean_Measured_Height(i) = 20.1 + mean(D_comped.Measured_Length(abs(D_comped.Nominal_Length+20.1-S_comped.Target_Height(i))<0.00001));
end

% Calculate mean absolute error
S_comped.Mean_Abs_Err = S_comped.Mean_Measured_Height-S_comped.Target_Height;

% Calculate total error in each block
S_comped.Error_This_Block(1) = 0;
S_comped.Num_Layers_This_Block(1) = 20.1/0.3;
for i = 2:size(S_comped,1)
    S_comped.Error_This_Block(i) = S_comped.Mean_Abs_Err(i)-S_comped.Mean_Abs_Err(i-1);
    S_comped.Num_Layers_This_Block(i) = (S_comped.Target_Height(i)-S_comped.Target_Height(i-1))/0.3;
end
S_comped.Error_Per_Layer_This_Block = S_comped.Error_This_Block./S_comped.Num_Layers_This_Block;