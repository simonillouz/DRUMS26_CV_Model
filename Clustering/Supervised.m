function Supervised(file_name, sheet_name, Norm_Flag)
close all;
% This function conducts supervised machine learning algorithms with
% k-nearest neighbors. 
% Treat EF as known data and parameter values as features

% Input the name of your file and specify the sheet you want data read from
% Input as a string ie, 'Patient Optimized Parameters','Sheet1'

% This function applies supervised machine learning to build a k-Nearest
% Neighbor classifier

% Norm_Flag determines what type of data normalization is conducted 
% 0 for no normalization, 1 for normalization by the mean, and 2 for normalization by the standard deviation

% Goes back into directory to find spreadsheet
cd ../Core/

data_table = readtable(file_name, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve', 'NumHeaderLines', 1, 'ReadVariableNames', true);

% Extract list of patient IDs from 1st row excel header
PatIDs = fillmissing(string(readcell(file_name, 'Range', '1:1')), 'previous');
patient_strings = unique(PatIDs(contains(PatIDs, 'Patient')), 'stable');
PIDs = str2double(regexp(patient_strings, '\d+', 'match', 'once'));

patientdata_table = readtable('AllPatsMRI_2_REU.xlsx','Sheet','Sheet1','VariableNamingRule','preserve');
[~, PID_indices] = ismember(PIDs, patientdata_table.Pat_Num);
EF = patientdata_table.EF(PID_indices);
cd ../Clustering/

% Preparing data for normalization - taking from Cluster.m code
mean_column_flags = contains(data_table.Properties.VariableNames, 'Mean');
mean_matrix_raw = table2array(data_table(:,mean_column_flags));
A_Optp = mean_matrix_raw';
[NumPats_Optp, Num_Optp] = size(A_Optp);

% Taking the mean of each optimized parameter and centering the matrix
AMean_Optp = mean(A_Optp,1);                                % Measure mean across all patients
ACent_Optp = zeros(NumPats_Optp, Num_Optp);                 % Preallocate
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ACent_Optp(i,j) = A_Optp(i,j) - AMean_Optp(j);      % Centering the data matrix    
        end
    end

% Now normalizing (or not) depending on normalization flag
% Norm_Flag = 0 - no normalization 
%           = 1 - normalize by mean
%           = 2 - normalize by standard deviation
AStd_Optp = std(A_Optp);                                    % Measure std dev across all pats
ANorm_Optp = zeros(NumPats_Optp, Num_Optp);                 % Preallocate
if (Norm_Flag == 0)                                         % No normailzation
    ANorm_Optp = ACent_Optp;
elseif (Norm_Flag == 1)                                     % Normalizing by the mean
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ANorm_Optp(i,j) = ...
                ACent_Optp(i,j) / AMean_Optp(j);
        end
    end
elseif (Norm_Flag == 2)                                                        % Normalizing by standard deviation
    for i = 1:NumPats_Optp
        for j = 1:Num_Optp
            ANorm_Optp(i,j) = ...
                ACent_Optp(i,j) / AStd_Optp(j);
        end
    end
end

% Features = parameter values from model
ANorm_Optp = ANorm_Optp(:, 1:end-2); % removes last two columns since all zero

% HF phenotypes:
%   HFpEF: EF > 50% (0.5)
%   HFrEF: EF < 50% (0.5)
edges = [0 .50 1.00];
labels = {'HFrEF', 'HFpEF'};

EF_class = discretize(EF, edges, labels); % separates EF data into HFpEF/HFrEF classes
EF_class = categorical(EF_class); % ensures that data is categorical

%% K-NEAREST NEIGHBORS CLASSIFICATION
% Split into training and test sets
rng(1); % for reproducibility
cv = cvpartition(EF_class, 'HoldOut', 0.3); % uses 70% to train, 30% to test

ParTrain = ANorm_Optp(cv.training, :);
EFTrain = EF_class(cv.training);
ParTest = ANorm_Optp(cv.test, :);
EFTest = EF_class(cv.test);

% Train the KNN Model
k = 5; % number of neighbors
knnModel = fitcknn(ParTrain, EFTrain,...
    'NumNeighbors', k, ...
    'Distance', 'euclidean', ...
    'Standardize', false); % already standarized/normalized above

% Predict on the test set and evaluate
EFPred = predict(knnModel, ParTest);
accuracy = sum(EFPred == EFTest)/numel(EFTest);
fprintf('Test Accuracy: %.2f%%\n', accuracy*100);

% Construct confusion matrix
figure(1);
c = confusionchart(EFTest, EFPred);
title('KNN Classification: Confusion Matrix');
set(c, 'FontSize', 18);

% Tune k using cross-validation
kValues = 1:2:21;
cvAccuracy = zeros(size(kValues));
for i = 1:length(kValues)
    knnCV = fitcknn(ParTrain, EFTrain, 'NumNeighbors', kValues(i));
    cvLoss = kfoldLoss(crossval(knnCV, 'KFold', 5));
    cvAccuracy(i) = 1 - cvLoss;
end

figure(2);
plot(kValues, cvAccuracy, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Neighbors (k)');
ylabel('Cross-Validated Accuracy');
title('Choosing Optimal k for KNN');
grid on;
set(gca, 'FontSize', 18);

%% LOGISTIC REGESSION
% Split into training and test sets
rng(1); % for reproducibility and to reseed from before
cv_log = cvpartition(size(ANorm_Optp, 1), 'HoldOut', 0.3);
ParTrain_log = ANorm_Optp(cv_log.training, :);
EFTrain_log = EF_class(cv_log.training, :);
ParTest_log = ANorm_Optp(cv_log.test, :);
EFTest_log = EF_class(cv_log.test, :);

% Fit logistic regression model
mdl = fitglm(ParTrain_log, EFTrain_log, 'Distribution', 'binomial', 'Link', 'logit');

% Predict probablities on test set and convert to class predictions
EFProb = predict(mdl, ParTest_log);
EFPred_log = double(EFProb >= 0.5);
catNames = categories(EF_class);
EFPred_log = categorical(EFPred_log, [0 1], catNames);

% Evaluate accuracy and create confusion matrix
accuracy_log = sum(EFPred_log == EFTest_log)/length(EFTest_log);
fprintf('Test Accuracy: %.2f%%\n', accuracy_log*100);

figure(3);
b = confusionchart(EFTest_log, EFPred_log);
title('Logistic Regression Classification: Confusion Matrix');
set(b, 'FontSize', 18)

