function ParsBox(PIDs)
% This code makes a box and whisker plots for each parameter, putting
% calculated and measured volume parameters side by side
close all;

% First extract optimized parameters from the MRI Volumes
filename = 'Multistart Patient Optimized Data.xlsx';
sheet = 'Optimized Parameters Only';
% headerRow = readcell(filename, 'Sheet', sheet, 'Range', 'C1:YR1');
headerRow = readcell(filename, 'Sheet', sheet, 'Range', 'C1:YR1');
patientCol = NaN(1, numel(headerRow));
currentPID = NaN;

for i = 1:numel(headerRow)
    val = headerRow{i};
    isEmptyCell = ( isa(val, 'missing') || (ischar(val) && isempty(val)));
    if ~isEmptyCell
        currentPID = str2double(regexp(char(val), '\d+', 'match', 'once'));
    end
    patientCol(i) = currentPID;
end

paramNames = readcell(filename, 'Sheet', sheet, 'Range', 'B3:B11');
dataBlock = readmatrix(filename, 'Sheet', sheet, 'Range', 'C3:YR11');

meanCols = 1:2:size(dataBlock, 2);
meanData = dataBlock(:, meanCols);
patAvail = patientCol(meanCols);

nParams = size(meanData, 1);
selData_MRI = NaN(numel(PIDs), nParams);

for k = 1:numel(PIDs)
    colIdx = find(patAvail == PIDs(k), 1);
    if isempty(colIdx)
        % warning('PID %d not found in MRI data.', PIDs(k))
        continue;
    end
    
    selData_MRI(k,:) = meanData(:, colIdx)';
end

% Then extract optimized parameters from the TTE Volumes
filename = 'Multistart Teichholz Patient Optimized Data.xlsx';
sheet = 'Teichholz Opt Pars Only';
headerRow = readcell(filename, 'Sheet', sheet, 'Range', 'C1:YN1');
patientCol = NaN(1, numel(headerRow));
currentPID = NaN;

for i = 1:numel(headerRow)
    val = headerRow{i};
    isEmptyCell = ( isa(val, 'missing') || (ischar(val) && isempty(val)));
    if ~isEmptyCell
        currentPID = str2double(regexp(char(val), '\d+', 'match', 'once'));
    end
    patientCol(i) = currentPID;
end

%  paramNames_TTE = readcell(filename, 'Sheet', sheet, 'Range', 'B3:B11');
dataBlock = readmatrix(filename, 'Sheet', sheet, 'Range', 'C3:YN11');

meanCols = 1:2:size(dataBlock, 2);
meanData = dataBlock(:, meanCols);
patAvail = patientCol(meanCols);

nParams = size(meanData, 1);
selData_TTE = NaN(numel(PIDs), nParams);

for k = 1:numel(PIDs)
    colIdx = find(patAvail == PIDs(k), 1);
    if isempty(colIdx)
        % warning('PID %d not found in MRI data.', PIDs(k))
        continue;
    end
    
    selData_TTE(k,:) = meanData(:, colIdx)';
end

% Graph box and whisker plots for parameters
warning('off', 'stats:boxplot:BadObjectType');

figure(1);
t = tiledlayout('flow');
title(t, 'Parameter Values Across Patients');
% For comparison of optimized parameters from both MRI and TTE measurements
for c = 1:numel(paramNames)
    nexttile;
    combinedData = [selData_MRI(:,c); selData_TTE(:,c)];
    groupLabels = [repmat({'MRI'}, numel(selData_MRI(:,c)), 1); repmat({'TTE'}, numel(selData_TTE(:,c)), 1)];
    b = boxplot(combinedData, groupLabels, 'Widths', 0.6);
    set(b, 'LineWidth', 2);
    xlabel(paramNames(c), 'FontSize', 20);
    grid on;
end

% For only parameters optimized from MRI measurements
% for c = 1:numel(paramNames)
%     nexttile;
%     b = boxplot(selData_MRI(:,c), 'Widths', 0.6);
%     set(b, 'LineWidth', 2);
%     xlabel(paramNames(c), 'FontSize', 20);
%     grid on;
% end

warning('on', 'stats:boxplot:BadObjectType');

% Check for statistical differences
% First perform t-test
disp('Results from t-test for each parameter:')
for i = 1:numel(paramNames)
    [h, p, ci, stats] = ttest2(selData_MRI(:, i), selData_TTE(:,i));
    fprintf('The p-value for %s is %.4f. The h-value is %d.\n', paramNames{i}, p, h);
end

% Then perform KL Divergence
disp('Results from KL Divergence for each parameter:')
for i = 1:numel(paramNames)
    MRI_data = selData_MRI(:, i);
    TTE_data = selData_TTE(:, i);
    combined = [MRI_data; TTE_data];

    % Using Freedman-Diaconis bin width
    binWidth = 2*iqr(combined)*numel(combined)^(-1/3);
    nbins = max(ceil(range(combined)/binWidth), 5); % floor of 5 bins as a sanity check

    edges = linspace(min(combined), max(combined), nbins+1);

    p = histcounts(MRI_data, edges, 'Normalization', 'probability') + eps;
    q = histcounts(TTE_data, edges, 'Normalization', 'probability') + eps;
    p = p/sum(p);
    q = q/sum(q);

    KL_PQ = sum(p.*log(p./q));
    KL_QP = sum(q.*log(q./p));
    fprintf('%s: KL(P||Q) = %.4f, KL(Q||P) = %.4f, using %d bins\n', paramNames{i}, KL_PQ, KL_QP, nbins);

    % Plots
    figure('Name', paramNames{i}, 'Position', [100, 100, 1400, 500]);
    xGrid = linspace(min(combined), max(combined), 500);
    fP = ksdensity(MRI_data, xGrid);
    fQ = ksdensity(TTE_data, xGrid);

    % Subplot 1: KL Divergence Plot
    subplot(1,2,1);
    cla; hold on;
    upperP = max(fP, fQ);
    lowerP = min(fP, fQ);
    fill([xGrid, fliplr(xGrid)], [fP, fliplr(fQ)], ...
        [1 0.6 0.6], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    plot(xGrid, fP, 'r', 'LineWidth', 2, 'DisplayName', 'P(x)');
    plot(xGrid, fQ, 'b', 'LineWidth', 2, 'DisplayName', 'Q(x)');
    xlabel('x');
    ylabel('Probability Density');
    title(sprintf('%s: KL(P||Q)=%.4f, KL(Q||P)=%.4f', paramNames{i}, KL_PQ, KL_QP));
    legend('Divergence region', 'P(x)', 'Q(x)', 'Location', 'best');
    hold off;

    % Subplot 2: Histogram
    subplot(1,2,2);
    cla; hold on;
    histogram(MRI_data, edges, 'Normalization', 'probability', ...
        'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'MRI Pars');
    histogram(TTE_data, edges, 'Normalization', 'probability', ...
        'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'TTE Pars');
    xlabel(paramNames{i});
    ylabel('Probability');
    title(sprintf('%s: Histogram Comparison', paramNames{i}));
    legend('Location', 'best');
    hold off;
end


end