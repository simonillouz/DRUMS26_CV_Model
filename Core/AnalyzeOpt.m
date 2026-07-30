function AnalyzeOpt(PIDs)
% Provides statistics on the outputs of the optimization runs
% Identifies outliers, assesses boundary conditions and convergence
% Optional input: PIDs vector
% If no input, analyzes all folders in Optimization_Matrices

arguments
    PIDs (1,:) double = []   % empty = analyze all patients
end

baseFolder = 'MultiOptPatients';

% If no PIDs given, discover all patient folders
if isempty(PIDs)
    folders = dir(fullfile(baseFolder, 'Patient_*'));
    folders = folders([folders.isdir]);
    PIDs = arrayfun(@(f) sscanf(f.name, 'Patient_%d'), folders);
end

for k = 1:length(PIDs)
    PID       = PIDs(k);
    folderName = fullfile(baseFolder, sprintf('Patient_%d', PID));
    files      = dir(fullfile(folderName, '*.mat'));

    if isempty(files)
        fprintf('No results found for Patient %d, skipping.\n', PID);
        continue
    end

    % Load all runs
    allCosts = zeros(1, length(files));
    allPars  = [];
    for i = 1:length(files)
        d = load(fullfile(folderName, files(i).name));
        allCosts(i)   = d.cost;
        allPars(:, i) = d.optPars;
    end

    % Best run
    [bestCost, bestIdx] = min(allCosts);

    % Parameter statistics across runs
    meanPars = mean(allPars, 2);
    stdPars  = std(allPars,  0, 2);
    cvPars   = stdPars ./ abs(meanPars);  % coefficient of variation

    fprintf('\n=== Patient %d ===\n', PID);
    fprintf('Best run: %d  (cost = %.4f)\n', bestIdx, bestCost);
    fprintf('Cost range: [%.4f, %.4f]\n', min(allCosts), max(allCosts));
    fprintf('Parameters with highest variance (CV > 0.1):\n');
    highVar = find(cvPars > 0.1);
    for j = 1:length(highVar)
        fprintf('  Par %d: mean=%.4f  std=%.4f  CV=%.4f\n', ...
                highVar(j), meanPars(highVar(j)), stdPars(highVar(j)), cvPars(highVar(j)));
    end
end
end