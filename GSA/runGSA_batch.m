%% Batch Global Sensitivity Analysis (GSA) Execution Function
% This function automates the execution of the GSA (run_cardio_GSA) across
% multiple patients. It supports running a specific list of patient IDs or
% randomly selecting a sebset of patients from the cohort
function runGSA_batch(PatientBatch, N_samples, AnalysisType, Target)

%% 1. Determine Cohort & Select Patients to Run
% Define the full cohort range (total number of available patients in the
% dataset)
allPatients = 1:343;
% If the input PatientBatch is a single scalar number instead of an array,
% interpret it as the requested count of patients to select randomly
if isscalar(PatientBatch)
    numRandom = PatientBatch;
    % Ensure the requested random count does not eceed the total cohort
    % size
    if numRandom > length(allPatients)
        error('Requested %d patients, but only %d exist.', ...
            numRandom, length(allPatients));
    end
    % Randomly sample patient IDs without replacement, then sort them in
    % ascending order
    PatientBatch = sort(randperm(length(allPatients), numRandom));
    fprintf('Randomly selected patients:\n');
    disp(PatientBatch)
end

%% 2. Display Batch Configuration Details
fprintf('Starting Batch Global Sensitivity Analysis... \n');
fprintf('Analysis Type: %s\n', AnalysisType);
fprintf('Patients: %s\n', mat2str(PatientBatch));
fprintf('Target: %s\n', Target);
fprintf('Samples: %d\n', N_samples);

numPatients = length(PatientBatch);

% Start the global timer to track the total execution time of the entire
% batch run
batchStart = tic;

%% 3. Loop Through Patient Batch
for k = 1:numPatients
    PID = PatientBatch(k);
    fprintf('Patient %d (%d of %d)\n', PID, k, numPatients);
    % Start a local timer to track the evaluation time for this specific
    % patient
    patientStart = tic;
    try
        % Execute the underlying GSA function for the current patient
        run_cardio_GSA(PID, N_samples, AnalysisType, Target);

        fprintf('Patient %d completed in %.2f seconds.\n', ...
            PID, toc(patientStart));
    catch ME
        % Print an error to prevent a single patient's solver failure from
        % terminating the entire batch queue
        fprintf('Patient %d failed:\n', PID);
        fprintf('%s\n', ME.message);
    end
end

% Print the total execution duration in minutes
fprintf('Batch completed in %.2f minutes.\n', ...
    toc(batchStart)/60);
end