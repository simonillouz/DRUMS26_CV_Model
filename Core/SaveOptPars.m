function SaveOptPars(PIDs, teich)

BatchTableOpt = {};
BatchTableOptP = {};
validPIDs = [];
meanOpt = [];
stdOpt = [];

close all;
cd ../
if teich == 0
    cd(sprintf('MultiOptPatients'));
elseif teich == 1
    % cd(sprintf('MultiOptPatientsTeich'));
    cd(sprintf("MultiOptPatientsCalc")); % so that we have a different file name
end

badset = [];
for i = 1:length(PIDs)
    patientDir = sprintf('Patient_%d', i);

    if exist(patientDir, 'dir')
        cd(patientDir)
        NumRuns = numel(dir('*.mat'));
        if NumRuns == 0
            badset(end+1) = i;
            warning('Optimization matrices do not exist for Patient %d', i)
            cd ..
            continue;
        end
       
        for k = 1:NumRuns
            runFile = sprintf('Opt_Run%d.mat', k-1);
            if ~exist(runFile, 'file')
                badset(end+1) = i;
                warning('Opt_Run%d.mat does not exist for Patient %d', k, i)
                continue;
            end
        end
        cd ..
    else
        badset(end+1) = i;
        warning('Optimization matrices do not exist for Patient %d', i)
        continue;
    end
end
badset = unique(badset);

for k = 1:length(PIDs)
    PID = PIDs(k);
    if ~(ismember(PID, badset))
        parsOpt = [];
        cd ..
        if teich == 0
            for i = 0:15
                folderName = 'MultiOptPatients';
                folderName1 = sprintf('Patient_%d', PID);
                fileName = sprintf('Opt_Run%d.mat', i);
                filePath = fullfile(folderName, folderName1, fileName);
        
                if exist(filePath, 'file')
                    load(filePath);
                    parsOpt(:, i+1) = exp(optPars);
                else
                    error('Optimized data not found for Patient %d at: %s', PID, filePath);
                end
            end
        elseif teich == 1
            for i = 0:15
                % folderName = 'MultiOptPatientsTeich';
                folderName = 'MultiOptPatientsCalc';
                folderName1 = sprintf('Patient_%d', PID);
                fileName = sprintf('Opt_Run%d.mat', i);
                filePath = fullfile(folderName, folderName1, fileName);
                
                if exist(filePath, 'file')
                    load(filePath);
                    parsOpt(:, i+1) = exp(optPars);
                else
                    error('Optimized data not found for Patient %d at: %s', PID, filePath);
                end
            end
        end
        cd Core/
        nomVals = parsOpt(:,1);
        meanOpt = mean(parsOpt(:,2:16), 2);
        stdOpt  = std(parsOpt(:, 2:16), 0, 2);
        
        % List of optimized parameters pulled from covariances
        meanOptP = meanOpt([1 4 8:14],:); % contains only optimized parameters
        stdOptP = stdOpt([1 4 8:14], :); % contains only optimized parameters

        patientTableOpt = table(nomVals, meanOpt, stdOpt, ...
            'VariableNames', {'Nom Start', 'Mean Opt', 'Std Opt'});

        patientTableOptP = table(meanOptP, stdOptP,'VariableNames', {'Mean Opt', 'Std Opt'});

        BatchTableOpt{end+1} = patientTableOpt;
        BatchTableOptP{end+1} = patientTableOptP;
        validPIDs(end+1) = PID;
    end
end

if teich == 0
    filename = 'Multistart Patient Optimized Data.xlsx';
elseif teich == 1
    % filename = 'Multistart Teichholz Patient Optimized Data.xlsx';
    filename = 'Multistart Calculated Patient Optimized Data.xlsx';
end
if exist(filename, 'file')
    delete(filename);
end

nParams = size(parsOpt, 1);
paramNums = strcat("Param", string(1:nParams))';
paramNames = {'Rpul', 'Rmval', 'Raval', 'Rsys', 'Rtval', 'Rpval', 'Cpv', 'Csa', 'Csv', 'Cpa', 'ElvM', 'Elvm', 'ErvM', 'Ervm', 'TC', 'TR'}';

paramNumsP = strcat("Param", string([1 4 8:14]))';
paramNamesP = {'Rpul', 'Rsys', 'Csa', 'Csv', 'Cpa', 'ElvM', 'Elvm', 'ErvM', 'Ervm'}';

if teich == 0
    paramTable = table(paramNums, paramNames, 'VariableNames', {'Parameter Num', 'Parameter Name'});
    writetable(paramTable, filename, 'Sheet', 'AllPatients', 'Range', 'A2:B18')
    
    paramTableP = table(paramNumsP, paramNamesP, 'VariableNames', {'Parameter Num', 'Parameter Name'});
    writetable(paramTableP, filename, 'Sheet', 'Optimized Parameters Only', 'Range', 'A2:B11')

    for k = 1:numel(BatchTableOpt)
        startCol = (k-1)*3+3;
        colLetter = num2colLetter(startCol);
    
        writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
            'Sheet', 'AllPatients', 'Range', sprintf('%s1', colLetter));
    
        writetable(BatchTableOpt{k}, filename, ...
            'Sheet', 'AllPatients', 'Range', sprintf('%s2', colLetter));
    end

    for k = 1:numel(BatchTableOptP)
        startColP = (k-1)*2+3;
        colLetterP = num2colLetter(startColP);
    
        writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
            'Sheet', 'Optimized Parameters Only', 'Range', sprintf('%s1', colLetterP));
    
        writetable(BatchTableOptP{k}, filename, ...
            'Sheet', 'Optimized Parameters Only', 'Range', sprintf('%s2', colLetterP));
    end
elseif teich == 1
    paramTable = table(paramNums, paramNames, 'VariableNames', {'Parameter Num', 'Parameter Name'});
    % writetable(paramTable, filename, 'Sheet', 'AllPatientsTeich', 'Range', 'A2:B18')
    writetable(paramTable, filename, 'Sheet', 'AllPatientsCalc', 'Range', 'A2:B18')

    paramTableP = table(paramNumsP, paramNamesP, 'VariableNames', {'Parameter Num', 'Parameter Name'});
    % writetable(paramTableP, filename, 'Sheet', 'Teichholz Opt Pars Only', 'Range', 'A2:B11')
    writetable(paramTableP, filename, 'Sheet', 'Calculated Opt Pars Only', 'Range', 'A2:B11')


    for k = 1:numel(BatchTableOpt)
        startCol = (k-1)*3+3;
        colLetter = num2colLetter(startCol);
    
        % writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
        %     'Sheet', 'AllPatientsTeich', 'Range', sprintf('%s1', colLetter));
        % 
        % writetable(BatchTableOpt{k}, filename, ...
        %     'Sheet', 'AllPatientsTeich', 'Range', sprintf('%s2', colLetter));

        writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
            'Sheet', 'AllPatientsCalc', 'Range', sprintf('%s1', colLetter));
    
        writetable(BatchTableOpt{k}, filename, ...
            'Sheet', 'AllPatientsCalc', 'Range', sprintf('%s2', colLetter));
    end

    for k = 1:numel(BatchTableOptP)
        startColP = (k-1)*2+3;
        colLetterP = num2colLetter(startColP);
    
        % writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
        %     'Sheet', 'Teichholz Opt Pars Only', 'Range', sprintf('%s1', colLetterP));
        % 
        % writetable(BatchTableOptP{k}, filename, ...
        %     'Sheet', 'Teichholz Opt Pars Only', 'Range', sprintf('%s2', colLetterP));

        writecell({sprintf('Patient %d', validPIDs(k))}, filename, ...
            'Sheet', 'Calculated Opt Pars Only', 'Range', sprintf('%s1', colLetterP));
    
        writetable(BatchTableOptP{k}, filename, ...
            'Sheet', 'Calculated Opt Pars Only', 'Range', sprintf('%s2', colLetterP));
    end
end
end

function col = num2colLetter(n)
    col = '';
    while n > 0
        rem = mod(n-1, 26);
        col = [char(65 + rem), col];
        n = floor((n-1)/26);
    end
end
