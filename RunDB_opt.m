 function RunDB_opt(PIDs)

badset = [3, 111, 178, 228];
failedPIDs = [];

for k = 1:length(PIDs)
    if ~(ismember(k, badset))
        PID = PIDs(k);
        try
            sprintf(' -------Optimization for Patient %d-------',PID)
            DriverBasic_opt(PID,0); 
            sprintf(' -------Patient %d: MRI Volume-Done-------',PID)
            DriverBasic_opt(PID,1);
            sprintf(' ----Patient %d: Teichholz Volume-Done----',PID)
            % DriverBasic_opt(PID, 1); % for teichholz
            % RunDriverBasic(PID);
        catch ME
            fprintf('Patient %d failed. Error: %s\n', PID, ME.message);
            failedPIDs(end+1) = PID; %#ok<AGROW>
            continue;
        end
    end
    % cd ..
end
sprintf(' ----Done With %d MultiStart Optimizations----',length(setdiff(PIDs, badset)) - length(failedPIDs))

if ~isempty(failedPIDs)
    fprintf('The following %d PID(s) failed and were skipped:\n', length(failedPIDs));
    disp(failedPIDs);
else
    fprintf('All PIDs completed successfully.\n');
end

% Write failed PIDs to a text file
failedFile = fullfile(pwd, 'failedPIDs.txt');
fid = fopen(failedFile, 'w');
if fid == -1
    warning('Could not open %s for writing.', failedFile);
else
    for i = 1:length(failedPIDs)
        fprintf(fid, '%d\n', failedPIDs(i));
    end
    fclose(fid);
    fprintf('Failed PIDs written to %s\n', failedFile);
end

end